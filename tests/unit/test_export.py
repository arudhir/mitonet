"""
Unit tests for network export and filtering functionality
"""

import pytest
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch
import pandas as pd
import networkx as nx

from mitonet.export import NetworkFilter, NetworkExporter, export_predefined_networks
from mitonet.database import MitoNetDatabase


@pytest.mark.database
class TestNetworkFilter:
    """Test NetworkFilter class and its factory methods"""
    
    def test_empty_filter(self):
        """Test empty filter creation"""
        filter_obj = NetworkFilter()
        assert filter_obj.include_mitochondrial is None
        assert filter_obj.include_muscle_expressed is None
        assert filter_obj.gene_symbols is None
        assert filter_obj.min_confidence is None
        assert filter_obj.include_neighbors == 0
    
    def test_mitochondrial_filter(self):
        """Test mitochondrial filter factory"""
        filter_obj = NetworkFilter.mitochondrial_network(min_confidence=0.5)
        assert filter_obj.include_mitochondrial is True
        assert filter_obj.min_confidence == 0.5
        assert filter_obj.include_muscle_expressed is None
    
    def test_muscle_filter(self):
        """Test muscle filter factory"""
        filter_obj = NetworkFilter.muscle_network(min_confidence=0.6)
        assert filter_obj.include_muscle_expressed is True
        assert filter_obj.min_confidence == 0.6
        assert filter_obj.include_mitochondrial is None
    
    def test_gene_set_filter(self):
        """Test gene set filter factory"""
        genes = ['ATP1A1', 'MYOD1', 'CYC1']
        filter_obj = NetworkFilter.gene_set_network(genes, include_neighbors=2, min_confidence=0.4)
        assert filter_obj.gene_symbols == set(genes)
        assert filter_obj.include_neighbors == 2
        assert filter_obj.min_confidence == 0.4
    
    def test_high_confidence_filter(self):
        """Test high confidence filter factory"""
        filter_obj = NetworkFilter.high_confidence_network(min_confidence=0.8)
        assert filter_obj.min_confidence == 0.8
        assert filter_obj.include_mitochondrial is None


@pytest.mark.database
class TestNetworkExporter:
    """Test NetworkExporter class"""
    
    @pytest.fixture
    def test_exporter(self, temp_db, tmp_path):
        """Create test exporter with temporary database and output directory"""
        return NetworkExporter(temp_db, tmp_path / "outputs")
    
    @pytest.fixture
    def populated_export_db(self, temp_db):
        """Create a database with sample data for export testing"""
        # Create proteins
        protein1 = temp_db.get_or_create_protein(
            uniprot_id="P12345",
            gene_symbol="ATP1A1",
            is_mitochondrial=True,
            is_muscle_expressed=True,
            muscle_tpm=45.6
        )
        
        protein2 = temp_db.get_or_create_protein(
            uniprot_id="Q67890",
            gene_symbol="MYOD1",
            is_mitochondrial=False,
            is_muscle_expressed=True,
            muscle_tpm=123.4
        )
        
        protein3 = temp_db.get_or_create_protein(
            uniprot_id="P00123",
            gene_symbol="CYC1",
            is_mitochondrial=True,
            is_muscle_expressed=False,
            muscle_tpm=0.0
        )
        
        # Create data source
        source = temp_db.get_or_create_data_source(
            name="TEST_SOURCE",
            version="1.0",
            file_path="/test/path"
        )
        
        # Add interactions
        temp_db.add_interaction(
            protein1=protein1,
            protein2=protein2,
            source=source,
            confidence_score=0.8,
            evidence_type="experimental"
        )
        
        temp_db.add_interaction(
            protein1=protein1,
            protein2=protein3,
            source=source,
            confidence_score=0.6,
            evidence_type="database"
        )
        
        temp_db.add_interaction(
            protein1=protein2,
            protein2=protein3,
            source=source,
            confidence_score=0.4,
            evidence_type="textmining"
        )
        
        return temp_db
    
    def test_get_filtered_proteins_all(self, test_exporter, populated_export_db):
        """Test getting all proteins with no filter"""
        filter_obj = NetworkFilter()
        proteins = test_exporter._get_filtered_proteins(filter_obj)
        assert len(proteins) == 3
        uniprot_ids = {p.uniprot_id for p in proteins}
        assert uniprot_ids == {"P12345", "Q67890", "P00123"}
    
    def test_get_filtered_proteins_mitochondrial(self, test_exporter, populated_export_db):
        """Test filtering for mitochondrial proteins only"""
        filter_obj = NetworkFilter()
        filter_obj.include_mitochondrial = True
        proteins = test_exporter._get_filtered_proteins(filter_obj)
        assert len(proteins) == 2
        uniprot_ids = {p.uniprot_id for p in proteins}
        assert uniprot_ids == {"P12345", "P00123"}
    
    def test_get_filtered_proteins_muscle(self, test_exporter, populated_export_db):
        """Test filtering for muscle-expressed proteins only"""
        filter_obj = NetworkFilter()
        filter_obj.include_muscle_expressed = True
        proteins = test_exporter._get_filtered_proteins(filter_obj)
        assert len(proteins) == 2
        uniprot_ids = {p.uniprot_id for p in proteins}
        assert uniprot_ids == {"P12345", "Q67890"}
    
    def test_get_filtered_proteins_gene_symbols(self, test_exporter, populated_export_db):
        """Test filtering by specific gene symbols"""
        filter_obj = NetworkFilter()
        filter_obj.gene_symbols = {"ATP1A1", "CYC1"}
        proteins = test_exporter._get_filtered_proteins(filter_obj)
        assert len(proteins) == 2
        gene_symbols = {p.gene_symbol for p in proteins}
        assert gene_symbols == {"ATP1A1", "CYC1"}
    
    def test_get_filtered_interactions(self, test_exporter, populated_export_db):
        """Test filtering interactions by confidence"""
        filter_obj = NetworkFilter()
        
        # Get all proteins first
        proteins = test_exporter._get_filtered_proteins(filter_obj)
        
        # Test with minimum confidence filter
        filter_obj.min_confidence = 0.5
        interactions = test_exporter._get_filtered_interactions(filter_obj, proteins)
        assert len(interactions) == 2  # 0.8 and 0.6 confidence interactions
        
        # Test with higher confidence filter
        filter_obj.min_confidence = 0.7
        interactions = test_exporter._get_filtered_interactions(filter_obj, proteins)
        assert len(interactions) == 1  # Only 0.8 confidence interaction
    
    def test_apply_degree_filters(self, test_exporter, populated_export_db):
        """Test applying degree filters"""
        filter_obj = NetworkFilter()
        proteins = test_exporter._get_filtered_proteins(filter_obj)
        interactions = test_exporter._get_filtered_interactions(filter_obj, proteins)
        
        # Test minimum degree filter
        filter_obj.min_degree = 2
        filtered_proteins, filtered_interactions = test_exporter._apply_degree_filters(
            proteins, interactions, filter_obj
        )
        
        # All proteins have degree >= 2 in our test data
        assert len(filtered_proteins) == 3
        
        # Test maximum degree filter
        filter_obj.min_degree = None
        filter_obj.max_degree = 1
        filtered_proteins, filtered_interactions = test_exporter._apply_degree_filters(
            proteins, interactions, filter_obj
        )
        
        # No proteins should have degree <= 1 in our test data
        assert len(filtered_proteins) == 0
    
    def test_build_networkx_graph(self, test_exporter, populated_export_db):
        """Test building NetworkX graph from proteins and interactions"""
        filter_obj = NetworkFilter()
        proteins = test_exporter._get_filtered_proteins(filter_obj)
        interactions = test_exporter._get_filtered_interactions(filter_obj, proteins)
        
        graph = test_exporter._build_networkx_graph(proteins, interactions)
        
        # Check graph structure
        assert len(graph.nodes) == 3
        assert len(graph.edges) == 3
        
        # Check node attributes
        node_data = graph.nodes(data=True)
        p12345_data = graph.nodes["P12345"]
        assert p12345_data["gene_symbol"] == "ATP1A1"
        assert p12345_data["is_mitochondrial"] is True
        assert p12345_data["is_muscle_expressed"] is True
        assert p12345_data["muscle_tpm"] == 45.6
        
        # Check edge attributes
        edge_data = graph.edges(data=True)
        assert len(edge_data) == 3
        
        # Find specific edge
        for u, v, data in edge_data:
            if {u, v} == {"P12345", "Q67890"}:
                assert data["confidence_score"] == 0.8
                assert data["evidence_type"] == "experimental"
                break
        else:
            pytest.fail("Expected edge P12345-Q67890 not found")
    
    @patch('networkx.write_graphml')
    def test_export_graphml(self, mock_write_graphml, test_exporter):
        """Test GraphML export"""
        # Create a simple graph for testing
        graph = nx.Graph()
        graph.add_node("P12345", gene_symbol="ATP1A1")
        graph.add_edge("P12345", "Q67890", confidence_score=0.8)
        
        output_file = test_exporter._export_graphml(graph, "test_network")
        
        # Check that NetworkX write function was called
        mock_write_graphml.assert_called_once()
        
        # Check output file path
        assert output_file.name == "test_network.graphml"
        assert output_file.parent == test_exporter.output_dir
    
    def test_export_json(self, test_exporter):
        """Test JSON export"""
        # Create a simple graph for testing
        graph = nx.Graph()
        graph.add_node("P12345", gene_symbol="ATP1A1")
        graph.add_edge("P12345", "Q67890", confidence_score=0.8)
        
        output_file = test_exporter._export_json(graph, "test_network")
        
        # Check output file path
        assert output_file.name == "test_network.json"
        assert output_file.exists()
        
        # Check that file contains valid JSON
        import json
        with open(output_file) as f:
            data = json.load(f)
        
        assert "nodes" in data
        assert "links" in data
        assert len(data["nodes"]) == 2
        assert len(data["links"]) == 1
    
    def test_export_csv(self, test_exporter, populated_export_db):
        """Test CSV export"""
        filter_obj = NetworkFilter()
        proteins = test_exporter._get_filtered_proteins(filter_obj)
        interactions = test_exporter._get_filtered_interactions(filter_obj, proteins)
        
        csv_files = test_exporter._export_csv(proteins, interactions, "test_network")
        
        # Check files were created
        assert "nodes_csv" in csv_files
        assert "edges_csv" in csv_files
        
        nodes_file = csv_files["nodes_csv"]
        edges_file = csv_files["edges_csv"]
        
        assert nodes_file.exists()
        assert edges_file.exists()
        
        # Check CSV content
        nodes_df = pd.read_csv(nodes_file)
        edges_df = pd.read_csv(edges_file)
        
        assert len(nodes_df) == 3
        assert len(edges_df) == 3
        
        # Check node columns
        expected_node_cols = ['uniprot_id', 'gene_symbol', 'is_mitochondrial', 'is_muscle_expressed']
        for col in expected_node_cols:
            assert col in nodes_df.columns
        
        # Check edge columns
        expected_edge_cols = ['protein1', 'protein2', 'confidence_score', 'evidence_type']
        for col in expected_edge_cols:
            assert col in edges_df.columns


@pytest.mark.database
class TestExportIntegration:
    """Integration tests for complete export workflow"""
    
    def test_export_network_complete_workflow(self, temp_db, tmp_path):
        """Test complete network export workflow"""
        # Setup database with test data
        protein1 = temp_db.get_or_create_protein(
            uniprot_id="P12345",
            gene_symbol="ATP1A1",
            is_mitochondrial=True
        )
        
        protein2 = temp_db.get_or_create_protein(
            uniprot_id="Q67890",
            gene_symbol="MYOD1",
            is_muscle_expressed=True
        )
        
        source = temp_db.get_or_create_data_source("TEST", "1.0", "/test")
        temp_db.add_interaction(protein1, protein2, source, 0.8)
        
        # Test export
        exporter = NetworkExporter(temp_db, tmp_path / "outputs")
        filter_obj = NetworkFilter.mitochondrial_network()
        
        output_files = exporter.export_network(
            network_filter=filter_obj,
            format_types=['json', 'csv'],
            filename_prefix='mito_test'
        )
        
        # Check outputs
        assert 'json' in output_files
        assert 'nodes_csv' in output_files
        assert 'edges_csv' in output_files
        
        # Verify files exist
        for file_path in output_files.values():
            assert file_path.exists()
    
    def test_export_predefined_networks(self, temp_db, tmp_path):
        """Test export of predefined network types"""
        # Setup database with diverse test data
        proteins = []
        
        # Create mitochondrial protein
        proteins.append(temp_db.get_or_create_protein(
            uniprot_id="P12345",
            gene_symbol="ATP1A1",
            is_mitochondrial=True,
            is_muscle_expressed=False
        ))
        
        # Create muscle protein
        proteins.append(temp_db.get_or_create_protein(
            uniprot_id="Q67890",
            gene_symbol="MYOD1",
            is_mitochondrial=False,
            is_muscle_expressed=True
        ))
        
        # Create both mitochondrial and muscle protein
        proteins.append(temp_db.get_or_create_protein(
            uniprot_id="P00123",
            gene_symbol="CYC1",
            is_mitochondrial=True,
            is_muscle_expressed=True
        ))
        
        # Add interactions
        source = temp_db.get_or_create_data_source("TEST", "1.0", "/test")
        temp_db.add_interaction(proteins[0], proteins[1], source, 0.8)
        temp_db.add_interaction(proteins[1], proteins[2], source, 0.6)
        temp_db.add_interaction(proteins[0], proteins[2], source, 0.9)
        
        # Test predefined export
        networks_exported = export_predefined_networks(temp_db, tmp_path / "outputs")
        
        # Check that all predefined networks were created
        expected_networks = ['mitochondrial', 'muscle', 'high_confidence']
        for network_type in expected_networks:
            assert network_type in networks_exported
            
            # Check files for each network
            files = networks_exported[network_type]
            assert 'json' in files
            assert 'graphml' in files
            assert 'nodes_csv' in files
            assert 'edges_csv' in files
            
            # Verify files exist
            for file_path in files.values():
                assert file_path.exists()


@pytest.mark.database
class TestNetworkFilterValidation:
    """Test filter validation and edge cases"""
    
    def test_filter_with_no_proteins(self, temp_db, tmp_path):
        """Test filter that results in no proteins"""
        exporter = NetworkExporter(temp_db, tmp_path / "outputs")
        
        # Filter for mitochondrial proteins when database is empty
        filter_obj = NetworkFilter.mitochondrial_network()
        
        output_files = exporter.export_network(
            network_filter=filter_obj,
            format_types=['json'],
            filename_prefix='empty_test'
        )
        
        # Should still create files, but with empty networks
        assert 'json' in output_files
        assert output_files['json'].exists()
    
    def test_filter_complex_criteria(self, temp_db, tmp_path):
        """Test filter with multiple complex criteria"""
        # Setup test data
        protein1 = temp_db.get_or_create_protein(
            uniprot_id="P12345",
            gene_symbol="ATP1A1",
            is_mitochondrial=True,
            is_muscle_expressed=True
        )
        
        protein2 = temp_db.get_or_create_protein(
            uniprot_id="Q67890",
            gene_symbol="MYOD1",
            is_mitochondrial=False,
            is_muscle_expressed=True
        )
        
        source = temp_db.get_or_create_data_source("TEST", "1.0", "/test")
        temp_db.add_interaction(protein1, protein2, source, 0.8, evidence_type="experimental")
        
        # Test complex filter
        exporter = NetworkExporter(temp_db, tmp_path / "outputs")
        filter_obj = NetworkFilter()
        filter_obj.include_mitochondrial = True
        filter_obj.include_muscle_expressed = True
        filter_obj.min_confidence = 0.7
        filter_obj.evidence_types = {"experimental"}
        
        output_files = exporter.export_network(
            network_filter=filter_obj,
            format_types=['json'],
            filename_prefix='complex_test'
        )
        
        # Should create valid output
        assert 'json' in output_files
        assert output_files['json'].exists()