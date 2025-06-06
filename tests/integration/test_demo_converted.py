"""
Converted demo functionality as integration tests
"""

import pytest
from pathlib import Path
from unittest.mock import patch

from mitonet.database import MitoNetDatabase
from mitonet.ingestion import DataIngestionManager


@pytest.mark.integration
class TestDemoScenarios:
    """Test scenarios from the original demo script"""
    
    def test_demo_workflow_comprehensive(self, tmp_path):
        """Test the complete demo workflow as a comprehensive integration test"""
        # Initialize database (equivalent to demo step 1)
        db_path = tmp_path / "demo_test.db"
        db = MitoNetDatabase(str(db_path))
        db.initialize_db()
        
        # Check initial state
        stats = db.get_statistics()
        assert stats['num_proteins'] == 0
        assert stats['num_interactions'] == 0
        assert stats['num_sources'] == 0
        
        # Add genes manually (equivalent to demo step 2)
        protein1 = db.get_or_create_protein("P12345", gene_symbol="ATP1A1", is_muscle_expressed=True)
        protein2 = db.get_or_create_protein("Q67890", gene_symbol="MYOD1", is_muscle_expressed=True)
        protein3 = db.get_or_create_protein("P00123", gene_symbol="CYC1", is_mitochondrial=True)
        
        # Verify manual additions
        stats = db.get_statistics()
        assert stats['num_proteins'] == 3
        assert stats['num_mitochondrial'] == 1
        assert stats['num_muscle_expressed'] == 2
        
        # Setup data directory
        data_dir = tmp_path / "networks"
        self._create_test_data_files(data_dir)
        
        # Initialize ingestion manager (equivalent to demo step 3)
        ingestion = DataIngestionManager(db, data_dir)
        
        # Test change detection
        source_files = {
            'MitoCarta': data_dir / 'mitocarta/Human.MitoCarta3.0.xls',
            'HPA_muscle': data_dir / 'hpa/hpa_skm.tsv',
        }
        
        for source_name, file_path in source_files.items():
            if file_path.exists():
                needs_update = ingestion.needs_update(source_name, file_path)
                assert needs_update is True  # New files should need update
        
        # Test incremental updates (equivalent to demo step 4)
        # Add aliases first for proper linking
        db.add_protein_alias(protein1, 'symbol', 'ATP1A1', None)
        db.add_protein_alias(protein2, 'symbol', 'MYOD1', None)
        db.add_protein_alias(protein3, 'symbol', 'CYC1', None)
        
        # Update HPA data
        if source_files['HPA_muscle'].exists():
            hpa_source = ingestion.ingest_hpa_muscle(source_files['HPA_muscle'])
            assert hpa_source.name == "HPA_muscle"
            assert hpa_source.version is not None
        
        # Update MitoCarta data (with mocking since Excel is complex)
        with patch('pandas.read_excel') as mock_excel:
            mock_excel.return_value = self._get_mock_mitocarta_data()
            mitocarta_source = ingestion.ingest_mitocarta(source_files['MitoCarta'])
            assert mitocarta_source.name == "MitoCarta"
            assert mitocarta_source.version is not None
        
        # Verify final state (equivalent to demo step 5)
        final_stats = db.get_statistics()
        assert final_stats['num_proteins'] == 3
        assert final_stats['num_mitochondrial'] >= 1  # Should have mitochondrial proteins
        assert final_stats['num_muscle_expressed'] >= 1  # Should have muscle proteins
        assert final_stats['num_sources'] >= 2  # Should have both sources
        
        # Verify data sources are tracked
        assert len(final_stats['data_sources']) >= 2
        source_names = [source.split()[0] for source in final_stats['data_sources']]
        assert 'HPA_muscle' in source_names
        assert 'MitoCarta' in source_names
        
        # Test protein attribute updates
        updated_protein1 = db.get_protein_by_uniprot("P12345")
        updated_protein3 = db.get_protein_by_uniprot("P00123")
        
        # Verify proteins exist and have expected attributes
        assert updated_protein1 is not None
        assert updated_protein3 is not None
        assert updated_protein3.is_mitochondrial is True
    
    def test_demo_gene_addition_scenarios(self, tmp_path):
        """Test various gene addition scenarios from demo"""
        db = MitoNetDatabase(str(tmp_path / "gene_test.db"))
        db.initialize_db()
        
        # Test adding genes by UniProt ID
        protein_by_uniprot = db.get_or_create_protein("P12345", gene_symbol="TEST1")
        assert protein_by_uniprot.uniprot_id == "P12345"
        assert protein_by_uniprot.gene_symbol == "TEST1"
        
        # Test adding genes with attributes
        protein_with_attrs = db.get_or_create_protein(
            "Q67890",
            gene_symbol="TEST2",
            is_mitochondrial=True,
            muscle_tpm=50.0,
            priority_score=3.5
        )
        assert protein_with_attrs.is_mitochondrial is True
        assert protein_with_attrs.muscle_tpm == 50.0
        assert protein_with_attrs.priority_score == 3.5
        
        # Test updating existing protein
        updated_protein = db.get_or_create_protein(
            "P12345",  # Same UniProt ID
            gene_symbol="TEST1_UPDATED",
            is_muscle_expressed=True
        )
        assert updated_protein.id == protein_by_uniprot.id  # Same protein
        assert updated_protein.gene_symbol == "TEST1_UPDATED"  # Updated attribute
        assert updated_protein.is_muscle_expressed is True  # New attribute
        
        # Verify database consistency
        stats = db.get_statistics()
        assert stats['num_proteins'] == 2  # Only 2 unique proteins
    
    def test_demo_data_source_versioning(self, tmp_path):
        """Test data source versioning scenarios from demo"""
        db = MitoNetDatabase(str(tmp_path / "version_test.db"))
        db.initialize_db()
        
        # Create multiple versions of the same source
        source_v1 = db.get_or_create_data_source(
            name="DEMO_SOURCE",
            version="1.0",
            file_path="/demo/path/v1",
            file_hash="hash_v1"
        )
        
        source_v2 = db.get_or_create_data_source(
            name="DEMO_SOURCE", 
            version="2.0",
            file_path="/demo/path/v2",
            file_hash="hash_v2"
        )
        
        # Should be different sources
        assert source_v1.id != source_v2.id
        assert source_v1.version == "1.0"
        assert source_v2.version == "2.0"
        
        # Test getting existing source
        existing_source = db.get_or_create_data_source(
            name="DEMO_SOURCE",
            version="1.0", 
            file_path="/demo/path/v1"
        )
        assert existing_source.id == source_v1.id
        
        # Verify statistics
        stats = db.get_statistics()
        assert stats['num_sources'] == 2
        data_source_strings = stats['data_sources']
        assert 'DEMO_SOURCE v1.0' in data_source_strings
        assert 'DEMO_SOURCE v2.0' in data_source_strings
    
    def test_demo_checkpoint_functionality(self, tmp_path):
        """Test checkpoint functionality as shown in demo"""
        db = MitoNetDatabase(str(tmp_path / "checkpoint_test.db"))
        db.initialize_db()
        
        # Test saving various checkpoint scenarios
        checkpoint_scenarios = [
            {
                "name": "protein_ingestion_phase",
                "phase": "data_ingestion", 
                "data": {"proteins_processed": 1000, "aliases_added": 2500}
            },
            {
                "name": "interaction_building",
                "phase": "network_construction",
                "data": {"interactions_processed": 50000, "edges_added": 45000}
            },
            {
                "name": "quality_control",
                "phase": "validation",
                "data": {"nodes_validated": 3000, "edges_validated": 45000, "errors_found": 5}
            }
        ]
        
        # Save all checkpoints
        for checkpoint in checkpoint_scenarios:
            db.save_checkpoint(
                checkpoint["name"],
                checkpoint["phase"], 
                checkpoint["data"]
            )
        
        # Test loading checkpoints
        for checkpoint in checkpoint_scenarios:
            loaded_data = db.load_checkpoint(checkpoint["name"])
            assert loaded_data == checkpoint["data"]
        
        # Test non-existent checkpoint
        missing_checkpoint = db.load_checkpoint("nonexistent_checkpoint")
        assert missing_checkpoint is None
    
    def test_demo_protein_lookup_scenarios(self, tmp_path):
        """Test protein lookup scenarios from demo"""
        db = MitoNetDatabase(str(tmp_path / "lookup_test.db"))
        db.initialize_db()
        
        # Create proteins with various aliases
        protein1 = db.get_or_create_protein("P12345", gene_symbol="ATP1A1")
        protein2 = db.get_or_create_protein("Q67890", gene_symbol="MYOD1")
        
        # Add various types of aliases
        db.add_protein_alias(protein1, 'symbol', 'ATP1A1', None)
        db.add_protein_alias(protein1, 'string', '9606.ENSP00000001', None)
        db.add_protein_alias(protein1, 'ensembl', 'ENSP00000001', None)
        
        db.add_protein_alias(protein2, 'symbol', 'MYOD1', None)
        db.add_protein_alias(protein2, 'string', '9606.ENSP00000002', None)
        
        # Test lookup by UniProt ID
        found_by_uniprot = db.get_protein_by_uniprot("P12345")
        assert found_by_uniprot.id == protein1.id
        
        # Test lookup by gene symbol
        found_by_symbol = db.find_protein_by_alias("ATP1A1", "symbol")
        assert found_by_symbol.id == protein1.id
        
        # Test lookup by STRING ID
        found_by_string = db.find_protein_by_alias("9606.ENSP00000001", "string")
        assert found_by_string.id == protein1.id
        
        # Test lookup by any alias type
        found_by_any = db.find_protein_by_alias("MYOD1")
        assert found_by_any.id == protein2.id
        
        # Test lookup of non-existent protein
        not_found = db.find_protein_by_alias("NONEXISTENT_GENE")
        assert not_found is None
        
        # Verify alias statistics
        stats = db.get_statistics()
        assert stats['num_aliases'] == 5  # 3 for protein1, 2 for protein2
    
    def _create_test_data_files(self, data_dir):
        """Create test data files for demo scenarios"""
        import pandas as pd
        
        # Create directory structure
        (data_dir / "mitocarta").mkdir(parents=True)
        (data_dir / "hpa").mkdir(parents=True)
        
        # Create HPA file
        hpa_data = pd.DataFrame([
            {
                "Gene": "ATP1A1",
                "Gene description": "ATPase Na+/K+ transporting subunit alpha 1",
                "Evidence": "Evidence at protein level",
                "Tissue RNA - skeletal muscle [nTPM]": 45.6,
                "RNA tissue specificity": "Not detected",
                "Subcellular main location": "Plasma membrane"
            },
            {
                "Gene": "MYOD1",
                "Gene description": "Myogenic differentiation 1",
                "Evidence": "Evidence at protein level",
                "Tissue RNA - skeletal muscle [nTPM]": 123.4,
                "RNA tissue specificity": "Tissue enhanced", 
                "Subcellular main location": "Nucleoplasm"
            }
        ])
        
        hpa_file = data_dir / "hpa/hpa_skm.tsv"
        hpa_data.to_csv(hpa_file, sep='\t', index=False)
        
        # Create MitoCarta file (empty, will be mocked)
        mitocarta_file = data_dir / "mitocarta/Human.MitoCarta3.0.xls"
        mitocarta_file.touch()
    
    def _get_mock_mitocarta_data(self):
        """Get mock MitoCarta data for testing"""
        import pandas as pd
        
        return pd.DataFrame([
            {
                "Symbol": "CYC1",
                "Description": "Cytochrome c1, heme protein",
                "MitoCarta3.0_List": "1",
                "MitoCarta3.0_Evidence": "Experimental",
                "MitoCarta3.0_SubMitoLocalization": "IMS",
                "MitoCarta3.0_MitoPathways": "Respiratory chain"
            },
            {
                "Symbol": "ATP5F1A",
                "Description": "ATP synthase F1 subunit alpha", 
                "MitoCarta3.0_List": "1",
                "MitoCarta3.0_Evidence": "Experimental",
                "MitoCarta3.0_SubMitoLocalization": "Matrix",
                "MitoCarta3.0_MitoPathways": "ATP synthesis"
            }
        ])


@pytest.mark.integration
class TestDemoErrorScenarios:
    """Test error scenarios and edge cases from demo context"""
    
    def test_demo_with_missing_data_files(self, tmp_path):
        """Test demo scenarios when data files are missing"""
        db = MitoNetDatabase(str(tmp_path / "missing_files.db"))
        db.initialize_db()
        
        data_dir = tmp_path / "empty_networks"
        data_dir.mkdir()
        
        ingestion = DataIngestionManager(db, data_dir)
        
        # Test that missing files are handled gracefully
        missing_file = data_dir / "nonexistent.txt"
        needs_update = ingestion.needs_update("TEST_SOURCE", missing_file)
        assert needs_update is False  # Missing files don't need update
        
        # Verify database remains stable
        stats = db.get_statistics()
        assert stats['num_proteins'] == 0
        assert stats['num_sources'] == 0
    
    def test_demo_with_corrupted_data(self, tmp_path):
        """Test demo scenarios with corrupted or invalid data"""
        db = MitoNetDatabase(str(tmp_path / "corrupted_data.db"))
        db.initialize_db()
        
        data_dir = tmp_path / "corrupted_networks"
        data_dir.mkdir()
        
        # Create corrupted file
        corrupted_file = data_dir / "corrupted.txt"
        corrupted_file.write_text("This is not valid CSV data\nRandom text\n")
        
        ingestion = DataIngestionManager(db, data_dir)
        
        # Test that corrupted files are handled gracefully
        try:
            # This might raise an exception, which is acceptable
            ingestion.ingest_string_aliases(corrupted_file, version="test")
        except Exception:
            # Exception handling is acceptable for corrupted data
            pass
        
        # Database should remain in consistent state
        stats = db.get_statistics()
        # Should either have no data or should have handled errors gracefully
        assert isinstance(stats['num_proteins'], int)
        assert isinstance(stats['num_sources'], int)
    
    def test_demo_database_recovery(self, tmp_path):
        """Test database recovery scenarios"""
        db_path = tmp_path / "recovery_test.db"
        
        # Create and populate database
        db1 = MitoNetDatabase(str(db_path))
        db1.initialize_db()
        
        protein = db1.get_or_create_protein("P12345", gene_symbol="TEST_RECOVERY")
        db1.save_checkpoint("recovery_test", "test_phase", {"data": "test"})
        
        original_stats = db1.get_statistics()
        
        # Close and reopen database (simulate restart)
        del db1
        
        db2 = MitoNetDatabase(str(db_path))
        # No need to initialize - should already exist
        
        # Verify data persisted
        recovered_stats = db2.get_statistics()
        assert recovered_stats['num_proteins'] == original_stats['num_proteins']
        
        # Verify specific data
        recovered_protein = db2.get_protein_by_uniprot("P12345")
        assert recovered_protein is not None
        assert recovered_protein.gene_symbol == "TEST_RECOVERY"
        
        # Verify checkpoint persisted
        recovered_checkpoint = db2.load_checkpoint("recovery_test")
        assert recovered_checkpoint == {"data": "test"}