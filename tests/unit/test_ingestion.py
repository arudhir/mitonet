"""
Unit tests for data ingestion functionality
"""

import pytest
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch, mock_open
import pandas as pd

from mitonet.ingestion import DataIngestionManager
from mitonet.database import MitoNetDatabase


@pytest.mark.database
class TestDataIngestionManager:
    """Test cases for DataIngestionManager"""
    
    def test_calculate_file_hash(self, ingestion_manager, test_data_dir):
        """Test file hash calculation"""
        # Create a test file
        test_file = test_data_dir / "test.txt"
        test_file.write_text("Hello, World!")
        
        hash1 = ingestion_manager.calculate_file_hash(test_file)
        hash2 = ingestion_manager.calculate_file_hash(test_file)
        
        # Same file should produce same hash
        assert hash1 == hash2
        assert len(hash1) == 64  # SHA256 produces 64-character hex string
        
        # Different content should produce different hash
        test_file.write_text("Hello, Universe!")
        hash3 = ingestion_manager.calculate_file_hash(test_file)
        assert hash3 != hash1
    
    def test_extract_version_from_filename(self, ingestion_manager):
        """Test version extraction from filenames"""
        # STRING files
        assert ingestion_manager._extract_version_from_filename(
            Path("9606.protein.aliases.v12.0.txt.gz")
        ) == "v12.0"
        
        # BioGRID files
        assert ingestion_manager._extract_version_from_filename(
            Path("BIOGRID-ALL-4.4.246.tab3.txt")
        ) == "4.4.246"
        
        # MitoCarta files
        assert ingestion_manager._extract_version_from_filename(
            Path("Human.MitoCarta3.0.xls")
        ) == "3.0"
    
    def test_needs_update_new_file(self, ingestion_manager, test_data_dir):
        """Test needs_update for a new file"""
        test_file = test_data_dir / "new_file.txt"
        test_file.write_text("New content")
        
        # New file should always need update
        assert ingestion_manager.needs_update("TEST_SOURCE", test_file, "1.0") is True
    
    def test_needs_update_unchanged_file(self, ingestion_manager, test_data_dir):
        """Test needs_update for unchanged file"""
        test_file = test_data_dir / "unchanged.txt"
        test_file.write_text("Unchanged content")
        
        # First time - needs update
        assert ingestion_manager.needs_update("TEST_SOURCE", test_file, "1.0") is True
        
        # Register the file as a data source
        ingestion_manager.db.get_or_create_data_source(
            name="TEST_SOURCE",
            version="1.0", 
            file_path=str(test_file),
            file_size=test_file.stat().st_size,
            file_hash=ingestion_manager.calculate_file_hash(test_file)
        )
        
        # Second time - no update needed
        assert ingestion_manager.needs_update("TEST_SOURCE", test_file, "1.0") is False
    
    def test_needs_update_changed_file(self, ingestion_manager, test_data_dir):
        """Test needs_update for changed file"""
        test_file = test_data_dir / "changed.txt"
        test_file.write_text("Original content")
        
        # Register original file
        ingestion_manager.db.get_or_create_data_source(
            name="TEST_SOURCE",
            version="1.0",
            file_path=str(test_file),
            file_size=test_file.stat().st_size,
            file_hash=ingestion_manager.calculate_file_hash(test_file)
        )
        
        # Change the file
        test_file.write_text("Modified content")
        
        # Should need update now
        assert ingestion_manager.needs_update("TEST_SOURCE", test_file, "1.0") is True


@pytest.mark.database
@pytest.mark.slow
class TestStringAliasIngestion:
    """Test STRING alias ingestion"""
    
    def test_ingest_string_aliases(self, ingestion_manager, test_data_dir, sample_string_aliases_data):
        """Test STRING aliases ingestion"""
        # Create test file
        aliases_file = test_data_dir / "string_aliases.txt.gz"
        sample_string_aliases_data.to_csv(aliases_file, sep='\t', index=False, compression='gzip')
        
        # Ingest data
        source = ingestion_manager.ingest_string_aliases(aliases_file, version="test")
        
        # Verify data source was created
        assert source.name == "STRING_aliases"
        assert source.version == "test"
        
        # Verify proteins were created
        stats = ingestion_manager.db.get_statistics()
        assert stats['num_proteins'] == 3  # P12345, Q67890, P00123
        
        # Verify aliases were created
        assert stats['num_aliases'] > 0
        
        # Test specific protein lookup
        protein = ingestion_manager.db.get_protein_by_uniprot("P12345")
        assert protein is not None
        assert protein.gene_symbol == "ATP1A1"


@pytest.mark.database
class TestMitoCartaIngestion:
    """Test MitoCarta data ingestion"""
    
    @patch('pandas.read_excel')
    def test_ingest_mitocarta(self, mock_read_excel, ingestion_manager, test_data_dir, sample_mitocarta_data):
        """Test MitoCarta data ingestion"""
        # Mock Excel file reading
        mock_read_excel.return_value = sample_mitocarta_data
        
        # Create test file
        mitocarta_file = test_data_dir / "mitocarta.xls"
        mitocarta_file.touch()
        
        # Pre-populate database with proteins that can be found by symbol
        ingestion_manager.db.get_or_create_protein(uniprot_id="P00123", gene_symbol="CYC1")
        ingestion_manager.db.add_protein_alias(
            ingestion_manager.db.get_protein_by_uniprot("P00123"),
            'symbol', 'CYC1', None
        )
        
        # Ingest MitoCarta data
        source = ingestion_manager.ingest_mitocarta(mitocarta_file, version="test")
        
        # Verify data source
        assert source.name == "MitoCarta"
        assert source.version == "test"
        
        # Verify mitochondrial annotation
        protein = ingestion_manager.db.get_protein_by_uniprot("P00123")
        assert protein.is_mitochondrial is True
        assert protein.mitocarta_sub_localization == "IMS"


@pytest.mark.database
class TestHPAIngestion:
    """Test HPA data ingestion"""
    
    def test_ingest_hpa_muscle(self, ingestion_manager, test_data_dir, sample_hpa_data):
        """Test HPA skeletal muscle data ingestion"""
        # Create test file
        hpa_file = test_data_dir / "hpa_muscle.tsv"
        sample_hpa_data.to_csv(hpa_file, sep='\t', index=False)
        
        # Pre-populate database with proteins
        for _, row in sample_hpa_data.iterrows():
            if row['Tissue RNA - skeletal muscle [nTPM]'] > 0:  # Only muscle-expressed
                protein = ingestion_manager.db.get_or_create_protein(
                    uniprot_id=f"TEST_{row['Gene']}",
                    gene_symbol=row['Gene']
                )
                ingestion_manager.db.add_protein_alias(protein, 'symbol', row['Gene'], None)
        
        # Ingest HPA data
        source = ingestion_manager.ingest_hpa_muscle(hpa_file, version="test")
        
        # Verify data source
        assert source.name == "HPA_muscle"
        assert source.version == "test"
        
        # Verify muscle expression annotation
        atp1a1_protein = ingestion_manager.db.find_protein_by_alias("ATP1A1", "symbol")
        assert atp1a1_protein.is_muscle_expressed is True
        assert atp1a1_protein.muscle_tpm == 45.6
        
        myod1_protein = ingestion_manager.db.find_protein_by_alias("MYOD1", "symbol")
        assert myod1_protein.is_muscle_expressed is True
        assert myod1_protein.muscle_tpm == 123.4


@pytest.mark.database
class TestStringInteractionIngestion:
    """Test STRING interaction ingestion"""
    
    def test_ingest_string_interactions(self, ingestion_manager, test_data_dir, sample_string_interactions_data):
        """Test STRING interaction data ingestion"""
        # Create test file
        interactions_file = test_data_dir / "string_interactions.txt.gz"
        sample_string_interactions_data.to_csv(interactions_file, sep=' ', index=False, compression='gzip')
        
        # Pre-populate database with proteins and their STRING aliases
        protein1 = ingestion_manager.db.get_or_create_protein(uniprot_id="P12345")
        protein2 = ingestion_manager.db.get_or_create_protein(uniprot_id="Q67890") 
        protein3 = ingestion_manager.db.get_or_create_protein(uniprot_id="P00123")
        
        ingestion_manager.db.add_protein_alias(protein1, 'string', '9606.ENSP00000000233', None)
        ingestion_manager.db.add_protein_alias(protein2, 'string', '9606.ENSP00000001234', None)
        ingestion_manager.db.add_protein_alias(protein3, 'string', '9606.ENSP00000005678', None)
        
        # Ingest interactions
        source = ingestion_manager.ingest_string_interactions(interactions_file, "STRING_test", version="test")
        
        # Verify data source
        assert source.name == "STRING_test"
        assert source.version == "test"
        
        # Verify interactions were created
        stats = ingestion_manager.db.get_statistics()
        assert stats['num_interactions'] == 2
        
        # Test specific interaction
        session = ingestion_manager.db.get_session()
        try:
            from mitonet.database import Interaction
            interaction = session.query(Interaction).first()
            assert interaction.confidence_score == 0.65  # 650/1000
            assert interaction.evidence_type in ['experimental', 'database', 'text_mining']
        finally:
            session.close()


@pytest.mark.database
class TestIngestionErrorHandling:
    """Test error handling in data ingestion"""
    
    def test_missing_file(self, ingestion_manager):
        """Test handling of missing files"""
        missing_file = Path("/nonexistent/file.txt")
        assert ingestion_manager.needs_update("TEST", missing_file) is False
    
    def test_invalid_file_format(self, ingestion_manager, test_data_dir):
        """Test handling of invalid file formats"""
        # Create invalid file
        invalid_file = test_data_dir / "invalid.txt"
        invalid_file.write_text("This is not valid data format")
        
        # Should handle gracefully (implementation dependent)
        try:
            source = ingestion_manager.ingest_string_aliases(invalid_file, version="test")
        except Exception:
            # Exception is acceptable for invalid format
            pass
    
    def test_empty_file(self, ingestion_manager, test_data_dir):
        """Test handling of empty files"""
        empty_file = test_data_dir / "empty.txt.gz"
        pd.DataFrame().to_csv(empty_file, compression='gzip', index=False)
        
        # Should handle empty files gracefully
        source = ingestion_manager.ingest_string_aliases(empty_file, version="test")
        assert source.name == "STRING_aliases"


@pytest.mark.database
class TestClassifyStringEvidence:
    """Test STRING evidence classification"""
    
    def test_classify_string_evidence(self, ingestion_manager):
        """Test STRING evidence type classification"""
        # Experimental dominates
        row = pd.Series({
            'experimental': 500,
            'database': 200,
            'textmining': 100
        })
        assert ingestion_manager._classify_string_evidence(row) == 'experimental'
        
        # Database dominates over textmining
        row = pd.Series({
            'experimental': 100,
            'database': 400,
            'textmining': 200
        })
        assert ingestion_manager._classify_string_evidence(row) == 'database'
        
        # Textmining by default
        row = pd.Series({
            'experimental': 100,
            'database': 150,
            'textmining': 300
        })
        assert ingestion_manager._classify_string_evidence(row) == 'text_mining'