"""
Unit tests for database functionality
"""

import pytest
import tempfile
from pathlib import Path
from datetime import datetime

from mitonet.database import MitoNetDatabase, Protein, DataSource, Interaction, ProteinAlias


@pytest.mark.database
class TestMitoNetDatabase:
    """Test cases for MitoNetDatabase class"""
    
    def test_initialize_database(self, temp_db):
        """Test database initialization"""
        # Database should be initialized and accessible
        stats = temp_db.get_statistics()
        assert stats['num_proteins'] == 0
        assert stats['num_interactions'] == 0
        assert stats['num_sources'] == 0
    
    def test_create_protein(self, temp_db):
        """Test creating a new protein"""
        protein = temp_db.get_or_create_protein(
            uniprot_id="P12345",
            gene_symbol="TEST1",
            gene_description="Test protein 1",
            is_mitochondrial=True,
            priority_score=5.0
        )
        
        # Re-fetch to avoid session issues
        retrieved = temp_db.get_protein_by_uniprot("P12345")
        assert retrieved.uniprot_id == "P12345"
        assert retrieved.gene_symbol == "TEST1"
        assert retrieved.is_mitochondrial is True
        assert retrieved.priority_score == 5.0
    
    def test_get_existing_protein(self, temp_db):
        """Test retrieving an existing protein"""
        # Create protein
        original = temp_db.get_or_create_protein(uniprot_id="P12345", gene_symbol="TEST1")
        
        # Retrieve the same protein
        retrieved = temp_db.get_protein_by_uniprot("P12345")
        
        assert retrieved is not None
        assert retrieved.uniprot_id == original.uniprot_id
        assert retrieved.gene_symbol == original.gene_symbol
    
    def test_update_existing_protein(self, temp_db):
        """Test updating an existing protein"""
        # Create protein
        temp_db.get_or_create_protein(uniprot_id="P12345", gene_symbol="TEST1")
        
        # Update with new attributes
        updated = temp_db.get_or_create_protein(
            uniprot_id="P12345",
            gene_symbol="TEST1",
            is_mitochondrial=True,
            muscle_tpm=10.5
        )
        
        assert updated.is_mitochondrial is True
        assert updated.muscle_tpm == 10.5
    
    def test_add_protein_alias(self, temp_db):
        """Test adding protein aliases"""
        protein = temp_db.get_or_create_protein(uniprot_id="P12345", gene_symbol="TEST1")
        
        # Create a data source
        source = temp_db.get_or_create_data_source(
            name="TEST_SOURCE",
            version="1.0",
            file_path="/test/path"
        )
        
        # Add alias
        temp_db.add_protein_alias(protein, 'string', '9606.ENSP12345', source)
        
        # Verify alias was added
        found_protein = temp_db.find_protein_by_alias('9606.ENSP12345', 'string')
        assert found_protein is not None
        assert found_protein.uniprot_id == "P12345"
    
    def test_find_protein_by_alias(self, temp_db):
        """Test finding protein by alias"""
        protein = temp_db.get_or_create_protein(uniprot_id="P12345", gene_symbol="TEST1")
        temp_db.add_protein_alias(protein, 'symbol', 'TEST1', None)
        
        # Find by symbol
        found = temp_db.find_protein_by_alias('TEST1', 'symbol')
        assert found is not None
        assert found.uniprot_id == "P12345"
        
        # Find by any alias type
        found_any = temp_db.find_protein_by_alias('TEST1')
        assert found_any is not None
        assert found_any.uniprot_id == "P12345"
        
        # Test non-existent alias
        not_found = temp_db.find_protein_by_alias('NONEXISTENT')
        assert not_found is None
    
    def test_add_interaction(self, temp_db):
        """Test adding protein interactions"""
        # Create proteins
        protein1 = temp_db.get_or_create_protein(uniprot_id="P12345", gene_symbol="TEST1")
        protein2 = temp_db.get_or_create_protein(uniprot_id="Q67890", gene_symbol="TEST2")
        
        # Create data source
        source = temp_db.get_or_create_data_source(
            name="TEST_SOURCE",
            version="1.0", 
            file_path="/test/path"
        )
        
        # Add interaction
        interaction = temp_db.add_interaction(
            protein1=protein1,
            protein2=protein2,
            source=source,
            confidence_score=0.8,
            evidence_type="experimental",
            interaction_type="physical"
        )
        
        assert interaction.confidence_score == 0.8
        assert interaction.evidence_type == "experimental"
        assert interaction.interaction_type == "physical"
    
    def test_interaction_ordering(self, temp_db):
        """Test that interactions maintain consistent protein ordering"""
        protein1 = temp_db.get_or_create_protein(uniprot_id="P12345", gene_symbol="TEST1")
        protein2 = temp_db.get_or_create_protein(uniprot_id="Q67890", gene_symbol="TEST2")
        
        source = temp_db.get_or_create_data_source("TEST", "1.0", "/test")
        
        # Add interaction in both directions - should be stored consistently
        interaction1 = temp_db.add_interaction(protein1, protein2, source, 0.8)
        interaction2 = temp_db.add_interaction(protein2, protein1, source, 0.9)
        
        # Should be the same interaction (updated)
        assert interaction1.id == interaction2.id
        assert interaction2.confidence_score == 0.9  # Updated to higher score
    
    def test_data_source_creation(self, temp_db):
        """Test data source creation and retrieval"""
        source = temp_db.get_or_create_data_source(
            name="STRING",
            version="v12.0",
            file_path="/path/to/string.txt",
            file_size=1024,
            file_hash="abc123"
        )
        
        assert source.name == "STRING"
        assert source.version == "v12.0"
        assert source.file_size == 1024
        assert source.file_hash == "abc123"
        
        # Get the same source again
        same_source = temp_db.get_or_create_data_source("STRING", "v12.0", "/path/to/string.txt")
        assert same_source.id == source.id
    
    def test_save_and_load_checkpoint(self, temp_db):
        """Test checkpoint functionality"""
        test_data = {"processed": 1000, "errors": 5}
        
        # Save checkpoint
        temp_db.save_checkpoint("test_checkpoint", "test_phase", test_data)
        
        # Load checkpoint
        loaded_data = temp_db.load_checkpoint("test_checkpoint")
        assert loaded_data == test_data
        
        # Test non-existent checkpoint
        not_found = temp_db.load_checkpoint("nonexistent")
        assert not_found is None
    
    def test_get_statistics(self, populated_db):
        """Test database statistics"""
        stats = populated_db.get_statistics()
        
        assert stats['num_proteins'] == 3
        assert stats['num_mitochondrial'] == 1
        assert stats['num_muscle_expressed'] == 2
        assert stats['num_aliases'] == 6  # 2 aliases per protein


@pytest.mark.database
class TestDatabaseModels:
    """Test database model behavior"""
    
    def test_protein_model(self):
        """Test Protein model attributes"""
        protein = Protein(
            uniprot_id="P12345",
            gene_symbol="TEST",
            is_mitochondrial=True,
            muscle_tpm=15.5
        )
        
        assert protein.uniprot_id == "P12345"
        assert protein.gene_symbol == "TEST" 
        assert protein.is_mitochondrial is True
        assert protein.muscle_tpm == 15.5
        assert protein.created_at is not None
    
    def test_data_source_model(self):
        """Test DataSource model attributes"""
        source = DataSource(
            name="TEST_SOURCE",
            version="1.0",
            file_path="/test/path",
            file_size=1024
        )
        
        assert source.name == "TEST_SOURCE"
        assert source.version == "1.0"
        assert source.file_path == "/test/path"
        assert source.file_size == 1024
        assert source.last_updated is not None
    
    def test_interaction_model(self):
        """Test Interaction model attributes"""
        interaction = Interaction(
            protein1_id=1,
            protein2_id=2,
            confidence_score=0.85,
            evidence_type="experimental",
            source_id=1
        )
        
        assert interaction.protein1_id == 1
        assert interaction.protein2_id == 2
        assert interaction.confidence_score == 0.85
        assert interaction.evidence_type == "experimental"
        assert interaction.created_at is not None


@pytest.mark.database
class TestDatabaseIntegrity:
    """Test database integrity and constraints"""
    
    def test_unique_uniprot_constraint(self, temp_db):
        """Test that UniProt IDs must be unique"""
        temp_db.get_or_create_protein(uniprot_id="P12345", gene_symbol="TEST1")
        
        # Creating another protein with same UniProt ID should return the existing one
        protein2 = temp_db.get_or_create_protein(uniprot_id="P12345", gene_symbol="TEST2")
        
        # Should be the same protein, updated with new gene symbol
        stats = temp_db.get_statistics()
        assert stats['num_proteins'] == 1
        assert protein2.gene_symbol == "TEST2"  # Should be updated
    
    def test_unique_data_source_constraint(self, temp_db):
        """Test that data source name+version combinations must be unique"""
        source1 = temp_db.get_or_create_data_source("TEST", "1.0", "/path1")
        source2 = temp_db.get_or_create_data_source("TEST", "1.0", "/path2")
        
        # Should be the same source
        assert source1.id == source2.id
        
        # Different version should create new source
        source3 = temp_db.get_or_create_data_source("TEST", "2.0", "/path3")
        assert source3.id != source1.id
    
    def test_interaction_uniqueness_per_source(self, temp_db):
        """Test that interactions are unique per source"""
        protein1 = temp_db.get_or_create_protein(uniprot_id="P12345")
        protein2 = temp_db.get_or_create_protein(uniprot_id="Q67890")
        source = temp_db.get_or_create_data_source("TEST", "1.0", "/path")
        
        # Add same interaction twice
        interaction1 = temp_db.add_interaction(protein1, protein2, source, 0.5)
        interaction2 = temp_db.add_interaction(protein1, protein2, source, 0.8)
        
        # Should be the same interaction, with updated confidence
        assert interaction1.id == interaction2.id
        assert interaction2.confidence_score == 0.8