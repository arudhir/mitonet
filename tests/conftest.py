"""
Pytest configuration and fixtures for MitoNet tests
"""

import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock
import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from mitonet.database import MitoNetDatabase, Base, Protein, DataSource, Interaction, ProteinAlias
from mitonet.ingestion import DataIngestionManager


@pytest.fixture(scope="session")
def test_data_dir():
    """Create a temporary directory with test data files"""
    temp_dir = Path(tempfile.mkdtemp())
    
    # Create directory structure
    (temp_dir / "string").mkdir()
    (temp_dir / "mitocarta").mkdir()
    (temp_dir / "hpa").mkdir()
    (temp_dir / "biogrid").mkdir()
    (temp_dir / "reactome").mkdir()
    (temp_dir / "corum").mkdir()
    
    yield temp_dir
    
    # Cleanup
    shutil.rmtree(temp_dir)


@pytest.fixture
def temp_db():
    """Create a temporary in-memory SQLite database for testing"""
    db = MitoNetDatabase(":memory:")
    db.initialize_db()
    yield db


@pytest.fixture
def temp_db_file():
    """Create a temporary file-based SQLite database for testing"""
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as f:
        db_path = f.name
    
    db = MitoNetDatabase(db_path)
    db.initialize_db()
    yield db
    
    # Cleanup
    Path(db_path).unlink(missing_ok=True)


@pytest.fixture
def sample_proteins():
    """Sample protein data for testing"""
    return [
        {
            "uniprot_id": "P12345",
            "gene_symbol": "ATP1A1",
            "gene_description": "ATPase Na+/K+ transporting subunit alpha 1",
            "is_mitochondrial": False,
            "is_muscle_expressed": True,
            "muscle_tpm": 45.6,
            "priority_score": 2.3
        },
        {
            "uniprot_id": "Q67890",
            "gene_symbol": "MYOD1", 
            "gene_description": "Myogenic differentiation 1",
            "is_mitochondrial": False,
            "is_muscle_expressed": True,
            "muscle_tpm": 123.4,
            "priority_score": 3.1
        },
        {
            "uniprot_id": "P00123",
            "gene_symbol": "CYC1",
            "gene_description": "Cytochrome c1, heme protein",
            "is_mitochondrial": True,
            "is_muscle_expressed": False,
            "muscle_tpm": 0.0,
            "priority_score": 4.5
        }
    ]


@pytest.fixture
def populated_db(temp_db, sample_proteins):
    """Database populated with sample data"""
    for protein_data in sample_proteins:
        protein = temp_db.get_or_create_protein(**protein_data)
        
        # Add some aliases
        temp_db.add_protein_alias(protein, 'symbol', protein_data['gene_symbol'], None)
        temp_db.add_protein_alias(protein, 'string', f"9606.ENSP{protein_data['uniprot_id']}", None)
    
    yield temp_db


@pytest.fixture
def sample_string_aliases_data():
    """Sample STRING aliases data"""
    return pd.DataFrame([
        {"#string_protein_id": "9606.ENSP00000000233", "alias": "P12345", "source": "UniProt_AC"},
        {"#string_protein_id": "9606.ENSP00000000233", "alias": "ATP1A1", "source": "UniProt_GN"},
        {"#string_protein_id": "9606.ENSP00000001234", "alias": "Q67890", "source": "UniProt_AC"},
        {"#string_protein_id": "9606.ENSP00000001234", "alias": "MYOD1", "source": "BLAST_UniProt_GN"},
        {"#string_protein_id": "9606.ENSP00000005678", "alias": "P00123", "source": "UniProt_AC"},
        {"#string_protein_id": "9606.ENSP00000005678", "alias": "CYC1", "source": "UniProt_GN"},
    ])


@pytest.fixture
def sample_string_interactions_data():
    """Sample STRING interactions data"""
    return pd.DataFrame([
        {
            "protein1": "9606.ENSP00000000233",
            "protein2": "9606.ENSP00000001234", 
            "neighborhood": 0,
            "fusion": 0,
            "cooccurence": 150,
            "coexpression": 200,
            "experimental": 450,
            "database": 300,
            "textmining": 100,
            "combined_score": 650
        },
        {
            "protein1": "9606.ENSP00000001234",
            "protein2": "9606.ENSP00000005678",
            "neighborhood": 0,
            "fusion": 0, 
            "cooccurence": 0,
            "coexpression": 350,
            "experimental": 200,
            "database": 150,
            "textmining": 80,
            "combined_score": 500
        }
    ])


@pytest.fixture
def sample_mitocarta_data():
    """Sample MitoCarta data"""
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


@pytest.fixture
def sample_hpa_data():
    """Sample HPA skeletal muscle data"""
    return pd.DataFrame([
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
        },
        {
            "Gene": "CYC1",
            "Gene description": "Cytochrome c1, heme protein",
            "Evidence": "Evidence at protein level",
            "Tissue RNA - skeletal muscle [nTPM]": 0.0,  # Not muscle expressed
            "RNA tissue specificity": "Not detected",
            "Subcellular main location": "Mitochondria"
        }
    ])


@pytest.fixture
def mock_files(test_data_dir, sample_string_aliases_data, sample_mitocarta_data, sample_hpa_data):
    """Create mock data files for testing"""
    files = {}
    
    # STRING aliases file
    aliases_file = test_data_dir / "string/9606.protein.aliases.v12.0.txt.gz"
    sample_string_aliases_data.to_csv(aliases_file, sep='\t', index=False, compression='gzip')
    files['string_aliases'] = aliases_file
    
    # MitoCarta file (Excel format is complex, so we'll mock this differently in specific tests)
    mitocarta_file = test_data_dir / "mitocarta/Human.MitoCarta3.0.xls"
    files['mitocarta'] = mitocarta_file
    
    # HPA file
    hpa_file = test_data_dir / "hpa/hpa_skm.tsv"
    sample_hpa_data.to_csv(hpa_file, sep='\t', index=False)
    files['hpa'] = hpa_file
    
    return files


@pytest.fixture
def ingestion_manager(temp_db, test_data_dir):
    """Data ingestion manager with temporary database and test data"""
    return DataIngestionManager(temp_db, test_data_dir)


class DatabaseAssertions:
    """Helper class for database-related assertions"""
    
    @staticmethod
    def assert_protein_exists(db: MitoNetDatabase, uniprot_id: str):
        """Assert that a protein exists in the database"""
        protein = db.get_protein_by_uniprot(uniprot_id)
        assert protein is not None, f"Protein {uniprot_id} not found in database"
        return protein
    
    @staticmethod
    def assert_protein_count(db: MitoNetDatabase, expected_count: int):
        """Assert the total number of proteins in database"""
        stats = db.get_statistics()
        assert stats['num_proteins'] == expected_count, \
            f"Expected {expected_count} proteins, found {stats['num_proteins']}"
    
    @staticmethod
    def assert_interaction_exists(db: MitoNetDatabase, uniprot1: str, uniprot2: str):
        """Assert that an interaction exists between two proteins"""
        session = db.get_session()
        try:
            protein1 = db.get_protein_by_uniprot(uniprot1)
            protein2 = db.get_protein_by_uniprot(uniprot2)
            assert protein1 and protein2, f"One or both proteins not found: {uniprot1}, {uniprot2}"
            
            # Check both directions
            interaction = session.query(Interaction).filter(
                ((Interaction.protein1_id == protein1.id) & (Interaction.protein2_id == protein2.id)) |
                ((Interaction.protein1_id == protein2.id) & (Interaction.protein2_id == protein1.id))
            ).first()
            
            assert interaction is not None, f"Interaction not found between {uniprot1} and {uniprot2}"
            return interaction
        finally:
            session.close()
    
    @staticmethod
    def assert_data_source_exists(db: MitoNetDatabase, name: str, version: str):
        """Assert that a data source exists with given name and version"""
        session = db.get_session()
        try:
            source = session.query(DataSource).filter_by(name=name, version=version).first()
            assert source is not None, f"Data source {name} v{version} not found"
            return source
        finally:
            session.close()


@pytest.fixture
def db_assertions():
    """Provide database assertion helpers"""
    return DatabaseAssertions


# Custom pytest markers for organizing tests
pytest.mark.database = pytest.mark.database
pytest.mark.cli = pytest.mark.cli
pytest.mark.integration = pytest.mark.integration
pytest.mark.slow = pytest.mark.slow