"""
Integration tests for the complete incremental update workflow
"""

import pytest
import tempfile
import shutil
from pathlib import Path
import pandas as pd
from unittest.mock import patch

from mitonet.database import MitoNetDatabase
from mitonet.ingestion import DataIngestionManager


@pytest.mark.integration
@pytest.mark.slow
class TestIncrementalWorkflow:
    """Test complete incremental update workflow"""
    
    @pytest.fixture
    def workflow_setup(self, tmp_path):
        """Setup for workflow testing"""
        # Create database
        db_path = tmp_path / "workflow.db"
        db = MitoNetDatabase(str(db_path))
        db.initialize_db()
        
        # Create data directory
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        
        # Create ingestion manager
        ingestion = DataIngestionManager(db, data_dir)
        
        return {
            'db': db,
            'data_dir': data_dir,
            'ingestion': ingestion,
            'db_path': db_path
        }
    
    def test_empty_to_populated_workflow(self, workflow_setup):
        """Test workflow from empty database to populated with real data"""
        db = workflow_setup['db']
        ingestion = workflow_setup['ingestion']
        data_dir = workflow_setup['data_dir']
        
        # 1. Start with empty database
        stats = db.get_statistics()
        assert stats['num_proteins'] == 0
        assert stats['num_interactions'] == 0
        
        # 2. Add some manual genes
        protein1 = db.get_or_create_protein(
            uniprot_id="P12345",
            gene_symbol="ATP1A1",
            is_muscle_expressed=True,
            muscle_tpm=45.6
        )
        
        protein2 = db.get_or_create_protein(
            uniprot_id="Q67890", 
            gene_symbol="MYOD1",
            is_muscle_expressed=True,
            muscle_tpm=123.4
        )
        
        protein3 = db.get_or_create_protein(
            uniprot_id="P00123",
            gene_symbol="CYC1",
            is_mitochondrial=True
        )
        
        # Add aliases for data source linking
        db.add_protein_alias(protein1, 'symbol', 'ATP1A1', None)
        db.add_protein_alias(protein2, 'symbol', 'MYOD1', None)
        db.add_protein_alias(protein3, 'symbol', 'CYC1', None)
        
        # 3. Verify manual additions
        stats = db.get_statistics()
        assert stats['num_proteins'] == 3
        assert stats['num_mitochondrial'] == 1
        assert stats['num_muscle_expressed'] == 2
        
        # 4. Create mock data files
        self._create_mock_hpa_file(data_dir)
        self._create_mock_mitocarta_file(data_dir)
        
        # 5. Test incremental updates
        # HPA update
        hpa_file = data_dir / "hpa/hpa_skm.tsv"
        hpa_source = ingestion.ingest_hpa_muscle(hpa_file)
        assert hpa_source.name == "HPA_muscle"
        
        # MitoCarta update (mocked)
        with patch('pandas.read_excel') as mock_excel:
            mock_excel.return_value = self._get_mock_mitocarta_data()
            mitocarta_file = data_dir / "mitocarta/Human.MitoCarta3.0.xls"
            mitocarta_file.touch()
            mitocarta_source = ingestion.ingest_mitocarta(mitocarta_file)
            assert mitocarta_source.name == "MitoCarta"
        
        # 6. Verify data was integrated
        stats = db.get_statistics()
        assert stats['num_proteins'] == 3
        assert stats['num_sources'] == 2
        
        # Check that protein attributes were updated
        updated_protein1 = db.get_protein_by_uniprot("P12345")
        assert updated_protein1.is_muscle_expressed is True
        assert updated_protein1.protein_evidence_level == "Evidence at protein level"
        
        updated_protein3 = db.get_protein_by_uniprot("P00123")
        assert updated_protein3.is_mitochondrial is True
        assert updated_protein3.mitocarta_sub_localization == "IMS"
    
    def test_incremental_update_detection(self, workflow_setup):
        """Test that incremental updates only process changed files"""
        ingestion = workflow_setup['ingestion']
        data_dir = workflow_setup['data_dir']
        
        # Create initial file
        test_file = data_dir / "test_file.txt"
        test_file.write_text("Initial content")
        
        # First check - should need update (new file)
        assert ingestion.needs_update("TEST_SOURCE", test_file, "1.0") is True
        
        # Register the file
        source = ingestion.db.get_or_create_data_source(
            name="TEST_SOURCE",
            version="1.0",
            file_path=str(test_file),
            file_size=test_file.stat().st_size,
            file_hash=ingestion.calculate_file_hash(test_file)
        )
        
        # Second check - should not need update (unchanged)
        assert ingestion.needs_update("TEST_SOURCE", test_file, "1.0") is False
        
        # Modify file
        test_file.write_text("Modified content")
        
        # Third check - should need update (changed)
        assert ingestion.needs_update("TEST_SOURCE", test_file, "1.0") is True
    
    def test_checkpoint_recovery(self, workflow_setup):
        """Test checkpoint saving and recovery"""
        db = workflow_setup['db']
        
        # Save a checkpoint
        test_data = {
            'proteins_processed': 100,
            'interactions_added': 500,
            'current_phase': 'string_ingestion'
        }
        
        db.save_checkpoint("test_recovery", "phase_4", test_data, status="completed")
        
        # Verify checkpoint can be loaded
        loaded_data = db.load_checkpoint("test_recovery")
        assert loaded_data == test_data
        
        # Test checkpoint history
        session = db.get_session()
        try:
            from mitonet.database import ProcessingCheckpoint
            checkpoints = session.query(ProcessingCheckpoint).filter_by(
                checkpoint_name="test_recovery"
            ).all()
            assert len(checkpoints) == 1
            assert checkpoints[0].status == "completed"
        finally:
            session.close()
    
    def test_version_tracking(self, workflow_setup):
        """Test that data source versions are tracked correctly"""
        db = workflow_setup['db']
        
        # Create multiple versions of the same source
        source_v1 = db.get_or_create_data_source(
            name="STRING",
            version="v11.5",
            file_path="/path/to/v11.5",
            file_hash="hash_v11_5"
        )
        
        source_v2 = db.get_or_create_data_source(
            name="STRING", 
            version="v12.0",
            file_path="/path/to/v12.0",
            file_hash="hash_v12_0"
        )
        
        # Should be different sources
        assert source_v1.id != source_v2.id
        assert source_v1.version == "v11.5"
        assert source_v2.version == "v12.0"
        
        # Check statistics
        stats = db.get_statistics()
        assert stats['num_sources'] == 2
        assert 'STRING vv11.5' in stats['data_sources']
        assert 'STRING vv12.0' in stats['data_sources']
    
    def test_data_integrity_across_updates(self, workflow_setup):
        """Test that data integrity is maintained across multiple updates"""
        db = workflow_setup['db']
        
        # Create initial protein
        protein = db.get_or_create_protein(
            uniprot_id="P12345",
            gene_symbol="TEST_GENE",
            is_mitochondrial=False,
            muscle_tpm=10.0
        )
        initial_created_at = protein.created_at
        
        # Update the same protein with new information
        updated_protein = db.get_or_create_protein(
            uniprot_id="P12345",  # Same UniProt ID
            gene_symbol="TEST_GENE",
            is_mitochondrial=True,  # Changed
            muscle_tpm=15.0,  # Changed
            protein_evidence_level="Evidence at protein level"  # New
        )
        
        # Should be the same protein object, updated
        assert updated_protein.id == protein.id
        assert updated_protein.uniprot_id == "P12345"
        assert updated_protein.is_mitochondrial is True
        assert updated_protein.muscle_tpm == 15.0
        assert updated_protein.protein_evidence_level == "Evidence at protein level"
        assert updated_protein.created_at == initial_created_at  # Should not change
        assert updated_protein.updated_at > initial_created_at  # Should be updated
        
        # Verify database still has only one protein
        stats = db.get_statistics()
        assert stats['num_proteins'] == 1
    
    def test_bulk_operations_performance(self, workflow_setup):
        """Test performance with bulk operations"""
        db = workflow_setup['db']
        
        # Add many proteins
        num_proteins = 1000
        for i in range(num_proteins):
            db.get_or_create_protein(
                uniprot_id=f"P{i:05d}",
                gene_symbol=f"GENE_{i}",
                is_mitochondrial=(i % 2 == 0),
                muscle_tpm=float(i)
            )
        
        # Verify all were created
        stats = db.get_statistics()
        assert stats['num_proteins'] == num_proteins
        assert stats['num_mitochondrial'] == num_proteins // 2
        
        # Test bulk alias addition
        session = db.get_session()
        try:
            from mitonet.database import Protein, ProteinAlias
            proteins = session.query(Protein).limit(100).all()
            
            for protein in proteins:
                db.add_protein_alias(protein, 'symbol', protein.gene_symbol, None)
            
            # Verify aliases were added
            stats = db.get_statistics()
            assert stats['num_aliases'] >= 100
        finally:
            session.close()
    
    def _create_mock_hpa_file(self, data_dir):
        """Create mock HPA data file"""
        hpa_dir = data_dir / "hpa"
        hpa_dir.mkdir()
        
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
        
        hpa_file = hpa_dir / "hpa_skm.tsv"
        hpa_data.to_csv(hpa_file, sep='\t', index=False)
    
    def _create_mock_mitocarta_file(self, data_dir):
        """Create mock MitoCarta data file"""
        mitocarta_dir = data_dir / "mitocarta"
        mitocarta_dir.mkdir()
        
        # Just create the file - content will be mocked
        mitocarta_file = mitocarta_dir / "Human.MitoCarta3.0.xls"
        mitocarta_file.touch()
    
    def _get_mock_mitocarta_data(self):
        """Get mock MitoCarta DataFrame"""
        return pd.DataFrame([
            {
                "Symbol": "CYC1",
                "Description": "Cytochrome c1, heme protein",
                "MitoCarta3.0_List": "1",
                "MitoCarta3.0_Evidence": "Experimental",
                "MitoCarta3.0_SubMitoLocalization": "IMS",
                "MitoCarta3.0_MitoPathways": "Respiratory chain"
            }
        ])


@pytest.mark.integration
class TestRealDataIngestion:
    """Integration tests with real data file formats (using sample data)"""
    
    def test_string_aliases_format_handling(self, temp_db, tmp_path):
        """Test handling of real STRING aliases file format"""
        # Create realistic STRING aliases data
        aliases_data = pd.DataFrame([
            {"#string_protein_id": "9606.ENSP00000000233", "alias": "P31946", "source": "UniProt_AC"},
            {"#string_protein_id": "9606.ENSP00000000233", "alias": "1433B_HUMAN", "source": "UniProt_ID"},
            {"#string_protein_id": "9606.ENSP00000000233", "alias": "YWHAB", "source": "UniProt_GN"},
            {"#string_protein_id": "9606.ENSP00000000233", "alias": "YWHAB", "source": "BLAST_UniProt_GN"},
        ])
        
        # Save as compressed file
        aliases_file = tmp_path / "string_aliases.txt.gz"
        aliases_data.to_csv(aliases_file, sep='\t', index=False, compression='gzip')
        
        # Test ingestion
        ingestion = DataIngestionManager(temp_db, tmp_path)
        source = ingestion.ingest_string_aliases(aliases_file, version="test")
        
        # Verify results
        assert source.name == "STRING_aliases"
        protein = temp_db.get_protein_by_uniprot("P31946")
        assert protein is not None
        assert protein.gene_symbol == "YWHAB"
    
    def test_string_interactions_format_handling(self, temp_db, tmp_path):
        """Test handling of real STRING interactions file format"""
        # Create realistic STRING interactions data
        interactions_data = pd.DataFrame([
            {
                "protein1": "9606.ENSP00000000233",
                "protein2": "9606.ENSP00000002125",
                "neighborhood": 0,
                "fusion": 0,
                "cooccurence": 332,
                "coexpression": 59,
                "experimental": 190,
                "database": 900,
                "textmining": 597,
                "combined_score": 944
            }
        ])
        
        # Save as space-separated compressed file (STRING format)
        interactions_file = tmp_path / "string_interactions.txt.gz"
        interactions_data.to_csv(interactions_file, sep=' ', index=False, compression='gzip')
        
        # Pre-populate with proteins
        protein1 = temp_db.get_or_create_protein(uniprot_id="P31946")
        protein2 = temp_db.get_or_create_protein(uniprot_id="P25007")
        temp_db.add_protein_alias(protein1, 'string', '9606.ENSP00000000233', None)
        temp_db.add_protein_alias(protein2, 'string', '9606.ENSP00000002125', None)
        
        # Test ingestion
        ingestion = DataIngestionManager(temp_db, tmp_path)
        source = ingestion.ingest_string_interactions(interactions_file, "STRING_test", version="test")
        
        # Verify results
        assert source.name == "STRING_test"
        stats = temp_db.get_statistics()
        assert stats['num_interactions'] == 1


@pytest.mark.integration
@pytest.mark.slow
class TestLargeDatasetHandling:
    """Test handling of large datasets (simulated)"""
    
    def test_chunked_processing_memory_efficiency(self, temp_db, tmp_path):
        """Test that chunked processing handles large files efficiently"""
        # Create a larger dataset
        large_data = []
        for i in range(10000):  # 10k rows
            large_data.append({
                "#string_protein_id": f"9606.ENSP{i:08d}",
                "alias": f"P{i:05d}", 
                "source": "UniProt_AC"
            })
        
        large_df = pd.DataFrame(large_data)
        large_file = tmp_path / "large_aliases.txt.gz"
        large_df.to_csv(large_file, sep='\t', index=False, compression='gzip')
        
        # Test ingestion with small chunk size
        ingestion = DataIngestionManager(temp_db, tmp_path)
        source = ingestion.ingest_string_aliases(large_file, version="large_test")
        
        # Verify all data was processed
        stats = temp_db.get_statistics()
        assert stats['num_proteins'] == 10000
        assert source.name == "STRING_aliases"
    
    def test_error_recovery_during_ingestion(self, temp_db, tmp_path):
        """Test error recovery and checkpoint functionality"""
        # Create data with some problematic rows
        mixed_data = pd.DataFrame([
            {"#string_protein_id": "9606.ENSP00000001", "alias": "P12345", "source": "UniProt_AC"},
            {"#string_protein_id": "9606.ENSP00000002", "alias": "Q67890", "source": "UniProt_AC"},
            # This would cause issues if we had strict validation
            {"#string_protein_id": "", "alias": "", "source": ""},  # Empty row
            {"#string_protein_id": "9606.ENSP00000004", "alias": "P00123", "source": "UniProt_AC"},
        ])
        
        mixed_file = tmp_path / "mixed_data.txt.gz"
        mixed_data.to_csv(mixed_file, sep='\t', index=False, compression='gzip')
        
        # Test ingestion - should handle problematic rows gracefully
        ingestion = DataIngestionManager(temp_db, tmp_path)
        source = ingestion.ingest_string_aliases(mixed_file, version="mixed_test")
        
        # Should have processed the valid rows
        stats = temp_db.get_statistics()
        assert stats['num_proteins'] >= 3  # At least the valid ones