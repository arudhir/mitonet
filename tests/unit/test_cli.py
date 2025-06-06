"""
Unit tests for CLI functionality
"""

import pytest
from click.testing import CliRunner
from unittest.mock import Mock, patch
from pathlib import Path

from mitonet.cli import cli
from mitonet.database import MitoNetDatabase


@pytest.mark.cli
class TestCLICommands:
    """Test CLI command functionality"""
    
    @pytest.fixture
    def runner(self):
        """Click testing runner"""
        return CliRunner()
    
    @pytest.fixture 
    def mock_db_path(self, tmp_path):
        """Temporary database path for testing"""
        return str(tmp_path / "test.db")
    
    def test_init_command(self, runner, mock_db_path):
        """Test database initialization command"""
        result = runner.invoke(cli, ['--db-path', mock_db_path, 'init'])
        
        assert result.exit_code == 0
        assert "Database initialized" in result.output
        assert Path(mock_db_path).exists()
    
    def test_status_command_empty_db(self, runner, mock_db_path):
        """Test status command on empty database"""
        # Initialize database first
        runner.invoke(cli, ['--db-path', mock_db_path, 'init'])
        
        # Check status
        result = runner.invoke(cli, ['--db-path', mock_db_path, 'status'])
        
        assert result.exit_code == 0
        assert "Proteins: 0" in result.output
        assert "Interactions: 0" in result.output
        assert "Data sources: 0" in result.output
    
    def test_export_network_empty_db(self, runner, mock_db_path):
        """Test export network command on empty database"""
        # Initialize database
        runner.invoke(cli, ['--db-path', mock_db_path, 'init'])
        
        # Try to export from empty database
        result = runner.invoke(cli, [
            '--db-path', mock_db_path,
            'export-network',
            '--filter-type', 'all'
        ])
        
        assert result.exit_code == 0
        assert "Database is empty" in result.output
    
    def test_export_network_command_help(self, runner, mock_db_path):
        """Test export network command help"""
        result = runner.invoke(cli, ['export-network', '--help'])
        
        assert result.exit_code == 0
        assert "Export filtered networks" in result.output
        assert "--filter-type" in result.output
    
    def test_export_predefined_empty_db(self, runner, mock_db_path):
        """Test export predefined command on empty database"""
        # Initialize database
        runner.invoke(cli, ['--db-path', mock_db_path, 'init'])
        
        # Try to export predefined networks from empty database
        result = runner.invoke(cli, [
            '--db-path', mock_db_path,
            'export-predefined'
        ])
        
        assert result.exit_code == 0
        assert "Database is empty" in result.output
    
    def test_checkpoints_command(self, runner, mock_db_path):
        """Test checkpoints command"""
        # Initialize database
        runner.invoke(cli, ['--db-path', mock_db_path, 'init'])
        
        # Add a checkpoint manually
        db = MitoNetDatabase(mock_db_path)
        db.save_checkpoint("test_checkpoint", "test_phase", {"data": "test"})
        
        # List checkpoints
        result = runner.invoke(cli, ['--db-path', mock_db_path, 'checkpoints'])
        
        assert result.exit_code == 0
        assert "test_checkpoint" in result.output
    
    @patch('mitonet.cli.DataIngestionManager')
    def test_update_command_specific_source(self, mock_ingestion_class, runner, mock_db_path, tmp_path):
        """Test update command for specific source"""
        # Setup mock
        mock_ingestion = Mock()
        mock_ingestion_class.return_value = mock_ingestion
        mock_ingestion.needs_update.return_value = True
        mock_ingestion.ingest_string_aliases.return_value = Mock()
        
        # Create fake data file
        data_dir = tmp_path / "networks" / "string"
        data_dir.mkdir(parents=True)
        test_file = data_dir / "9606.protein.aliases.v12.0.txt.gz"
        test_file.touch()
        
        # Initialize database
        runner.invoke(cli, ['--db-path', mock_db_path, 'init'])
        
        # Test update command
        result = runner.invoke(cli, [
            '--db-path', mock_db_path,
            '--data-dir', str(tmp_path / "networks"),
            'update',
            '--source', 'STRING_aliases'
        ])
        
        assert result.exit_code == 0
        assert "Updating STRING_aliases" in result.output
    
    @patch('mitonet.cli.DataIngestionManager')
    def test_update_command_no_changes(self, mock_ingestion_class, runner, mock_db_path, tmp_path):
        """Test update command when no changes are needed"""
        # Setup mock
        mock_ingestion = Mock()
        mock_ingestion_class.return_value = mock_ingestion
        mock_ingestion.needs_update.return_value = False
        
        # Create fake data file
        data_dir = tmp_path / "networks" / "string"
        data_dir.mkdir(parents=True)
        test_file = data_dir / "9606.protein.aliases.v12.0.txt.gz"
        test_file.touch()
        
        # Initialize database
        runner.invoke(cli, ['--db-path', mock_db_path, 'init'])
        
        # Test update command
        result = runner.invoke(cli, [
            '--db-path', mock_db_path,
            '--data-dir', str(tmp_path / "networks"),
            'update',
            '--source', 'STRING_aliases'
        ])
        
        assert result.exit_code == 0
        assert "No updates needed" in result.output
    
    def test_update_command_unknown_source(self, runner, mock_db_path):
        """Test update command with unknown source"""
        # Initialize database
        runner.invoke(cli, ['--db-path', mock_db_path, 'init'])
        
        # Test unknown source
        result = runner.invoke(cli, [
            '--db-path', mock_db_path,
            'update',
            '--source', 'UNKNOWN_SOURCE'
        ])
        
        assert result.exit_code == 0
        assert "Unknown source" in result.output
    
    def test_update_command_missing_file(self, runner, mock_db_path, tmp_path):
        """Test update command with missing data file"""
        # Initialize database
        runner.invoke(cli, ['--db-path', mock_db_path, 'init'])
        
        # Test with missing file
        result = runner.invoke(cli, [
            '--db-path', mock_db_path,
            '--data-dir', str(tmp_path / "networks"),
            'update',
            '--source', 'STRING_aliases'
        ])
        
        assert result.exit_code == 0
        assert "File not found" in result.output
    
    def test_verbose_flag(self, runner, mock_db_path):
        """Test verbose logging flag"""
        result = runner.invoke(cli, [
            '--db-path', mock_db_path,
            '--verbose',
            'init'
        ])
        
        assert result.exit_code == 0
        # Verbose mode should work without errors
    
    def test_help_command(self, runner):
        """Test help command"""
        result = runner.invoke(cli, ['--help'])
        
        assert result.exit_code == 0
        assert "MitoNet Incremental Update System" in result.output
        assert "Commands:" in result.output
        
        # Test help for specific command
        result = runner.invoke(cli, ['init', '--help'])
        assert result.exit_code == 0
        assert "Initialize the database" in result.output


@pytest.mark.cli
class TestCLIErrorHandling:
    """Test CLI error handling"""
    
    @pytest.fixture
    def runner(self):
        """Click testing runner"""
        return CliRunner()
    
    def test_invalid_database_path(self, runner):
        """Test handling of invalid database path"""
        # Try to use a directory as database path
        result = runner.invoke(cli, ['--db-path', '/invalid/path/db.db', 'status'])
        
        # Should handle error gracefully
        # (specific behavior depends on implementation)
        assert result.exit_code != 0 or "error" in result.output.lower()
    
    def test_export_network_invalid_filter(self, runner, tmp_path):
        """Test export network command with invalid filter type"""
        db_path = str(tmp_path / "test.db")
        runner.invoke(cli, ['--db-path', db_path, 'init'])
        
        result = runner.invoke(cli, [
            '--db-path', db_path,
            'export-network',
            '--filter-type', 'invalid_filter'
        ])
        
        # Should fail with invalid choice
        assert result.exit_code == 2
        assert "Invalid value" in result.output
    
    def test_export_genes_no_genes_specified(self, runner, tmp_path):
        """Test export network with genes filter but no genes specified"""
        db_path = str(tmp_path / "test.db")
        runner.invoke(cli, ['--db-path', db_path, 'init'])
        
        # Add a fake protein to avoid empty database error
        from mitonet.database import MitoNetDatabase
        db = MitoNetDatabase(db_path)
        db.get_or_create_protein(uniprot_id="P12345", gene_symbol="TEST")
        
        result = runner.invoke(cli, [
            '--db-path', db_path,
            'export-network',
            '--filter-type', 'genes'
        ])
        
        assert result.exit_code == 0
        assert "For 'genes' filter, provide --genes or --uniprots" in result.output


@pytest.mark.cli
@pytest.mark.integration
class TestCLIIntegration:
    """Integration tests for CLI with real database operations"""
    
    @pytest.fixture
    def runner(self):
        """Click testing runner"""
        return CliRunner()
    
    def test_full_workflow(self, runner, tmp_path):
        """Test complete workflow: init -> status -> export (empty)"""
        db_path = str(tmp_path / "workflow.db")
        
        # Initialize database
        result = runner.invoke(cli, ['--db-path', db_path, 'init'])
        assert result.exit_code == 0
        
        # Check status of empty database
        result = runner.invoke(cli, ['--db-path', db_path, 'status'])
        assert result.exit_code == 0
        assert "Proteins: 0" in result.output
        
        # Try to export from empty database
        result = runner.invoke(cli, [
            '--db-path', db_path,
            'export-network',
            '--filter-type', 'all'
        ])
        assert result.exit_code == 0
        assert "Database is empty" in result.output
        
        # Check that database file exists
        assert Path(db_path).exists()
    
    def test_persistent_data(self, runner, tmp_path):
        """Test that database persists across CLI invocations"""
        db_path = str(tmp_path / "persistent.db")
        
        # Initialize database
        runner.invoke(cli, ['--db-path', db_path, 'init'])
        
        # Check initial status
        result = runner.invoke(cli, ['--db-path', db_path, 'status'])
        assert result.exit_code == 0
        assert "Proteins: 0" in result.output
        
        # In a new invocation, database should still exist
        result = runner.invoke(cli, ['--db-path', db_path, 'status'])
        assert result.exit_code == 0
        assert "Proteins: 0" in result.output