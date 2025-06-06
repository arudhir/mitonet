# MitoNet Testing Framework

This document describes the comprehensive testing framework for the MitoNet incremental update system.

## 🧪 Test Structure

```
tests/
├── conftest.py              # Pytest configuration and fixtures
├── unit/                    # Unit tests for individual components
│   ├── test_database.py     # Database layer tests
│   ├── test_ingestion.py    # Data ingestion tests
│   └── test_cli.py          # CLI functionality tests
└── integration/             # Integration and workflow tests
    ├── test_incremental_workflow.py  # Complete workflow tests
    └── test_demo_converted.py        # Converted demo scenarios
```

## 🚀 Running Tests

### Quick Start
```bash
# Run basic working tests
python test_runner.py basic

# Run all unit tests
python test_runner.py unit

# Run with coverage report
python test_runner.py coverage
```

### Manual pytest Commands
```bash
# Install test dependencies
uv sync --extra test

# Run specific test
uv run pytest tests/unit/test_database.py::TestMitoNetDatabase::test_initialize_database -v

# Run all database tests
uv run pytest tests/unit/test_database.py -v

# Run CLI tests
uv run pytest tests/unit/test_cli.py -v

# Run with coverage
uv run pytest tests/unit/test_database.py --cov=mitonet --cov-report=html

# Run integration tests
uv run pytest tests/integration/ -v -m integration

# Run tests excluding slow ones
uv run pytest tests/ -v -m "not slow"
```

## 📊 Test Categories

### Unit Tests (`tests/unit/`)

**Database Tests (`test_database.py`)**
- ✅ Database initialization
- ✅ Protein creation and retrieval
- ✅ Data source management
- ✅ Checkpoint functionality
- ✅ Statistics generation
- ⚠️ Session management (some tests need fixes)

**CLI Tests (`test_cli.py`)**
- ✅ Help commands
- ✅ Database initialization via CLI
- ✅ Gene addition commands
- ✅ Status reporting
- ✅ Error handling

**Ingestion Tests (`test_ingestion.py`)**
- ✅ File hash calculation
- ✅ Version extraction
- ✅ Change detection
- ✅ Data format handling
- ⚠️ Full ingestion workflow (needs session fixes)

### Integration Tests (`tests/integration/`)

**Workflow Tests (`test_incremental_workflow.py`)**
- ✅ Complete incremental update workflow
- ✅ Change detection and versioning
- ✅ Checkpoint recovery
- ✅ Data integrity across updates
- ✅ Performance with bulk operations

**Demo Scenarios (`test_demo_converted.py`)**
- ✅ Demo workflow converted to tests
- ✅ Gene addition scenarios
- ✅ Data source versioning
- ✅ Error scenarios and recovery

## 🔧 Test Fixtures and Utilities

### Key Fixtures (defined in `conftest.py`)

```python
@pytest.fixture
def temp_db():
    """In-memory SQLite database for testing"""

@pytest.fixture  
def temp_db_file():
    """File-based temporary database"""

@pytest.fixture
def sample_proteins():
    """Sample protein data for testing"""

@pytest.fixture
def populated_db(temp_db, sample_proteins):
    """Database pre-populated with test data"""

@pytest.fixture
def ingestion_manager(temp_db, test_data_dir):
    """Data ingestion manager for testing"""
```

### Test Utilities

**DatabaseAssertions** - Helper methods for database testing:
```python
db_assertions.assert_protein_exists(db, "P12345")
db_assertions.assert_protein_count(db, 10)
db_assertions.assert_interaction_exists(db, "P12345", "Q67890")
```

## 🏷️ Test Markers

Tests are organized with pytest markers:

```python
@pytest.mark.database     # Database-related tests
@pytest.mark.cli          # CLI functionality tests  
@pytest.mark.integration  # Integration tests
@pytest.mark.slow         # Slow-running tests
```

Run specific categories:
```bash
# Only database tests
uv run pytest -m database

# Exclude slow tests
uv run pytest -m "not slow"

# Only integration tests
uv run pytest -m integration
```

## 📈 Coverage Reports

Test coverage is automatically generated:

```bash
# Generate HTML coverage report
uv run pytest --cov=mitonet --cov-report=html

# View report
open htmlcov/index.html
```

Current coverage areas:
- ✅ Database initialization and basic operations
- ✅ CLI command parsing and help
- ✅ Configuration and fixtures
- ⚠️ Complex database operations (session management)
- ⚠️ Data ingestion workflows
- ⚠️ Error handling edge cases

## 🐛 Known Issues and Fixes Needed

### Session Management Issues
Several tests fail due to SQLAlchemy session management:
```
DetachedInstanceError: Instance <Protein> is not bound to a Session
```

**Fix Strategy:**
1. Modify database methods to return data instead of ORM objects
2. Use session merging for cross-session object access
3. Add proper session context managers

### Test Data Management
Some integration tests need:
1. Better mock data file creation
2. Cleanup of temporary files
3. More realistic data scenarios

## ✅ Verified Test Scenarios

### Working Test Cases
1. **Database initialization** - ✅ Passes
2. **Checkpoint save/load** - ✅ Passes  
3. **CLI help system** - ✅ Passes
4. **Basic statistics** - ✅ Passes
5. **Test fixtures** - ✅ Working

### Test Coverage Highlights
- Database schema creation and validation
- CLI command structure and help
- Configuration file parsing
- Test fixture setup and teardown
- Error handling for missing dependencies

## 🎯 Testing Best Practices

### Writing New Tests
1. **Use appropriate fixtures** - `temp_db` for database tests
2. **Follow naming conventions** - `test_*` functions
3. **Use descriptive docstrings** - Explain what is being tested
4. **Add appropriate markers** - `@pytest.mark.database`, etc.
5. **Clean up resources** - Tests should be independent

### Test Data
1. **Use fixtures for sample data** - Don't hardcode in tests
2. **Make tests deterministic** - Avoid random data
3. **Test edge cases** - Empty data, invalid formats, etc.
4. **Mock external dependencies** - File I/O, network calls

### Example Test Pattern
```python
@pytest.mark.database
def test_new_functionality(self, temp_db, db_assertions):
    """Test description of what this tests"""
    # Arrange
    protein = temp_db.get_or_create_protein("P12345", gene_symbol="TEST")
    
    # Act  
    result = temp_db.some_operation(protein)
    
    # Assert
    db_assertions.assert_protein_exists(temp_db, "P12345")
    assert result.expected_attribute == "expected_value"
```

## 🚀 Future Testing Enhancements

### Planned Improvements
1. **Fix session management** in database tests
2. **Add performance benchmarks** for large datasets  
3. **Integration with CI/CD** pipeline
4. **Property-based testing** with Hypothesis
5. **Database migration tests** for schema changes
6. **Load testing** for concurrent operations

### Additional Test Coverage
1. **Real data file parsing** tests
2. **Network error simulation** tests
3. **Concurrent access** tests
4. **Memory usage** validation
5. **Cross-platform** compatibility tests

This testing framework provides a solid foundation for ensuring the reliability and correctness of the MitoNet incremental update system while supporting future development and refactoring efforts.