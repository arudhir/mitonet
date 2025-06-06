# MitoNet Pipeline

A comprehensive mitochondrial protein network integration pipeline with incremental database updates, modular architecture, and comprehensive testing.

## Overview

MitoNet integrates multiple biological databases to create a comprehensive protein interaction database, then allows you to generate customized network views through powerful filtering. The pipeline features:

- **Complete Database Approach**: Ingest ALL protein interactions from all sources into one comprehensive database
- **Flexible Network Filtering**: Generate mitochondrial, muscle, gene-specific, or custom networks from the same data
- **Incremental Updates**: Only process changed data files, with intelligent change detection
- **Modular Architecture**: Clean separation of database, ingestion, filtering, and export components  
- **Comprehensive Testing**: Full test suite with unit, integration, and CLI tests
- **Multiple Data Sources**: STRING, BioGRID, MitoCarta, HPA, and more
- **Multiple Export Formats**: JSON, GraphML, and CSV network formats

## Quick Start

```bash
# Setup development environment
make dev-setup

# Ingest all available data sources
make update-all

# Export predefined networks (mitochondrial, muscle, high-confidence)
make export-predefined

# Check what's in the database
make show-stats
```

## Installation

### Using the Makefile (Recommended)

```bash
# Install dependencies
make install

# Or install with development dependencies
make install-dev
```

### Manual Installation

```bash
# Using uv (recommended)
uv sync

# Or using pip
pip install -e .
```

## Building the Complete Database

### 1. Initialize and Setup

```bash
# Setup development environment (installs deps + initializes DB)
make dev-setup

# Or do it step by step:
make install-dev
make db-init
```

### 2. Ingest Complete Datasets from All Sources

The new architecture ingests ALL proteins and interactions from each data source, creating a comprehensive database that you can filter later:

```bash
# Update all data sources at once
make update-all

# Or update specific sources individually:
make update-string    # STRING protein aliases and interactions (COMPLETE dataset)
make update-biogrid   # BioGRID interactions (COMPLETE dataset)
make update-mitocarta # MitoCarta mitochondrial annotations (ALL proteins)
make update-hpa       # Human Protein Atlas muscle expression (ALL proteins)
```

This step processes complete datasets and may take some time, but only needs to be done once per data source version.

### 3. Monitor Database Building

```bash
# Check database statistics
make show-stats

# View processing checkpoints
make checkpoints
```

You should see the database grow to contain hundreds of thousands of proteins and millions of interactions.

### 4. Export Filtered Networks

Once the complete database is built, you can quickly generate different network views:

```bash
# Export commonly used networks
make export-predefined

# Or export specific network types
make export-mitochondrial    # Only mitochondrial proteins
make export-muscle          # Only muscle-expressed proteins  
make export-all             # Complete network (all proteins)
make export-high-confidence # High-confidence interactions only
```

## Iterating on the Database

### Exploring Different Networks from the Same Data

The power of the complete database approach is that you can generate many different network views without re-ingesting data:

```bash
# Start with a mitochondrial network
make export-mitochondrial

# Then explore muscle-specific networks
make export-muscle

# Focus on specific genes with their neighbors
make export-genes GENES="ATP1A1,MYOD1,CYC1" NEIGHBORS=1

# Try different confidence thresholds
make export-high-confidence CONFIDENCE=0.8
```

### Updating Data Sources

The pipeline intelligently detects which files have changed:

```bash
# This will only process files that have been updated since last run
make update-all

# Force update specific source even if unchanged
make update-string

# Check what data sources are tracked
make show-stats
```

### Iterative Analysis Workflow

```bash
# Build the complete database once
make dev-setup
make update-all

# Then quickly generate different views for analysis
make export-mitochondrial                    # For mitochondrial analysis
make export-muscle                          # For muscle tissue analysis  
make export-genes GENES="ATP1A1,NDUFA1" NEIGHBORS=2  # For specific gene analysis

# When new data becomes available, update and re-export
make update-mitocarta                        # Add new data source
make export-predefined                       # Re-export all common networks

# Check what changed
make show-stats
```

### Custom Network Generation

Generate highly specific networks using advanced filtering:

```bash
# High-confidence mitochondrial network
make export-custom FILTER_TYPE=mitochondrial MIN_CONF=0.8

# Muscle network with experimental evidence only
make export-custom FILTER_TYPE=muscle MIN_CONF=0.5 EVIDENCE_TYPES=experimental

# Gene-specific network with 2-hop neighbors
make export-custom FILTER_TYPE=genes GENES="ATP1A1,COX1" NEIGHBORS=2 MIN_CONF=0.6
```

## Testing

### Running Tests

```bash
# Run all tests
make test

# Run specific test categories
make test-unit        # Unit tests only
make test-integration # Integration tests only
make test-cli         # CLI functionality tests
make test-database    # Database operation tests

# Run tests with coverage report
make test-coverage

# Quick development testing
make quick-test
```

### Test Categories

- **Unit Tests**: Test individual components in isolation
- **Integration Tests**: Test complete workflows and data processing
- **CLI Tests**: Test command-line interface functionality
- **Database Tests**: Test database operations and integrity
- **Slow Tests**: Performance and large dataset tests

### Development Testing Workflow

```bash
# Setup test environment
make dev-setup

# Add test data
make dev-test-data

# Run quick tests during development
make quick-test

# Run full test suite before committing
make test

# Generate coverage report
make test-coverage
# Open htmlcov/index.html to view detailed coverage
```

## Development Workflows

### Complete Development Setup

```bash
# Setup everything for development
make dev-setup

# Add sample data for testing
make dev-test-data

# Run a complete pipeline test
make full-pipeline
```

### Demo Environment

```bash
# Setup demo with sample proteins
make demo

# Check what's available
make show-stats

# Run some updates
make update-string
make generate-network
```

### Clean Development Environment

```bash
# Clean temporary files and caches
make clean

# Reset everything and start fresh
make dev-reset

# Clean output files
make clean-outputs
```

## Data Sources

The pipeline integrates data from:

- **STRING**: Protein-protein interactions and aliases
- **BioGRID**: Curated protein interactions
- **MitoCarta**: Mitochondrial protein annotations
- **Human Protein Atlas (HPA)**: Tissue-specific expression data
- **CORUM**: Protein complex data
- **Reactome**: Pathway interactions

## Output Files

Networks are exported to the `outputs/` directory:

- `mitonet_network.json` - NetworkX JSON format
- `mitonet_network.graphml` - GraphML for Cytoscape/Gephi
- `mitonet_nodes.csv` - Node attributes table
- `mitonet_edges.csv` - Edge list with attributes

## Command Reference

```bash
# View all available commands organized by category
make help

# Key commands:
make dev-setup              # Setup development environment
make update-all             # Ingest all data sources into complete database
make export-predefined      # Export common networks (mitochondrial, muscle, high-confidence)
make export-mitochondrial   # Export mitochondrial protein network
make export-muscle          # Export muscle-expressed protein network
make export-genes           # Export gene-specific network (use GENES="..." NEIGHBORS=1)
make export-custom          # Export custom filtered network (use FILTER_TYPE=...)
make test                   # Run test suite
make show-stats             # Show database statistics
make clean                  # Clean temporary files
```

## Advanced Usage

### Using the CLI Directly

```bash
# The Makefile wraps the CLI, but you can use it directly:
uv run python -m mitonet.cli --help
uv run python -m mitonet.cli init
uv run python -m mitonet.cli update --source STRING_aliases
uv run python -m mitonet.cli export-network --filter-type mitochondrial
uv run python -m mitonet.cli export-predefined
uv run python -m mitonet.cli status
```

### Legacy Pipeline

```bash
# Run the original pipeline scripts
make run-legacy          # Original main.py
make run-memory-optimized # Memory-optimized version
make run-simplified      # Simplified version
```

### Custom Database Path

```bash
# All commands support custom database paths
make db-init DB_PATH=/path/to/custom.db
make db-status DB_PATH=/path/to/custom.db
```

## Project Structure

```
mitonet/
├── mitonet/              # Main package
│   ├── database.py       # Database models and operations
│   ├── ingestion.py      # Data ingestion and processing
│   ├── export.py         # Network filtering and export
│   └── cli.py           # Command-line interface
├── tests/               # Test suite
│   ├── unit/            # Unit tests
│   ├── integration/     # Integration tests
│   └── conftest.py      # Test fixtures
├── networks/            # Raw data files
├── outputs/             # Generated network files
├── Makefile            # Development commands
└── pyproject.toml      # Project configuration
```

## Contributing

1. Setup development environment: `make dev-setup`
2. Add test data: `make dev-test-data`
3. Make your changes
4. Run tests: `make test`
5. Clean up: `make clean`

## Troubleshooting

### Common Issues

**Database locked errors**: Make sure no other processes are using the database
```bash
make clean  # Clean up any stale connections
```

**Memory issues with large files**: The pipeline uses chunked processing, but very large datasets may still cause issues
```bash
make run-memory-optimized  # Use memory-optimized version
```

**Test failures**: Make sure you have test dependencies installed
```bash
make install-dev  # Install test dependencies
make clean        # Clean test cache
make test         # Re-run tests
```

**Missing data files**: Ensure data files are in the correct `networks/` subdirectories
```bash
ls networks/string/     # Check STRING files are present
ls networks/mitocarta/  # Check MitoCarta files are present
```

### Getting Help

```bash
make help                           # Show all available commands
uv run python -m mitonet.cli --help # CLI help
make test --help                    # Pytest help
```
