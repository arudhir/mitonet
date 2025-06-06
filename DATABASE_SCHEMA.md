# MitoNet Database Schema

This document describes the complete database schema for the MitoNet protein interaction network pipeline. The database uses SQLAlchemy ORM with SQLite backend and implements a comprehensive approach where all protein data is ingested first, then filtered for specific network exports.

## Overview

The MitoNet database follows a **complete ingestion → selective export** architecture:

1. **Ingest all data sources** into a comprehensive database
2. **Apply filters** to generate specific networks (mitochondrial, muscle-expressed, gene sets, etc.)
3. **Export networks** in multiple formats (JSON, GraphML, CSV)

This approach enables flexible network analysis while maintaining data integrity and performance.

## Core Tables

### DataSource Table

Tracks versions and metadata for each data source to enable incremental updates.

```sql
CREATE TABLE data_sources (
    id INTEGER PRIMARY KEY,
    name VARCHAR(50) NOT NULL,           -- 'STRING', 'BioGRID', 'MitoCarta', etc.
    version VARCHAR(20) NOT NULL,        -- 'v12.0', '4.4.246', '3.0', etc.
    file_path VARCHAR(500) NOT NULL,     -- Path to source file
    file_size INTEGER,                   -- File size in bytes
    file_hash VARCHAR(64),               -- SHA256 hash for integrity checking
    last_updated DATETIME,               -- Last update timestamp
    source_metadata JSON,                -- Additional source-specific metadata
    
    CONSTRAINT _source_version_uc UNIQUE (name, version)
);
```

**Key Features:**
- **Incremental updates**: Uses file hash and size to detect changes
- **Version tracking**: Supports multiple versions of the same data source
- **Metadata storage**: Flexible JSON field for source-specific information

**Indexes:**
- `idx_source_name` on `name` column

### Protein Table

Comprehensive protein information combining data from multiple sources.

```sql
CREATE TABLE proteins (
    id INTEGER PRIMARY KEY,
    uniprot_id VARCHAR(20) UNIQUE NOT NULL,  -- Primary protein identifier
    gene_symbol VARCHAR(50),                 -- Gene symbol (HGNC)
    gene_description TEXT,                   -- Protein/gene description
    
    -- MitoCarta attributes
    is_mitochondrial BOOLEAN DEFAULT FALSE,
    mitocarta_list VARCHAR(20),              -- MitoCarta list membership
    mitocarta_evidence VARCHAR(100),         -- Evidence for mitochondrial localization
    mitocarta_sub_localization VARCHAR(100), -- Sub-mitochondrial localization
    mitocarta_pathways TEXT,                 -- Associated mitochondrial pathways
    
    -- HPA (Human Protein Atlas) attributes
    is_muscle_expressed BOOLEAN DEFAULT FALSE,
    muscle_tpm FLOAT DEFAULT 0.0,           -- Muscle tissue expression (nTPM)
    protein_evidence_level VARCHAR(50),      -- Protein evidence level
    main_localization VARCHAR(100),         -- Main subcellular localization
    
    -- Derived attributes
    priority_score FLOAT,                   -- Calculated priority score
    
    -- Tracking
    created_at DATETIME,
    updated_at DATETIME,
    
    -- Flexible attributes
    extra_attributes JSON                   -- Additional attributes as needed
);
```

**Key Features:**
- **UniProt-centric**: Uses UniProt ID as primary identifier
- **Multi-source integration**: Combines MitoCarta, HPA, and other data
- **Boolean flags**: Fast filtering for mitochondrial and muscle-expressed proteins
- **Extensible**: JSON field for future attributes

**Indexes:**
- `uniprot_id` (unique index)
- `gene_symbol`
- `is_mitochondrial`
- `is_muscle_expressed`
- `priority_score`

### ProteinAlias Table

Stores all protein identifiers and aliases from different sources.

```sql
CREATE TABLE protein_aliases (
    id INTEGER PRIMARY KEY,
    protein_id INTEGER NOT NULL,         -- Foreign key to proteins.id
    alias_type VARCHAR(20) NOT NULL,     -- 'string', 'ensembl', 'symbol', etc.
    alias_value VARCHAR(100) NOT NULL,   -- The actual alias/identifier
    source_id INTEGER,                   -- Foreign key to data_sources.id
    
    FOREIGN KEY (protein_id) REFERENCES proteins (id),
    FOREIGN KEY (source_id) REFERENCES data_sources (id)
);
```

**Key Features:**
- **Unified ID mapping**: Links all protein identifiers to UniProt IDs
- **Source tracking**: Records which data source provided each alias
- **Flexible types**: Supports STRING IDs, Ensembl IDs, gene symbols, etc.

**Indexes:**
- `idx_alias_value` on `alias_value`
- `idx_alias_type_value` on `(alias_type, alias_value)`

### Interaction Table

Stores protein-protein interactions with confidence scores and evidence.

```sql
CREATE TABLE interactions (
    id INTEGER PRIMARY KEY,
    protein1_id INTEGER NOT NULL,        -- Foreign key to proteins.id (smaller ID)
    protein2_id INTEGER NOT NULL,        -- Foreign key to proteins.id (larger ID)
    
    -- Confidence and evidence
    confidence_score FLOAT,              -- Normalized confidence (0.0-1.0)
    evidence_type VARCHAR(50),           -- 'experimental', 'database', 'text_mining'
    interaction_type VARCHAR(50),        -- 'physical', 'functional'
    
    -- Source tracking
    source_id INTEGER NOT NULL,          -- Foreign key to data_sources.id
    source_specific_id VARCHAR(100),     -- Original interaction ID from source
    
    -- Detailed scores
    source_scores JSON,                  -- Source-specific confidence scores
    
    -- Metadata
    created_at DATETIME,
    
    FOREIGN KEY (protein1_id) REFERENCES proteins (id),
    FOREIGN KEY (protein2_id) REFERENCES proteins (id),
    FOREIGN KEY (source_id) REFERENCES data_sources (id),
    
    CONSTRAINT _interaction_source_uc UNIQUE (protein1_id, protein2_id, source_id)
);
```

**Key Features:**
- **Consistent ordering**: `protein1_id < protein2_id` to avoid duplicates
- **Confidence scoring**: Normalized scores for cross-source comparison
- **Evidence classification**: Categorizes evidence types
- **Source tracking**: Links interactions to originating data sources
- **Detailed scores**: JSON storage for source-specific metrics

**Indexes:**
- `idx_proteins` on `(protein1_id, protein2_id)`
- `idx_confidence` on `confidence_score`

### NetworkVersion Table

Tracks different versions of the integrated network with metadata.

```sql
CREATE TABLE network_versions (
    id INTEGER PRIMARY KEY,
    version VARCHAR(20) UNIQUE NOT NULL,  -- Version identifier
    description TEXT,                     -- Human-readable description
    created_at DATETIME,
    
    -- Network statistics
    num_nodes INTEGER,                    -- Number of proteins
    num_edges INTEGER,                    -- Number of interactions
    num_mitochondrial INTEGER,            -- Number of mitochondrial proteins
    num_muscle_expressed INTEGER,         -- Number of muscle-expressed proteins
    
    -- Processing metadata
    parameters JSON,                      -- Processing parameters used
    data_sources_used JSON               -- List of data source IDs used
);
```

**Key Features:**
- **Version tracking**: Maintains history of network builds
- **Statistics**: Precomputed network metrics
- **Reproducibility**: Stores parameters and data sources used

### ProcessingCheckpoint Table

Stores checkpoints for resumable processing of large datasets.

```sql
CREATE TABLE processing_checkpoints (
    id INTEGER PRIMARY KEY,
    checkpoint_name VARCHAR(100) NOT NULL,  -- Unique checkpoint identifier
    phase VARCHAR(50) NOT NULL,             -- 'id_mapping', 'network_integration', etc.
    status VARCHAR(20) DEFAULT 'in_progress', -- 'completed', 'failed', 'in_progress'
    created_at DATETIME,
    completed_at DATETIME,
    
    data JSON,                              -- Intermediate results
    error_message TEXT                      -- Error details if failed
);
```

**Key Features:**
- **Resumable processing**: Enables recovery from interruptions
- **Progress tracking**: Monitors long-running operations
- **Error handling**: Captures failure details for debugging

## Database Management

### MitoNetDatabase Class

The `MitoNetDatabase` class provides high-level operations:

**Core Methods:**
- `initialize_db()`: Create all tables
- `get_session()`: Get database session
- `get_statistics()`: Compute database statistics

**Protein Management:**
- `get_or_create_protein(uniprot_id, **attributes)`: Create or update protein
- `get_protein_by_uniprot(uniprot_id)`: Find protein by UniProt ID
- `find_protein_by_alias(alias_value, alias_type)`: Find protein by any alias
- `add_protein_alias(protein, alias_type, alias_value, source)`: Add protein alias

**Interaction Management:**
- `add_interaction(protein1, protein2, source, confidence_score, **attributes)`: Add interaction
- Uses consistent protein ordering to prevent duplicates
- Implements session management to avoid DetachedInstanceError

**Data Source Management:**
- `get_or_create_data_source(name, version, file_path, **metadata)`: Manage data sources

**Processing Support:**
- `save_checkpoint(name, phase, data, status)`: Save processing checkpoint
- `load_checkpoint(name)`: Load processing checkpoint

## Data Ingestion Architecture

### Incremental Updates

The database supports incremental updates using:

1. **File hashing**: SHA256 checksums detect file changes
2. **Version tracking**: Multiple versions of the same source
3. **Chunked processing**: Large files processed in 50K row chunks
4. **Checkpointing**: Resumable processing for large datasets

### Session Management

Critical session management patterns prevent SQLAlchemy DetachedInstanceError:

```python
# Re-fetch objects within processing loops
fresh_source = self.db.get_or_create_data_source(...)
protein = session.merge(protein)  # Attach to current session
```

### Data Source Integration

**STRING Database:**
- Protein aliases (`9606.protein.aliases.v12.0.txt.gz`)
- Protein interactions (`9606.protein.links.detailed.v12.0.txt.gz`)
- Physical interactions (`9606.protein.physical.links.detailed.v12.0.txt.gz`)

**MitoCarta 3.0:**
- Mitochondrial protein annotations (`Human.MitoCarta3.0.xls`)
- Evidence codes and subcellular localization

**Human Protein Atlas (HPA):**
- Skeletal muscle expression data (`hpa_skm.tsv`)
- Protein evidence levels and localization

## Network Export Architecture

### Filtering System

The `NetworkFilter` class provides flexible filtering:

```python
# Predefined filters
NetworkFilter.mitochondrial_network(min_confidence=0.4)
NetworkFilter.muscle_network(min_confidence=0.4)
NetworkFilter.gene_set_network(genes=['ATP1A1', 'MYOD1'], include_neighbors=1)
NetworkFilter.high_confidence_network(min_confidence=0.7)

# Custom filters
filter = NetworkFilter()
filter.min_confidence = 0.5
filter.evidence_types = {'experimental', 'database'}
filter.min_degree = 3
```

### Export Formats

**JSON Format:**
```json
{
  "nodes": [{"id": "P12345", "gene_symbol": "ATP1A1", "is_mitochondrial": true}],
  "edges": [{"source": "P12345", "target": "Q67890", "confidence": 0.85}]
}
```

**GraphML Format:**
- NetworkX-compatible XML format
- Node and edge attributes preserved
- Compatible with Cytoscape, Gephi

**CSV Format:**
- `*_nodes.csv`: Protein information
- `*_edges.csv`: Interaction information

## Performance Considerations

### Indexing Strategy

Optimized indexes for common queries:
- Protein lookups by UniProt ID and gene symbol
- Alias-based protein finding
- Confidence-based interaction filtering
- Boolean flag filtering (mitochondrial, muscle-expressed)

### Memory Management

- **Chunked processing**: 50K row chunks for large files
- **Session cleanup**: Proper session closing prevents memory leaks
- **Lazy loading**: Relationships loaded on demand

### Query Optimization

- **Consistent protein ordering**: Avoids duplicate interactions
- **Bulk operations**: Batch processing for efficiency
- **Prepared statements**: SQLAlchemy ORM optimization

## Usage Examples

### Basic Operations

```python
# Initialize database
db = MitoNetDatabase("mitonet.db")
db.initialize_db()

# Get statistics
stats = db.get_statistics()
print(f"Proteins: {stats['num_proteins']}")
print(f"Interactions: {stats['num_interactions']}")

# Find protein
protein = db.find_protein_by_alias("ATP1A1", "symbol")
if protein:
    print(f"UniProt ID: {protein.uniprot_id}")
    print(f"Mitochondrial: {protein.is_mitochondrial}")
```

### Network Export

```python
from mitonet.export import NetworkExporter, NetworkFilter

# Create mitochondrial network
filter = NetworkFilter.mitochondrial_network(min_confidence=0.5)
exporter = NetworkExporter(db, Path("outputs"))

files = exporter.export_network(
    network_filter=filter,
    format_types=['json', 'graphml', 'csv'],
    filename_prefix='mitochondrial_network'
)
```

### CLI Usage

```bash
# Initialize and update database
mitonet init
mitonet update

# Export networks
mitonet export-network --filter-type mitochondrial --min-confidence 0.5
mitonet export-network --filter-type genes --genes "ATP1A1,MYOD1" --neighbors 1
mitonet export-predefined

# Check status
mitonet status
mitonet checkpoints
```

## Migration and Versioning

The database schema is designed for evolution:

1. **JSON fields**: Flexible attribute storage
2. **Version tracking**: Data source versioning
3. **Incremental updates**: Non-destructive data updates
4. **Checkpoint system**: Resumable migrations

Future schema changes can be implemented through:
- Additional JSON attributes in existing tables
- New tables with foreign key relationships
- Versioned ingestion pipelines

This schema provides a robust foundation for protein interaction network analysis while maintaining flexibility for future enhancements and data sources.