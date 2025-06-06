# MitoNet Incremental Update System

This branch implements a comprehensive incremental update system for MitoNet, enabling efficient management and updating of biological data without rebuilding the entire network from scratch.

## 🚀 Key Features

### 1. Persistent Database Storage
- **SQLite backend** with SQLAlchemy ORM
- **Versioned data sources** with integrity checking (SHA256 hashes)
- **Comprehensive protein model** with mitochondrial and muscle expression attributes
- **Flexible interaction storage** with multi-source support

### 2. Incremental Data Ingestion
- **Change detection** - only processes files that have been modified
- **Chunked processing** - handles large files without memory overflow
- **Source versioning** - tracks different versions of data sources
- **Checkpointing** - save progress and resume processing

### 3. Command-Line Interface
- **Easy management** of database updates
- **Gene addition** for custom proteins
- **Status monitoring** and statistics
- **Selective updates** by data source

## 📊 Database Schema

### Core Tables
- **`proteins`** - Comprehensive protein information
- **`protein_aliases`** - All protein identifiers and mappings
- **`interactions`** - Protein-protein interactions with confidence scores
- **`data_sources`** - Versioned data source tracking
- **`processing_checkpoints`** - Processing state management

### Key Features
- **Automatic versioning** of data sources
- **Conflict resolution** for data updates
- **Multi-source interaction merging**
- **Flexible attribute storage** with JSON fields

## 🛠 Installation

```bash
# Install dependencies
uv sync

# Initialize database
uv run mitonet init

# Check status
uv run mitonet status
```

## 📖 Usage Examples

### Initialize Database
```bash
uv run mitonet init
```

### Add Custom Genes
```bash
# Add genes by symbol
uv run mitonet add-genes --genes "ATP1A1,MYOD1,CYC1"

# Add proteins by UniProt ID
uv run mitonet add-genes --uniprots "P05023,P15531,P08574"
```

### Update Data Sources
```bash
# Update all sources (only if changed)
uv run mitonet update

# Force update all sources
uv run mitonet update --force

# Update specific source
uv run mitonet update --source STRING_aliases
uv run mitonet update --source MitoCarta
```

### Monitor Status
```bash
# Show database statistics
uv run mitonet status

# List processing checkpoints
uv run mitonet checkpoints

# Show checkpoints for specific phase
uv run mitonet checkpoints --phase "alias_ingestion"
```

## 🔄 Incremental Update Workflow

### 1. Change Detection
The system automatically detects changes by:
- Comparing file hashes (SHA256)
- Checking file sizes
- Tracking modification times
- Version string extraction from filenames

### 2. Data Processing
When updates are needed:
```python
# Check if update needed
if ingestion.needs_update('STRING_aliases', file_path):
    # Process incrementally
    source = ingestion.ingest_string_aliases(file_path)
```

### 3. Conflict Resolution
- **Version tracking** prevents duplicate ingestion
- **Merge strategies** for overlapping data
- **Confidence score optimization** for interactions
- **Attribute updates** with timestamp tracking

## 🧪 Demo Script

Run the demo to see incremental updates in action:

```bash
uv run python demo_incremental.py
```

This demonstrates:
- Database initialization
- Manual gene addition
- Automatic change detection
- Incremental data ingestion
- Statistics reporting

## 📈 Performance Benefits

### Memory Efficiency
- **Chunked processing** - 50K rows at a time
- **Streaming ingestion** - no full file loading
- **Garbage collection** - automatic memory cleanup

### Processing Speed
- **Skip unchanged files** - only process what's new
- **Checkpoint system** - resume from failures
- **Batch operations** - efficient database writes

### Storage Optimization
- **Normalized schema** - eliminate data duplication
- **Compressed JSON** - flexible attribute storage
- **Indexed lookups** - fast protein/interaction queries

## 🔗 Integration with Existing Pipeline

The incremental system is designed to complement the existing pipeline:

```python
# Traditional approach
integrator = MitoNetIntegrator()
integrator.phase2_build_id_mapping()

# Incremental approach
db = MitoNetDatabase()
ingestion = DataIngestionManager(db)
ingestion.ingest_all_sources()
```

## 🎯 Future Enhancements

### Planned Features
1. **Network export** from database to GraphML/JSON
2. **API integration** for real-time data fetching
3. **Automated scheduling** for periodic updates
4. **Conflict resolution UI** for manual curation
5. **Delta exports** for change tracking

### Advanced Capabilities
- **Distributed processing** for large datasets
- **Cloud storage integration** (S3, GCS)
- **Real-time collaboration** features
- **Machine learning** for data quality assessment

## 💡 Key Design Decisions

### Why SQLite?
- **Single-file deployment** - easy to share and backup
- **ACID compliance** - data integrity guarantees
- **Cross-platform** - works everywhere Python does
- **Performance** - sufficient for biological networks

### Why SQLAlchemy?
- **ORM benefits** - object-oriented data access
- **Migration support** - schema evolution
- **Query optimization** - automatic SQL generation
- **Type safety** - prevents data corruption

### Why Chunked Processing?
- **Memory constraints** - handle files larger than RAM
- **Progress tracking** - see processing status
- **Error recovery** - resume from checkpoints
- **Parallel processing** - future scalability

## 🏗 Architecture Overview

```
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│   CLI Interface │────│  Ingestion Mgr  │────│     Database    │
└─────────────────┘    └─────────────────┘    └─────────────────┘
         │                       │                       │
         ▼                       ▼                       ▼
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│  User Commands  │    │  Change Detection│    │  SQLite + ORM   │
│  - init         │    │  - File hashing  │    │  - Proteins     │
│  - update       │    │  - Version check │    │  - Interactions │
│  - add-genes    │    │  - Chunked read  │    │  - Sources      │
│  - status       │    │  - Checkpoints   │    │  - Checkpoints  │
└─────────────────┘    └─────────────────┘    └─────────────────┘
```

This system provides a robust foundation for incremental biological network updates while maintaining data integrity and processing efficiency.