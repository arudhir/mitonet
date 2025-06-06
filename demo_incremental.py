#!/usr/bin/env python3
"""
Demo script showing the incremental update capabilities
"""

import logging
from pathlib import Path
from mitonet.database import MitoNetDatabase
from mitonet.ingestion import DataIngestionManager

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def demo_incremental_updates():
    """Demonstrate incremental update functionality"""
    
    print("🚀 MitoNet Incremental Update Demo")
    print("=" * 50)
    
    # Initialize database
    db = MitoNetDatabase("demo.db")
    db.initialize_db()
    
    print("\n1. Initial database state:")
    stats = db.get_statistics()
    print(f"   Proteins: {stats['num_proteins']}")
    print(f"   Interactions: {stats['num_interactions']}")
    print(f"   Data sources: {stats['num_sources']}")
    
    # Add some genes manually
    print("\n2. Adding some genes manually...")
    db.get_or_create_protein("P12345", gene_symbol="ATP1A1", is_muscle_expressed=True)
    db.get_or_create_protein("Q67890", gene_symbol="MYOD1", is_muscle_expressed=True)
    db.get_or_create_protein("P00123", gene_symbol="CYC1", is_mitochondrial=True)
    
    stats = db.get_statistics()
    print(f"   After adding genes: {stats['num_proteins']} proteins")
    
    # Initialize ingestion manager
    ingestion = DataIngestionManager(db, Path("networks"))
    
    # Check what needs updating
    print("\n3. Checking which data sources need updating...")
    
    source_files = {
        'STRING_aliases': Path('networks/string/9606.protein.aliases.v12.0.txt.gz'),
        'MitoCarta': Path('networks/mitocarta/Human.MitoCarta3.0.xls'),
        'HPA_muscle': Path('networks/hpa/hpa_skm.tsv'),
    }
    
    for source_name, file_path in source_files.items():
        if file_path.exists():
            needs_update = ingestion.needs_update(source_name, file_path)
            print(f"   {source_name}: {'UPDATE NEEDED' if needs_update else 'up to date'}")
        else:
            print(f"   {source_name}: FILE NOT FOUND")
    
    # Simulate data source updates
    print("\n4. Simulating incremental updates...")
    
    if source_files['MitoCarta'].exists():
        print("   Updating MitoCarta data...")
        try:
            source = ingestion.ingest_mitocarta(source_files['MitoCarta'])
            print(f"   ✅ MitoCarta updated (version: {source.version})")
        except Exception as e:
            print(f"   ❌ MitoCarta update failed: {e}")
    
    if source_files['HPA_muscle'].exists():
        print("   Updating HPA muscle data...")
        try:
            source = ingestion.ingest_hpa_muscle(source_files['HPA_muscle'])
            print(f"   ✅ HPA muscle updated (version: {source.version})")
        except Exception as e:
            print(f"   ❌ HPA muscle update failed: {e}")
    
    # Show final statistics
    print("\n5. Final database state:")
    stats = db.get_statistics()
    print(f"   Proteins: {stats['num_proteins']}")
    print(f"   - Mitochondrial: {stats['num_mitochondrial']}")
    print(f"   - Muscle-expressed: {stats['num_muscle_expressed']}")
    print(f"   Interactions: {stats['num_interactions']}")
    print(f"   Protein aliases: {stats['num_aliases']}")
    print(f"   Data sources: {stats['num_sources']}")
    
    if stats['data_sources']:
        print("\n   Data Sources:")
        for source in stats['data_sources']:
            print(f"     - {source}")
    
    # Show some example proteins
    print("\n6. Example proteins in database:")
    session = db.get_session()
    try:
        from mitonet.database import Protein
        proteins = session.query(Protein).limit(5).all()
        for protein in proteins:
            print(f"   {protein.uniprot_id} ({protein.gene_symbol})")
            print(f"     Mitochondrial: {protein.is_mitochondrial}")
            print(f"     Muscle-expressed: {protein.is_muscle_expressed}")
    finally:
        session.close()
    
    print("\n🎉 Demo complete!")
    print(f"Database saved as: demo.db")
    print("You can inspect it using: mitonet --db-path demo.db status")

if __name__ == "__main__":
    demo_incremental_updates()