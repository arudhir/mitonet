"""
Command-line interface for MitoNet incremental updates
"""

import click
import logging
from pathlib import Path
from .database import MitoNetDatabase
from .ingestion import DataIngestionManager
from .export import NetworkExporter, NetworkFilter, export_predefined_networks

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

@click.group()
@click.option('--db-path', default='mitonet.db', help='Database file path')
@click.option('--data-dir', default='networks', help='Data directory path')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose logging')
@click.pass_context
def cli(ctx, db_path, data_dir, verbose):
    """MitoNet Incremental Update System"""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Initialize database and ingestion manager
    db = MitoNetDatabase(db_path)
    ctx.ensure_object(dict)
    ctx.obj['db'] = db
    ctx.obj['data_dir'] = Path(data_dir)
    ctx.obj['ingestion'] = DataIngestionManager(db, Path(data_dir))

@cli.command()
@click.pass_context
def init(ctx):
    """Initialize the database"""
    db = ctx.obj['db']
    db.initialize_db()
    click.echo(f"Database initialized at {db.db_path}")

@cli.command()
@click.pass_context
def status(ctx):
    """Show database status and statistics"""
    db = ctx.obj['db']
    stats = db.get_statistics()
    
    click.echo("=== MitoNet Database Status ===")
    click.echo(f"Proteins: {stats['num_proteins']:,}")
    click.echo(f"  - Mitochondrial: {stats['num_mitochondrial']:,}")
    click.echo(f"  - Muscle-expressed: {stats['num_muscle_expressed']:,}")
    click.echo(f"Interactions: {stats['num_interactions']:,}")
    click.echo(f"Protein aliases: {stats['num_aliases']:,}")
    click.echo(f"Data sources: {stats['num_sources']:,}")
    
    if stats['data_sources']:
        click.echo("\nData Sources:")
        for source in stats['data_sources']:
            click.echo(f"  - {source}")

@cli.command()
@click.option('--source', help='Specific source to check (STRING_aliases, MitoCarta, etc.)')
@click.option('--force', is_flag=True, help='Force update even if no changes detected')
@click.pass_context
def update(ctx, source, force):
    """Update data sources incrementally"""
    ingestion = ctx.obj['ingestion']
    
    if source:
        # Update specific source
        file_mapping = {
            'STRING_aliases': 'string/9606.protein.aliases.v12.0.txt.gz',
            'STRING_info': 'string/9606.protein.info.v12.0.txt.gz', 
            'STRING_full': 'string/9606.protein.links.detailed.v12.0.txt.gz',
            'STRING_physical': 'string/9606.protein.physical.links.detailed.v12.0.txt.gz',
            'MitoCarta': 'mitocarta/Human.MitoCarta3.0.xls',
            'HPA_muscle': 'hpa/hpa_skm.tsv',
        }
        
        if source not in file_mapping:
            click.echo(f"Unknown source: {source}")
            click.echo(f"Available sources: {', '.join(file_mapping.keys())}")
            return
            
        file_path = ctx.obj['data_dir'] / file_mapping[source]
        
        if not file_path.exists():
            click.echo(f"File not found: {file_path}")
            return
            
        if force or ingestion.needs_update(source, file_path):
            click.echo(f"Updating {source}...")
            
            try:
                if source == 'STRING_aliases':
                    ingestion.ingest_string_aliases(file_path)
                elif source in ['STRING_full', 'STRING_physical']:
                    ingestion.ingest_string_interactions(file_path, source)
                elif source == 'MitoCarta':
                    ingestion.ingest_mitocarta(file_path)
                elif source == 'HPA_muscle':
                    ingestion.ingest_hpa_muscle(file_path)
                    
                click.echo(f"✅ {source} updated successfully")
            except Exception as e:
                click.echo(f"❌ Failed to update {source}: {e}")
        else:
            click.echo(f"No updates needed for {source}")
    else:
        # Update all sources
        click.echo("Checking all sources for updates...")
        updated_sources = ingestion.ingest_all_sources(force_update=force)
        
        if updated_sources:
            click.echo(f"✅ Updated {len(updated_sources)} sources:")
            for source_name in updated_sources:
                click.echo(f"  - {source_name}")
        else:
            click.echo("No updates were needed")

@cli.command()
@click.option('--filter-type', type=click.Choice(['mitochondrial', 'muscle', 'genes', 'high_confidence', 'all']), 
              default='all', help='Type of network filter to apply')
@click.option('--genes', help='Comma-separated list of gene symbols (for genes filter)')
@click.option('--uniprots', help='Comma-separated list of UniProt IDs (for genes filter)')
@click.option('--neighbors', type=int, default=0, help='Include N-hop neighbors (0=direct only)')
@click.option('--min-confidence', type=float, default=0.4, help='Minimum confidence score')
@click.option('--max-confidence', type=float, help='Maximum confidence score')
@click.option('--evidence-types', help='Comma-separated evidence types (experimental,database,textmining)')
@click.option('--min-degree', type=int, help='Minimum protein degree (number of connections)')
@click.option('--max-degree', type=int, help='Maximum protein degree (number of connections)')
@click.option('--formats', default='json,graphml,csv', help='Output formats (json,graphml,csv)')
@click.option('--output-prefix', default='filtered_network', help='Output filename prefix')
@click.option('--output-dir', default='outputs', help='Output directory')
@click.pass_context
def export_network(ctx, filter_type, genes, uniprots, neighbors, min_confidence, max_confidence,
                  evidence_types, min_degree, max_degree, formats, output_prefix, output_dir):
    """Export filtered networks from the complete database"""
    db = ctx.obj['db']
    
    # Check if database has data
    stats = db.get_statistics()
    if stats['num_proteins'] == 0:
        click.echo("❌ Database is empty. Run 'update' to ingest data first.")
        return
    
    click.echo(f"Database contains {stats['num_proteins']:,} proteins and {stats['num_interactions']:,} interactions")
    
    # Create network filter based on type
    if filter_type == 'mitochondrial':
        network_filter = NetworkFilter.mitochondrial_network(min_confidence=min_confidence)
        click.echo("🧬 Creating mitochondrial protein network...")
        
    elif filter_type == 'muscle':
        network_filter = NetworkFilter.muscle_network(min_confidence=min_confidence)
        click.echo("💪 Creating muscle-expressed protein network...")
        
    elif filter_type == 'genes':
        if not genes and not uniprots:
            click.echo("❌ For 'genes' filter, provide --genes or --uniprots")
            return
        
        gene_list = []
        if genes:
            gene_list.extend([g.strip() for g in genes.split(',')])
        if uniprots:
            gene_list.extend([u.strip() for u in uniprots.split(',')])
            
        network_filter = NetworkFilter.gene_set_network(
            genes=gene_list, 
            include_neighbors=neighbors,
            min_confidence=min_confidence
        )
        click.echo(f"🎯 Creating network for {len(gene_list)} genes with {neighbors}-hop neighbors...")
        
    elif filter_type == 'high_confidence':
        network_filter = NetworkFilter.high_confidence_network(min_confidence=min_confidence)
        click.echo(f"⭐ Creating high-confidence network (min_confidence={min_confidence})...")
        
    else:  # filter_type == 'all'
        network_filter = NetworkFilter()
        network_filter.min_confidence = min_confidence
        click.echo("🌐 Creating complete network...")
    
    # Apply additional filters
    if max_confidence:
        network_filter.max_confidence = max_confidence
    if evidence_types:
        network_filter.evidence_types = set(evidence_types.split(','))
    if min_degree:
        network_filter.min_degree = min_degree
    if max_degree:
        network_filter.max_degree = max_degree
    
    # Setup exporter
    exporter = NetworkExporter(db, Path(output_dir))
    format_list = [f.strip() for f in formats.split(',')]
    
    try:
        # Export network
        output_files = exporter.export_network(
            network_filter=network_filter,
            format_types=format_list,
            filename_prefix=output_prefix
        )
        
        click.echo("✅ Network export completed!")
        click.echo("\nOutput files:")
        for format_type, file_path in output_files.items():
            click.echo(f"  📄 {format_type}: {file_path}")
            
    except Exception as e:
        click.echo(f"❌ Export failed: {e}")

@cli.command()
@click.option('--output-dir', default='outputs', help='Output directory')
@click.option('--min-confidence', type=float, default=0.4, help='Minimum confidence score')
@click.pass_context
def export_predefined(ctx, output_dir, min_confidence):
    """Export commonly used predefined networks (mitochondrial, muscle, high-confidence)"""
    db = ctx.obj['db']
    
    # Check if database has data
    stats = db.get_statistics()
    if stats['num_proteins'] == 0:
        click.echo("❌ Database is empty. Run 'update' to ingest data first.")
        return
    
    click.echo(f"Exporting predefined networks from database with {stats['num_proteins']:,} proteins...")
    
    try:
        # Update the global min_confidence for predefined networks
        from .export import NetworkFilter
        NetworkFilter.mitochondrial_network.__func__.__defaults__ = (min_confidence,)
        NetworkFilter.muscle_network.__func__.__defaults__ = (min_confidence,)
        
        networks_exported = export_predefined_networks(db, Path(output_dir))
        
        click.echo("✅ Predefined networks exported!")
        click.echo("\nNetworks created:")
        
        for network_type, files in networks_exported.items():
            click.echo(f"\n🔸 {network_type.upper()} network:")
            for format_type, file_path in files.items():
                click.echo(f"    📄 {format_type}: {file_path}")
                
    except Exception as e:
        click.echo(f"❌ Export failed: {e}")

@cli.command()
@click.option('--phase', help='Specific processing phase to check')
@click.pass_context
def checkpoints(ctx, phase):
    """List processing checkpoints"""
    db = ctx.obj['db']
    session = db.get_session()
    
    try:
        from .database import ProcessingCheckpoint
        query = session.query(ProcessingCheckpoint)
        
        if phase:
            query = query.filter_by(phase=phase)
            
        checkpoints = query.order_by(ProcessingCheckpoint.created_at.desc()).limit(20).all()
        
        if not checkpoints:
            click.echo("No checkpoints found")
            return
            
        click.echo("Recent Processing Checkpoints:")
        click.echo("-" * 60)
        
        for cp in checkpoints:
            status_icon = "✅" if cp.status == "completed" else "❌" if cp.status == "failed" else "🔄"
            click.echo(f"{status_icon} {cp.checkpoint_name}")
            click.echo(f"   Phase: {cp.phase}")
            click.echo(f"   Status: {cp.status}")
            click.echo(f"   Created: {cp.created_at}")
            if cp.completed_at:
                click.echo(f"   Completed: {cp.completed_at}")
            click.echo()
            
    finally:
        session.close()


if __name__ == '__main__':
    cli()