"""
Command-line interface for MitoNet incremental updates
"""

import click
import logging
from pathlib import Path
from .database import MitoNetDatabase
from .ingestion import DataIngestionManager

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
@click.option('--genes', help='Comma-separated list of gene symbols to add')
@click.option('--uniprots', help='Comma-separated list of UniProt IDs to add')
@click.pass_context
def add_genes(ctx, genes, uniprots):
    """Add new genes to the database"""
    db = ctx.obj['db']
    
    added_count = 0
    
    if genes:
        gene_list = [g.strip() for g in genes.split(',')]
        click.echo(f"Adding {len(gene_list)} genes by symbol...")
        
        for gene_symbol in gene_list:
            # Try to find existing protein
            existing = db.find_protein_by_alias(gene_symbol, 'symbol')
            if existing:
                click.echo(f"  - {gene_symbol}: already exists ({existing.uniprot_id})")
            else:
                # Create placeholder protein (UniProt ID will be resolved later)
                placeholder_uniprot = f"TEMP_{gene_symbol}_{added_count}"
                protein = db.get_or_create_protein(
                    uniprot_id=placeholder_uniprot,
                    gene_symbol=gene_symbol
                )
                
                # Add alias in the same transaction
                from .database import ProteinAlias
                session = db.get_session()
                try:
                    # Refresh protein in this session
                    session.add(protein)
                    alias = ProteinAlias(
                        protein_id=protein.id,
                        alias_type='symbol',
                        alias_value=gene_symbol,
                        source_id=None
                    )
                    session.add(alias)
                    session.commit()
                    click.echo(f"  + {gene_symbol}: added as {placeholder_uniprot}")
                    added_count += 1
                except Exception as e:
                    session.rollback()
                    click.echo(f"  - {gene_symbol}: failed to add ({e})")
                finally:
                    session.close()
    
    if uniprots:
        uniprot_list = [u.strip() for u in uniprots.split(',')]
        click.echo(f"Adding {len(uniprot_list)} proteins by UniProt ID...")
        
        for uniprot_id in uniprot_list:
            existing = db.get_protein_by_uniprot(uniprot_id)
            if existing:
                click.echo(f"  - {uniprot_id}: already exists")
            else:
                protein = db.get_or_create_protein(uniprot_id=uniprot_id)
                click.echo(f"  + {uniprot_id}: added")
                added_count += 1
    
    click.echo(f"\n✅ Added {added_count} new proteins")

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

@cli.command()
@click.option('--output', '-o', default='mitonet_current.graphml', help='Output file path')
@click.pass_context
def export_network(ctx, output):
    """Export current network from database"""
    # This would integrate with the existing network export functionality
    click.echo(f"Exporting network to {output}...")
    click.echo("Note: Network export functionality to be implemented")

if __name__ == '__main__':
    cli()