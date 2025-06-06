"""
Incremental data ingestion system for MitoNet
"""

import hashlib
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Iterator, Tuple, Any
import pandas as pd
from .database import MitoNetDatabase, DataSource, Protein, Interaction

logger = logging.getLogger(__name__)

class DataIngestionManager:
    """Manages incremental data ingestion with change detection"""
    
    def __init__(self, db: MitoNetDatabase, data_dir: Path = Path("networks")):
        self.db = db
        self.data_dir = data_dir
        
    def calculate_file_hash(self, file_path: Path) -> str:
        """Calculate SHA256 hash of a file"""
        sha256_hash = hashlib.sha256()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                sha256_hash.update(chunk)
        return sha256_hash.hexdigest()
    
    def needs_update(self, source_name: str, file_path: Path, 
                    version: Optional[str] = None) -> bool:
        """Check if a data source needs updating"""
        if not file_path.exists():
            logger.warning(f"File not found: {file_path}")
            return False
            
        current_hash = self.calculate_file_hash(file_path)
        current_size = file_path.stat().st_size
        
        # Auto-detect version if not provided
        if not version:
            version = self._extract_version_from_filename(file_path)
        
        session = self.db.get_session()
        try:
            existing_source = session.query(DataSource).filter_by(
                name=source_name, version=version
            ).first()
            
            if not existing_source:
                logger.info(f"New data source detected: {source_name} v{version}")
                return True
                
            # Check if file has changed
            if (existing_source.file_hash != current_hash or 
                existing_source.file_size != current_size):
                logger.info(f"File changes detected for {source_name} v{version}")
                return True
                
            logger.debug(f"No changes for {source_name} v{version}")
            return False
            
        finally:
            session.close()
    
    def _extract_version_from_filename(self, file_path: Path) -> str:
        """Extract version from filename (heuristic)"""
        name = file_path.name
        
        # STRING patterns
        if 'v12.0' in name:
            return 'v12.0'
        elif 'v11.5' in name:
            return 'v11.5'
            
        # BioGRID patterns  
        if '4.4.246' in name:
            return '4.4.246'
        elif 'BIOGRID-ALL-' in name:
            # Extract version from BIOGRID-ALL-X.X.X format
            parts = name.split('-')
            for part in parts:
                if part.replace('.', '').isdigit():
                    return part
                    
        # MitoCarta patterns
        if 'MitoCarta3.0' in name:
            return '3.0'
        elif 'MitoCarta2.0' in name:
            return '2.0'
            
        # Default to file modification time as version
        mtime = datetime.fromtimestamp(file_path.stat().st_mtime)
        return mtime.strftime('%Y%m%d')
    
    def ingest_string_aliases(self, file_path: Path, version: Optional[str] = None,
                             chunk_size: int = 50000) -> DataSource:
        """Ingest STRING protein aliases incrementally"""
        if not version:
            version = self._extract_version_from_filename(file_path)
            
        logger.info(f"Ingesting STRING aliases v{version} from {file_path}")
        
        # Create or get data source
        source = self.db.get_or_create_data_source(
            name='STRING_aliases',
            version=version,
            file_path=str(file_path),
            file_size=file_path.stat().st_size,
            file_hash=self.calculate_file_hash(file_path)
        )
        
        # Process in chunks to avoid memory issues
        total_processed = 0
        aliases_added = 0
        
        chunk_reader = pd.read_csv(file_path, sep='\t', chunksize=chunk_size)
        
        for chunk_num, chunk in enumerate(chunk_reader, 1):
            logger.info(f"Processing chunk {chunk_num} ({len(chunk):,} rows)")
            
            for _, row in chunk.iterrows():
                string_id = row['#string_protein_id']
                alias = row['alias']
                alias_source = row['source']
                
                # Process UniProt mappings
                if alias_source == 'UniProt_AC':
                    protein = self.db.get_or_create_protein(uniprot_id=alias)
                    self.db.add_protein_alias(protein, 'string', string_id, source)
                    aliases_added += 1
                    
                # Process gene symbols
                elif 'UniProt_GN' in alias_source or 'BLAST_UniProt_GN' in alias_source:
                    # Find protein by UniProt ID via string mapping
                    existing_protein = self.db.find_protein_by_alias(string_id, 'string')
                    if existing_protein:
                        self.db.add_protein_alias(existing_protein, 'symbol', alias, source)
                        # Update gene symbol if not set
                        if not existing_protein.gene_symbol:
                            session = self.db.get_session()
                            try:
                                existing_protein.gene_symbol = alias
                                existing_protein.updated_at = datetime.utcnow()
                                session.commit()
                            finally:
                                session.close()
                        aliases_added += 1
                
                total_processed += 1
                
            # Save checkpoint every chunk
            self.db.save_checkpoint(
                name=f"string_aliases_v{version}_chunk_{chunk_num}",
                phase="alias_ingestion",
                data={"total_processed": total_processed, "aliases_added": aliases_added}
            )
        
        logger.info(f"STRING aliases ingestion complete: {total_processed:,} rows processed, "
                   f"{aliases_added:,} aliases added")
        
        return source
    
    def ingest_string_interactions(self, file_path: Path, source_name: str,
                                  version: Optional[str] = None,
                                  chunk_size: int = 50000) -> DataSource:
        """Ingest STRING protein interactions incrementally"""
        if not version:
            version = self._extract_version_from_filename(file_path)
            
        logger.info(f"Ingesting {source_name} v{version} from {file_path}")
        
        # Create or get data source
        source = self.db.get_or_create_data_source(
            name=source_name,
            version=version,
            file_path=str(file_path),
            file_size=file_path.stat().st_size,
            file_hash=self.calculate_file_hash(file_path)
        )
        
        total_processed = 0
        interactions_added = 0
        
        # Determine separator
        sep = ' ' if 'links' in file_path.name else '\t'
        
        chunk_reader = pd.read_csv(file_path, sep=sep, chunksize=chunk_size)
        
        for chunk_num, chunk in enumerate(chunk_reader, 1):
            logger.info(f"Processing chunk {chunk_num} ({len(chunk):,} rows)")
            
            for _, row in chunk.iterrows():
                protein1_string = row['protein1']
                protein2_string = row['protein2']
                confidence = row.get('combined_score', 0) / 1000.0  # Normalize to 0-1
                
                # Find proteins by STRING ID
                protein1 = self.db.find_protein_by_alias(protein1_string, 'string')
                protein2 = self.db.find_protein_by_alias(protein2_string, 'string')
                
                if protein1 and protein2 and protein1.id != protein2.id:
                    # Extract additional scores
                    source_scores = {}
                    for score_col in ['experimental', 'database', 'textmining', 
                                    'neighborhood', 'fusion', 'cooccurence', 'coexpression']:
                        if score_col in row:
                            source_scores[score_col] = row[score_col]
                    
                    # Determine evidence type
                    evidence_type = self._classify_string_evidence(row)
                    interaction_type = 'physical' if 'physical' in source_name else 'functional'
                    
                    self.db.add_interaction(
                        protein1=protein1,
                        protein2=protein2,
                        source=source,
                        confidence_score=confidence,
                        evidence_type=evidence_type,
                        interaction_type=interaction_type,
                        source_specific_id=f"{protein1_string}___{protein2_string}",
                        source_scores=source_scores
                    )
                    
                    interactions_added += 1
                
                total_processed += 1
            
            # Save checkpoint
            self.db.save_checkpoint(
                name=f"{source_name}_v{version}_chunk_{chunk_num}",
                phase="interaction_ingestion",
                data={"total_processed": total_processed, "interactions_added": interactions_added}
            )
        
        logger.info(f"{source_name} ingestion complete: {total_processed:,} rows processed, "
                   f"{interactions_added:,} interactions added")
        
        return source
    
    def ingest_mitocarta(self, file_path: Path, version: Optional[str] = None) -> DataSource:
        """Ingest MitoCarta data"""
        if not version:
            version = self._extract_version_from_filename(file_path)
            
        logger.info(f"Ingesting MitoCarta v{version} from {file_path}")
        
        source = self.db.get_or_create_data_source(
            name='MitoCarta',
            version=version,
            file_path=str(file_path),
            file_size=file_path.stat().st_size,
            file_hash=self.calculate_file_hash(file_path)
        )
        
        # Load MitoCarta data
        df = pd.read_excel(file_path, sheet_name='A Human MitoCarta3.0', engine='xlrd')
        
        proteins_updated = 0
        
        for _, row in df.iterrows():
            symbol = row['Symbol']
            
            # Find protein by gene symbol
            protein = self.db.find_protein_by_alias(symbol, 'symbol')
            
            if protein:
                # Update mitochondrial attributes
                session = self.db.get_session()
                try:
                    protein.is_mitochondrial = True
                    protein.mitocarta_list = row.get('MitoCarta3.0_List', '')
                    protein.mitocarta_evidence = row.get('MitoCarta3.0_Evidence', '')
                    protein.mitocarta_sub_localization = row.get('MitoCarta3.0_SubMitoLocalization', '')
                    protein.mitocarta_pathways = row.get('MitoCarta3.0_MitoPathways', '')
                    protein.gene_description = row.get('Description', protein.gene_description)
                    protein.updated_at = datetime.utcnow()
                    
                    session.commit()
                    proteins_updated += 1
                finally:
                    session.close()
        
        logger.info(f"MitoCarta ingestion complete: {proteins_updated:,} proteins updated")
        return source
    
    def ingest_hpa_muscle(self, file_path: Path, version: Optional[str] = None) -> DataSource:
        """Ingest HPA skeletal muscle data"""
        if not version:
            version = self._extract_version_from_filename(file_path)
            
        logger.info(f"Ingesting HPA muscle v{version} from {file_path}")
        
        source = self.db.get_or_create_data_source(
            name='HPA_muscle',
            version=version,
            file_path=str(file_path),
            file_size=file_path.stat().st_size,
            file_hash=self.calculate_file_hash(file_path)
        )
        
        df = pd.read_csv(file_path, sep='\t')
        
        proteins_updated = 0
        
        for _, row in df.iterrows():
            symbol = row['Gene']
            muscle_tpm = float(row.get('Tissue RNA - skeletal muscle [nTPM]', 0.0))
            
            # Only process muscle-expressed genes
            if muscle_tpm > 0:
                protein = self.db.find_protein_by_alias(symbol, 'symbol')
                
                if protein:
                    session = self.db.get_session()
                    try:
                        protein.is_muscle_expressed = True
                        protein.muscle_tpm = muscle_tpm
                        protein.protein_evidence_level = row.get('Evidence', '')
                        protein.main_localization = row.get('Subcellular main location', '')
                        protein.gene_description = row.get('Gene description', protein.gene_description)
                        protein.updated_at = datetime.utcnow()
                        
                        session.commit()
                        proteins_updated += 1
                    finally:
                        session.close()
        
        logger.info(f"HPA muscle ingestion complete: {proteins_updated:,} proteins updated")
        return source
    
    def _classify_string_evidence(self, row) -> str:
        """Classify STRING evidence type based on score distribution"""
        experimental = row.get('experimental', 0)
        database = row.get('database', 0)
        textmining = row.get('textmining', 0)
        
        if experimental > max(database, textmining):
            return 'experimental'
        elif database > textmining:
            return 'database'
        else:
            return 'text_mining'
    
    def ingest_all_sources(self, force_update: bool = False) -> Dict[str, DataSource]:
        """Ingest all available data sources"""
        sources = {}
        
        # Define source files
        source_files = [
            ('STRING_aliases', self.data_dir / 'string/9606.protein.aliases.v12.0.txt.gz'),
            ('STRING_info', self.data_dir / 'string/9606.protein.info.v12.0.txt.gz'),
            ('STRING_full', self.data_dir / 'string/9606.protein.links.detailed.v12.0.txt.gz'),
            ('STRING_physical', self.data_dir / 'string/9606.protein.physical.links.detailed.v12.0.txt.gz'),
            ('MitoCarta', self.data_dir / 'mitocarta/Human.MitoCarta3.0.xls'),
            ('HPA_muscle', self.data_dir / 'hpa/hpa_skm.tsv'),
        ]
        
        for source_name, file_path in source_files:
            if not file_path.exists():
                logger.warning(f"File not found: {file_path}")
                continue
                
            if force_update or self.needs_update(source_name, file_path):
                try:
                    if source_name == 'STRING_aliases':
                        sources[source_name] = self.ingest_string_aliases(file_path)
                    elif source_name in ['STRING_full', 'STRING_physical']:
                        sources[source_name] = self.ingest_string_interactions(file_path, source_name)
                    elif source_name == 'MitoCarta':
                        sources[source_name] = self.ingest_mitocarta(file_path)
                    elif source_name == 'HPA_muscle':
                        sources[source_name] = self.ingest_hpa_muscle(file_path)
                    else:
                        logger.info(f"Skipping {source_name} - no ingestion method defined")
                        
                except Exception as e:
                    logger.error(f"Failed to ingest {source_name}: {e}")
                    
        return sources