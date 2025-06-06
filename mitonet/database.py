"""
Database models and schema for persistent storage of MitoNet data
"""

from datetime import datetime
from typing import Dict, List, Optional, Any
from sqlalchemy import (
    create_engine, Column, Integer, String, Float, Boolean, DateTime, 
    Text, JSON, ForeignKey, Index, UniqueConstraint
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy.engine import Engine
import json
import logging

logger = logging.getLogger(__name__)

Base = declarative_base()

class DataSource(Base):
    """Track versions and metadata for each data source"""
    __tablename__ = 'data_sources'
    
    id = Column(Integer, primary_key=True)
    name = Column(String(50), nullable=False)  # 'STRING', 'BioGRID', etc.
    version = Column(String(20), nullable=False)  # 'v12.0', '4.4.246', etc.
    file_path = Column(String(500), nullable=False)
    file_size = Column(Integer)
    file_hash = Column(String(64))  # SHA256 for integrity
    last_updated = Column(DateTime, default=datetime.utcnow)
    source_metadata = Column(JSON)  # Store additional source-specific info
    
    __table_args__ = (
        UniqueConstraint('name', 'version', name='_source_version_uc'),
        Index('idx_source_name', 'name'),
    )

class Protein(Base):
    """Comprehensive protein information"""
    __tablename__ = 'proteins'
    
    id = Column(Integer, primary_key=True)
    uniprot_id = Column(String(20), unique=True, nullable=False, index=True)
    gene_symbol = Column(String(50), index=True)
    gene_description = Column(Text)
    
    # MitoCarta attributes
    is_mitochondrial = Column(Boolean, default=False, index=True)
    mitocarta_list = Column(String(20))
    mitocarta_evidence = Column(String(100))
    mitocarta_sub_localization = Column(String(100))
    mitocarta_pathways = Column(Text)
    
    # HPA attributes
    is_muscle_expressed = Column(Boolean, default=False, index=True)
    muscle_tpm = Column(Float, default=0.0)
    protein_evidence_level = Column(String(50))
    main_localization = Column(String(100))
    
    # Derived attributes
    priority_score = Column(Float, index=True)
    
    # Tracking
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Additional attributes as JSON for flexibility
    extra_attributes = Column(JSON)

class ProteinAlias(Base):
    """Store all protein identifiers and aliases"""
    __tablename__ = 'protein_aliases'
    
    id = Column(Integer, primary_key=True)
    protein_id = Column(Integer, ForeignKey('proteins.id'), nullable=False)
    alias_type = Column(String(20), nullable=False)  # 'string', 'ensembl', 'symbol', etc.
    alias_value = Column(String(100), nullable=False, index=True)
    source_id = Column(Integer, ForeignKey('data_sources.id'))
    
    protein = relationship("Protein", backref="aliases")
    source = relationship("DataSource")
    
    __table_args__ = (
        Index('idx_alias_value', 'alias_value'),
        Index('idx_alias_type_value', 'alias_type', 'alias_value'),
    )

class Interaction(Base):
    """Store protein-protein interactions"""
    __tablename__ = 'interactions'
    
    id = Column(Integer, primary_key=True)
    protein1_id = Column(Integer, ForeignKey('proteins.id'), nullable=False)
    protein2_id = Column(Integer, ForeignKey('proteins.id'), nullable=False)
    
    # Confidence and evidence
    confidence_score = Column(Float, index=True)
    evidence_type = Column(String(50))
    interaction_type = Column(String(50))  # 'physical', 'functional', etc.
    
    # Source tracking
    source_id = Column(Integer, ForeignKey('data_sources.id'), nullable=False)
    source_specific_id = Column(String(100))  # Original ID from source
    
    # Source-specific scores stored as JSON
    source_scores = Column(JSON)
    
    # Metadata
    created_at = Column(DateTime, default=datetime.utcnow)
    
    protein1 = relationship("Protein", foreign_keys=[protein1_id])
    protein2 = relationship("Protein", foreign_keys=[protein2_id])
    source = relationship("DataSource")
    
    __table_args__ = (
        # Ensure no duplicate interactions per source
        UniqueConstraint('protein1_id', 'protein2_id', 'source_id', name='_interaction_source_uc'),
        Index('idx_proteins', 'protein1_id', 'protein2_id'),
        Index('idx_confidence', 'confidence_score'),
    )

class NetworkVersion(Base):
    """Track different versions of the integrated network"""
    __tablename__ = 'network_versions'
    
    id = Column(Integer, primary_key=True)
    version = Column(String(20), nullable=False, unique=True)
    description = Column(Text)
    created_at = Column(DateTime, default=datetime.utcnow)
    
    # Network statistics
    num_nodes = Column(Integer)
    num_edges = Column(Integer)
    num_mitochondrial = Column(Integer)
    num_muscle_expressed = Column(Integer)
    
    # Processing parameters used
    parameters = Column(JSON)
    
    # Data sources used (JSON list of source IDs)
    data_sources_used = Column(JSON)

class ProcessingCheckpoint(Base):
    """Store checkpoints for resuming processing"""
    __tablename__ = 'processing_checkpoints'
    
    id = Column(Integer, primary_key=True)
    checkpoint_name = Column(String(100), nullable=False)
    phase = Column(String(50), nullable=False)  # 'id_mapping', 'network_integration', etc.
    status = Column(String(20), default='in_progress')  # 'completed', 'failed', 'in_progress'
    created_at = Column(DateTime, default=datetime.utcnow)
    completed_at = Column(DateTime)
    
    # Store intermediate results
    data = Column(JSON)
    error_message = Column(Text)

class MitoNetDatabase:
    """Database manager for MitoNet data"""
    
    def __init__(self, db_path: str = "mitonet.db"):
        self.db_path = db_path
        self.engine = create_engine(f"sqlite:///{db_path}")
        self.SessionLocal = sessionmaker(bind=self.engine)
        
    def initialize_db(self):
        """Create all tables"""
        Base.metadata.create_all(bind=self.engine)
        logger.info(f"Database initialized at {self.db_path}")
    
    def get_session(self):
        """Get a database session"""
        return self.SessionLocal()
    
    def get_or_create_data_source(self, name: str, version: str, 
                                  file_path: str, **metadata) -> DataSource:
        """Get existing data source or create new one"""
        session = self.get_session()
        try:
            source = session.query(DataSource).filter_by(
                name=name, version=version
            ).first()
            
            if not source:
                source = DataSource(
                    name=name,
                    version=version,
                    file_path=file_path,
                    source_metadata=metadata
                )
                session.add(source)
                session.commit()
                logger.info(f"Created new data source: {name} v{version}")
            
            return source
        finally:
            session.close()
    
    def get_protein_by_uniprot(self, uniprot_id: str) -> Optional[Protein]:
        """Get protein by UniProt ID"""
        session = self.get_session()
        try:
            return session.query(Protein).filter_by(uniprot_id=uniprot_id).first()
        finally:
            session.close()
    
    def get_or_create_protein(self, uniprot_id: str, **attributes) -> Protein:
        """Get existing protein or create new one"""
        session = self.get_session()
        try:
            protein = session.query(Protein).filter_by(uniprot_id=uniprot_id).first()
            
            if not protein:
                protein = Protein(uniprot_id=uniprot_id, **attributes)
                session.add(protein)
                session.commit()
                logger.debug(f"Created new protein: {uniprot_id}")
            else:
                # Update existing protein with new attributes
                for key, value in attributes.items():
                    if hasattr(protein, key) and value is not None:
                        setattr(protein, key, value)
                protein.updated_at = datetime.utcnow()
                session.commit()
                
            return protein
        finally:
            session.close()
    
    def add_protein_alias(self, protein: Protein, alias_type: str, 
                         alias_value: str, source: Optional[DataSource]):
        """Add a protein alias"""
        session = self.get_session()
        try:
            # Merge protein into current session
            protein = session.merge(protein)
            
            # Check if alias already exists
            existing = session.query(ProteinAlias).filter_by(
                protein_id=protein.id,
                alias_type=alias_type,
                alias_value=alias_value
            ).first()
            
            if not existing:
                alias = ProteinAlias(
                    protein_id=protein.id,
                    alias_type=alias_type,
                    alias_value=alias_value,
                    source_id=source.id if source else None
                )
                session.add(alias)
                session.commit()
        finally:
            session.close()
    
    def find_protein_by_alias(self, alias_value: str, 
                             alias_type: Optional[str] = None) -> Optional[Protein]:
        """Find protein by any alias"""
        session = self.get_session()
        try:
            query = session.query(Protein).join(ProteinAlias).filter(
                ProteinAlias.alias_value == alias_value
            )
            
            if alias_type:
                query = query.filter(ProteinAlias.alias_type == alias_type)
                
            return query.first()
        finally:
            session.close()
    
    def add_interaction(self, protein1: Protein, protein2: Protein,
                       source: DataSource, confidence_score: float,
                       **attributes) -> Interaction:
        """Add or update an interaction"""
        session = self.get_session()
        try:
            # Ensure consistent ordering (smaller ID first)
            if protein1.id > protein2.id:
                protein1, protein2 = protein2, protein1
            
            # Check if interaction already exists for this source
            existing = session.query(Interaction).filter_by(
                protein1_id=protein1.id,
                protein2_id=protein2.id,
                source_id=source.id
            ).first()
            
            if existing:
                # Update existing interaction
                existing.confidence_score = max(existing.confidence_score or 0, confidence_score)
                for key, value in attributes.items():
                    if hasattr(existing, key):
                        setattr(existing, key, value)
                session.commit()
                return existing
            else:
                # Create new interaction
                interaction = Interaction(
                    protein1_id=protein1.id,
                    protein2_id=protein2.id,
                    source_id=source.id,
                    confidence_score=confidence_score,
                    **attributes
                )
                session.add(interaction)
                session.commit()
                return interaction
        finally:
            session.close()
    
    def save_checkpoint(self, name: str, phase: str, data: Dict[str, Any], 
                       status: str = 'completed'):
        """Save a processing checkpoint"""
        session = self.get_session()
        try:
            checkpoint = ProcessingCheckpoint(
                checkpoint_name=name,
                phase=phase,
                status=status,
                data=data,
                completed_at=datetime.utcnow() if status == 'completed' else None
            )
            session.add(checkpoint)
            session.commit()
            logger.info(f"Saved checkpoint: {name} ({phase})")
        finally:
            session.close()
    
    def load_checkpoint(self, name: str) -> Optional[Dict[str, Any]]:
        """Load a processing checkpoint"""
        session = self.get_session()
        try:
            checkpoint = session.query(ProcessingCheckpoint).filter_by(
                checkpoint_name=name,
                status='completed'
            ).order_by(ProcessingCheckpoint.created_at.desc()).first()
            
            return checkpoint.data if checkpoint else None
        finally:
            session.close()
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get database statistics"""
        session = self.get_session()
        try:
            stats = {
                'num_proteins': session.query(Protein).count(),
                'num_mitochondrial': session.query(Protein).filter_by(is_mitochondrial=True).count(),
                'num_muscle_expressed': session.query(Protein).filter_by(is_muscle_expressed=True).count(),
                'num_interactions': session.query(Interaction).count(),
                'num_sources': session.query(DataSource).count(),
                'num_aliases': session.query(ProteinAlias).count(),
            }
            
            # Source breakdown
            sources = session.query(DataSource.name, DataSource.version).all()
            stats['data_sources'] = [f"{name} v{version}" for name, version in sources]
            
            return stats
        finally:
            session.close()