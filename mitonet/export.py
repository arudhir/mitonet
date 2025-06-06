"""
Network export and filtering functionality for MitoNet
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Any, Union
import pandas as pd
import networkx as nx
from .database import MitoNetDatabase, Protein, Interaction

logger = logging.getLogger(__name__)

class NetworkFilter:
    """Defines filtering criteria for network export"""
    
    def __init__(self):
        self.include_mitochondrial: Optional[bool] = None
        self.include_muscle_expressed: Optional[bool] = None
        self.gene_symbols: Optional[Set[str]] = None
        self.uniprot_ids: Optional[Set[str]] = None
        self.min_confidence: Optional[float] = None
        self.max_confidence: Optional[float] = None
        self.evidence_types: Optional[Set[str]] = None
        self.interaction_types: Optional[Set[str]] = None
        self.data_sources: Optional[Set[str]] = None
        self.include_neighbors: int = 0  # 0 = direct only, 1 = 1-hop neighbors, etc.
        self.min_degree: Optional[int] = None
        self.max_degree: Optional[int] = None
        
    @classmethod
    def mitochondrial_network(cls, min_confidence: float = 0.4) -> 'NetworkFilter':
        """Create filter for mitochondrial proteins"""
        filter_obj = cls()
        filter_obj.include_mitochondrial = True
        filter_obj.min_confidence = min_confidence
        return filter_obj
    
    @classmethod
    def muscle_network(cls, min_confidence: float = 0.4) -> 'NetworkFilter':
        """Create filter for muscle-expressed proteins"""
        filter_obj = cls()
        filter_obj.include_muscle_expressed = True
        filter_obj.min_confidence = min_confidence
        return filter_obj
    
    @classmethod
    def gene_set_network(cls, genes: List[str], include_neighbors: int = 1,
                        min_confidence: float = 0.4) -> 'NetworkFilter':
        """Create filter for specific gene set with neighbors"""
        filter_obj = cls()
        filter_obj.gene_symbols = set(genes)
        filter_obj.include_neighbors = include_neighbors
        filter_obj.min_confidence = min_confidence
        return filter_obj
    
    @classmethod
    def high_confidence_network(cls, min_confidence: float = 0.7) -> 'NetworkFilter':
        """Create filter for high-confidence interactions only"""
        filter_obj = cls()
        filter_obj.min_confidence = min_confidence
        return filter_obj

class NetworkExporter:
    """Export filtered networks from the complete database"""
    
    def __init__(self, db: MitoNetDatabase, output_dir: Path = Path("outputs")):
        self.db = db
        self.output_dir = output_dir
        self.output_dir.mkdir(exist_ok=True)
        
    def export_network(self, network_filter: NetworkFilter, 
                      format_types: List[str] = ['json', 'graphml', 'csv'],
                      filename_prefix: str = 'filtered_network') -> Dict[str, Path]:
        """Export a filtered network in specified formats"""
        
        logger.info("Building filtered network from database...")
        
        # Get filtered proteins and interactions
        proteins = self._get_filtered_proteins(network_filter)
        interactions = self._get_filtered_interactions(network_filter, proteins)
        
        # Apply neighbor expansion if requested
        if network_filter.include_neighbors > 0:
            proteins, interactions = self._expand_neighbors(
                proteins, interactions, network_filter.include_neighbors
            )
        
        # Apply degree filters
        proteins, interactions = self._apply_degree_filters(
            proteins, interactions, network_filter
        )
        
        logger.info(f"Filtered network: {len(proteins)} proteins, {len(interactions)} interactions")
        
        # Build NetworkX graph
        graph = self._build_networkx_graph(proteins, interactions)
        
        # Export in requested formats
        output_files = {}
        
        if 'json' in format_types:
            json_file = self._export_json(graph, filename_prefix)
            output_files['json'] = json_file
            
        if 'graphml' in format_types:
            graphml_file = self._export_graphml(graph, filename_prefix)
            output_files['graphml'] = graphml_file
            
        if 'csv' in format_types:
            csv_files = self._export_csv(proteins, interactions, filename_prefix)
            output_files.update(csv_files)
            
        return output_files
    
    def _get_filtered_proteins(self, network_filter: NetworkFilter) -> List[Protein]:
        """Get proteins that match the filter criteria"""
        session = self.db.get_session()
        try:
            query = session.query(Protein)
            
            # Filter by mitochondrial status
            if network_filter.include_mitochondrial is not None:
                query = query.filter(Protein.is_mitochondrial == network_filter.include_mitochondrial)
            
            # Filter by muscle expression
            if network_filter.include_muscle_expressed is not None:
                query = query.filter(Protein.is_muscle_expressed == network_filter.include_muscle_expressed)
            
            # Filter by specific genes
            if network_filter.gene_symbols:
                query = query.filter(Protein.gene_symbol.in_(network_filter.gene_symbols))
            
            # Filter by specific UniProt IDs
            if network_filter.uniprot_ids:
                query = query.filter(Protein.uniprot_id.in_(network_filter.uniprot_ids))
            
            return query.all()
            
        finally:
            session.close()
    
    def _get_filtered_interactions(self, network_filter: NetworkFilter, 
                                  proteins: List[Protein]) -> List[Interaction]:
        """Get interactions that match the filter criteria"""
        session = self.db.get_session()
        try:
            protein_ids = {p.id for p in proteins}
            
            query = session.query(Interaction)
            
            # Only interactions between filtered proteins
            query = query.filter(
                Interaction.protein1_id.in_(protein_ids),
                Interaction.protein2_id.in_(protein_ids)
            )
            
            # Filter by confidence score
            if network_filter.min_confidence is not None:
                query = query.filter(Interaction.confidence_score >= network_filter.min_confidence)
            
            if network_filter.max_confidence is not None:
                query = query.filter(Interaction.confidence_score <= network_filter.max_confidence)
            
            # Filter by evidence types
            if network_filter.evidence_types:
                query = query.filter(Interaction.evidence_type.in_(network_filter.evidence_types))
            
            # Filter by interaction types
            if network_filter.interaction_types:
                query = query.filter(Interaction.interaction_type.in_(network_filter.interaction_types))
            
            # Filter by data sources
            if network_filter.data_sources:
                query = query.join(Interaction.source).filter(
                    Interaction.source.name.in_(network_filter.data_sources)
                )
            
            return query.all()
            
        finally:
            session.close()
    
    def _expand_neighbors(self, proteins: List[Protein], interactions: List[Interaction],
                         neighbor_levels: int) -> Tuple[List[Protein], List[Interaction]]:
        """Expand the network to include neighbors up to specified levels"""
        session = self.db.get_session()
        try:
            current_protein_ids = {p.id for p in proteins}
            all_interactions = list(interactions)
            
            for level in range(neighbor_levels):
                # Find interactions involving current proteins
                new_interactions = session.query(Interaction).filter(
                    (Interaction.protein1_id.in_(current_protein_ids)) |
                    (Interaction.protein2_id.in_(current_protein_ids))
                ).all()
                
                # Collect new protein IDs
                new_protein_ids = set()
                for interaction in new_interactions:
                    if interaction not in all_interactions:
                        all_interactions.append(interaction)
                        new_protein_ids.add(interaction.protein1_id)
                        new_protein_ids.add(interaction.protein2_id)
                
                # Add new proteins
                if new_protein_ids - current_protein_ids:
                    new_proteins = session.query(Protein).filter(
                        Protein.id.in_(new_protein_ids - current_protein_ids)
                    ).all()
                    proteins.extend(new_proteins)
                    current_protein_ids.update(new_protein_ids)
                else:
                    break  # No new proteins found
            
            return proteins, all_interactions
            
        finally:
            session.close()
    
    def _apply_degree_filters(self, proteins: List[Protein], interactions: List[Interaction],
                             network_filter: NetworkFilter) -> Tuple[List[Protein], List[Interaction]]:
        """Apply minimum/maximum degree filters"""
        if network_filter.min_degree is None and network_filter.max_degree is None:
            return proteins, interactions
        
        # Calculate degrees
        degree_count = {}
        for interaction in interactions:
            degree_count[interaction.protein1_id] = degree_count.get(interaction.protein1_id, 0) + 1
            degree_count[interaction.protein2_id] = degree_count.get(interaction.protein2_id, 0) + 1
        
        # Filter proteins by degree
        filtered_proteins = []
        for protein in proteins:
            degree = degree_count.get(protein.id, 0)
            
            if network_filter.min_degree is not None and degree < network_filter.min_degree:
                continue
            if network_filter.max_degree is not None and degree > network_filter.max_degree:
                continue
                
            filtered_proteins.append(protein)
        
        # Filter interactions to only include remaining proteins
        remaining_protein_ids = {p.id for p in filtered_proteins}
        filtered_interactions = [
            interaction for interaction in interactions
            if (interaction.protein1_id in remaining_protein_ids and 
                interaction.protein2_id in remaining_protein_ids)
        ]
        
        return filtered_proteins, filtered_interactions
    
    def _build_networkx_graph(self, proteins: List[Protein], 
                             interactions: List[Interaction]) -> nx.Graph:
        """Build a NetworkX graph from proteins and interactions"""
        graph = nx.Graph()
        
        # Add nodes with attributes
        for protein in proteins:
            graph.add_node(protein.uniprot_id, 
                          gene_symbol=protein.gene_symbol or '',
                          gene_description=protein.gene_description or '',
                          is_mitochondrial=protein.is_mitochondrial or False,
                          is_muscle_expressed=protein.is_muscle_expressed or False,
                          muscle_tpm=protein.muscle_tpm or 0.0,
                          priority_score=protein.priority_score or 0.0,
                          protein_evidence_level=protein.protein_evidence_level or '',
                          mitocarta_sub_localization=protein.mitocarta_sub_localization or '',
                          main_localization=protein.main_localization or '')
        
        # Add edges with attributes
        protein_id_to_uniprot = {p.id: p.uniprot_id for p in proteins}
        
        for interaction in interactions:
            protein1_uniprot = protein_id_to_uniprot[interaction.protein1_id]
            protein2_uniprot = protein_id_to_uniprot[interaction.protein2_id]
            
            graph.add_edge(protein1_uniprot, protein2_uniprot,
                          confidence_score=interaction.confidence_score or 0.0,
                          evidence_type=interaction.evidence_type or '',
                          interaction_type=interaction.interaction_type or '',
                          source_scores=interaction.source_scores or {})
        
        return graph
    
    def _export_json(self, graph: nx.Graph, filename_prefix: str) -> Path:
        """Export graph as JSON"""
        output_file = self.output_dir / f"{filename_prefix}.json"
        
        # Convert to JSON-serializable format
        data = nx.node_link_data(graph, edges="links")
        
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2, default=str)
        
        logger.info(f"Exported JSON network: {output_file}")
        return output_file
    
    def _export_graphml(self, graph: nx.Graph, filename_prefix: str) -> Path:
        """Export graph as GraphML"""
        output_file = self.output_dir / f"{filename_prefix}.graphml"
        
        # Create a copy of the graph with string-converted complex data types
        graphml_graph = graph.copy()
        
        # Convert dictionary attributes to strings for GraphML compatibility
        for node, data in graphml_graph.nodes(data=True):
            for key, value in data.items():
                if isinstance(value, dict):
                    data[key] = str(value)
        
        for u, v, data in graphml_graph.edges(data=True):
            for key, value in data.items():
                if isinstance(value, dict):
                    data[key] = str(value)
        
        nx.write_graphml(graphml_graph, output_file)
        logger.info(f"Exported GraphML network: {output_file}")
        return output_file
    
    def _export_csv(self, proteins: List[Protein], interactions: List[Interaction],
                   filename_prefix: str) -> Dict[str, Path]:
        """Export nodes and edges as CSV files"""
        
        # Export nodes
        nodes_file = self.output_dir / f"{filename_prefix}_nodes.csv"
        nodes_data = []
        
        for protein in proteins:
            nodes_data.append({
                'uniprot_id': protein.uniprot_id,
                'gene_symbol': protein.gene_symbol or '',
                'gene_description': protein.gene_description or '',
                'is_mitochondrial': protein.is_mitochondrial or False,
                'is_muscle_expressed': protein.is_muscle_expressed or False,
                'muscle_tpm': protein.muscle_tpm or 0.0,
                'priority_score': protein.priority_score or 0.0,
                'protein_evidence_level': protein.protein_evidence_level or '',
                'mitocarta_sub_localization': protein.mitocarta_sub_localization or '',
                'main_localization': protein.main_localization or ''
            })
        
        pd.DataFrame(nodes_data).to_csv(nodes_file, index=False)
        
        # Export edges
        edges_file = self.output_dir / f"{filename_prefix}_edges.csv"
        edges_data = []
        
        protein_id_to_uniprot = {p.id: p.uniprot_id for p in proteins}
        
        for interaction in interactions:
            edges_data.append({
                'protein1': protein_id_to_uniprot[interaction.protein1_id],
                'protein2': protein_id_to_uniprot[interaction.protein2_id],
                'confidence_score': interaction.confidence_score or 0.0,
                'evidence_type': interaction.evidence_type or '',
                'interaction_type': interaction.interaction_type or '',
                'source_scores': str(interaction.source_scores or {})
            })
        
        pd.DataFrame(edges_data).to_csv(edges_file, index=False)
        
        logger.info(f"Exported CSV files: {nodes_file}, {edges_file}")
        return {'nodes_csv': nodes_file, 'edges_csv': edges_file}

def export_predefined_networks(db: MitoNetDatabase, output_dir: Path = Path("outputs")) -> Dict[str, Any]:
    """Export commonly used predefined networks"""
    exporter = NetworkExporter(db, output_dir)
    
    networks_exported = {}
    
    # Mitochondrial network
    mito_filter = NetworkFilter.mitochondrial_network(min_confidence=0.4)
    mito_files = exporter.export_network(mito_filter, filename_prefix='mitochondrial_network')
    networks_exported['mitochondrial'] = mito_files
    
    # Muscle network
    muscle_filter = NetworkFilter.muscle_network(min_confidence=0.4)
    muscle_files = exporter.export_network(muscle_filter, filename_prefix='muscle_network')
    networks_exported['muscle'] = muscle_files
    
    # High confidence network (all proteins)
    high_conf_filter = NetworkFilter.high_confidence_network(min_confidence=0.7)
    high_conf_files = exporter.export_network(high_conf_filter, filename_prefix='high_confidence_network')
    networks_exported['high_confidence'] = high_conf_files
    
    return networks_exported