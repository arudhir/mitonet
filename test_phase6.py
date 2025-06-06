#!/usr/bin/env python3
"""
Test Phase 6 Quality Control with simulated network data
"""

import pandas as pd
import numpy as np
import networkx as nx
from pathlib import Path
import logging
import json
from typing import Dict, List, Set, Tuple, Optional

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class MockMitoNetIntegrator:
    """
    Mock version of MitoNetIntegrator for testing Phase 6
    """
    
    def __init__(self):
        self.data_dir = Path("networks")
        self.protein_reference = {}
        self.mitochondrial_uniprot_ids = set()
        self.muscle_expressed_uniprot_ids = set()
        self.included_uniprot_ids = set()
        self.id_mapping = {
            'string_to_uniprot': {},
            'symbol_to_uniprot': {},
            'uniprot_to_symbol': {},
            'uniprot_info': {}
        }
        self.network = nx.Graph()
        
        # Create mock data
        self._create_mock_network()
        
    def _create_mock_network(self):
        """Create a mock network with realistic structure"""
        logger.info("Creating mock network for Phase 6 testing...")
        
        # Create mock proteins
        mock_proteins = [
            ('P12345', 'ATP5F1A', True, True),  # Mito + Muscle
            ('Q67890', 'NDUFA1', True, False),  # Mito only
            ('R11111', 'ACTA1', False, True),   # Muscle only
            ('S22222', 'TPM1', False, True),    # Muscle only
            ('T33333', 'CYCS', True, True),     # Mito + Muscle
        ]
        
        for uniprot_id, symbol, is_mito, is_muscle in mock_proteins:
            # Add to reference
            self.protein_reference[uniprot_id] = {
                'gene_symbol': symbol,
                'gene_description': f'{symbol} protein description',
                'mitocarta_member': is_mito,
                'muscle_TPM': 100.0 if is_muscle else 0.0,
                'priority_score': 1.0 if (is_mito and is_muscle) else 0.5,
                'protein_evidence_level': 'Evidence at protein level',
                'main_localization': 'Mitochondria' if is_mito else 'Cytoplasm'
            }
            
            # Add to sets
            if is_mito:
                self.mitochondrial_uniprot_ids.add(uniprot_id)
            if is_muscle:
                self.muscle_expressed_uniprot_ids.add(uniprot_id)
            self.included_uniprot_ids.add(uniprot_id)
            
            # Add to ID mapping
            self.id_mapping['symbol_to_uniprot'][symbol] = uniprot_id
            self.id_mapping['uniprot_to_symbol'][uniprot_id] = symbol
            
            # Add to network with attributes
            self.network.add_node(uniprot_id, 
                                gene_symbol=symbol,
                                gene_description=f'{symbol} protein description',
                                mitocarta_member=is_mito,
                                muscle_expressed=is_muscle,
                                protein_evidence_level='Evidence at protein level',
                                main_localization='Mitochondria' if is_mito else 'Cytoplasm',
                                reactome_pathways=['Pathway1', 'Pathway2'] if is_mito else [],
                                protein_complexes=['Complex1'] if is_mito else [],
                                priority_score=1.0 if (is_mito and is_muscle) else 0.5)
        
        # Add edges with realistic attributes
        edges = [
            ('P12345', 'Q67890', 0.9, ['STRING_physical', 'BioGRID']),
            ('P12345', 'T33333', 0.8, ['STRING_full']),
            ('R11111', 'S22222', 0.7, ['STRING_physical']),
            ('T33333', 'Q67890', 0.85, ['CORUM', 'Reactome']),
        ]
        
        for u, v, confidence, sources in edges:
            self.network.add_edge(u, v,
                                sources=sources,
                                num_sources=len(sources),
                                composite_confidence=confidence,
                                max_confidence=confidence)
        
        logger.info(f"Created mock network: {self.network.number_of_nodes()} nodes, {self.network.number_of_edges()} edges")
    
    def map_to_uniprot(self, identifier: str, id_type: str = 'symbol') -> Optional[str]:
        """Mock ID mapping"""
        if id_type == 'symbol':
            return self.id_mapping['symbol_to_uniprot'].get(identifier)
        return None
    
    def is_mitochondrial(self, uniprot_id: str) -> bool:
        return uniprot_id in self.mitochondrial_uniprot_ids
        
    def is_muscle_expressed(self, uniprot_id: str) -> bool:
        return uniprot_id in self.muscle_expressed_uniprot_ids
        
    def get_protein_info(self, uniprot_id: str) -> Dict:
        return self.protein_reference.get(uniprot_id, {})

    # Include all the Phase 6 methods from the main class
    def _analyze_network_connectivity(self):
        """
        Analyze network connectivity and component structure
        """
        logger.info("\n" + "="*60)
        logger.info("NETWORK CONNECTIVITY ANALYSIS")
        logger.info("="*60)
        
        # Basic connectivity
        num_nodes = self.network.number_of_nodes()
        num_edges = self.network.number_of_edges()
        
        logger.info(f"Network size: {num_nodes:,} nodes, {num_edges:,} edges")
        
        # Connected components analysis
        components = list(nx.connected_components(self.network))
        num_components = len(components)
        
        logger.info(f"Connected components: {num_components}")
        
        # Analyze component sizes
        component_sizes = [len(comp) for comp in components]
        component_sizes.sort(reverse=True)
        
        logger.info(f"Largest component: {component_sizes[0]:,} nodes ({component_sizes[0]/num_nodes*100:.1f}%)")
        if len(component_sizes) > 1:
            logger.info(f"Second largest: {component_sizes[1]:,} nodes ({component_sizes[1]/num_nodes*100:.1f}%)")
            
        # Isolated nodes
        isolated_nodes = list(nx.isolates(self.network))
        logger.info(f"Isolated nodes: {len(isolated_nodes):,}")
        
        # Network density and efficiency
        density = nx.density(self.network)
        logger.info(f"Network density: {density:.6f}")
        
        logger.info("="*60)
        
    def _validate_mitochondrial_proteins(self):
        """
        Validate presence of key mitochondrial proteins
        """
        logger.info("\n" + "="*60)
        logger.info("MITOCHONDRIAL PROTEIN VALIDATION")
        logger.info("="*60)
        
        # Key mitochondrial proteins we expect to find
        key_proteins = {
            'OXPHOS': ['NDUFA1', 'ATP5F1A', 'CYCS'],
            'Test_proteins': ['ATP5F1A', 'NDUFA1']  # Use our mock proteins
        }
        
        validation_results = {}
        total_found = 0
        total_expected = 0
        
        for category, proteins in key_proteins.items():
            found = 0
            present_proteins = []
            
            for protein_symbol in proteins:
                # Find this protein in our network
                uniprot_id = self.map_to_uniprot(protein_symbol, 'symbol')
                if uniprot_id and uniprot_id in self.network.nodes():
                    found += 1
                    present_proteins.append(protein_symbol)
                    
            validation_results[category] = {
                'found': found,
                'total': len(proteins),
                'percentage': found / len(proteins) * 100,
                'present': present_proteins
            }
            
            total_found += found
            total_expected += len(proteins)
            
        # Print results
        for category, results in validation_results.items():
            logger.info(f"{category}: {results['found']}/{results['total']} ({results['percentage']:.1f}%)")
            
        overall_percentage = total_found / total_expected * 100
        logger.info(f"\nOverall validation: {total_found}/{total_expected} ({overall_percentage:.1f}%)")
        
        if overall_percentage >= 80:
            logger.info("✓ Good coverage of key mitochondrial proteins")
        elif overall_percentage >= 60:
            logger.info("⚠ Moderate coverage of key mitochondrial proteins")
        else:
            logger.info("✗ Low coverage of key mitochondrial proteins")
            
        logger.info("="*60)
        
    def _check_data_integrity(self):
        """
        Check data integrity and consistency
        """
        logger.info("\n" + "="*60)
        logger.info("DATA INTEGRITY CHECKS")
        logger.info("="*60)
        
        issues_found = []
        
        # Check 1: Verify all nodes have UniProt IDs
        nodes_without_uniprot = 0
        for node_id in self.network.nodes():
            if not node_id or len(node_id) < 6:  # UniProt IDs are typically 6+ chars
                nodes_without_uniprot += 1
                
        if nodes_without_uniprot == 0:
            logger.info("✓ All nodes have valid UniProt IDs")
        else:
            issues_found.append(f"Nodes with invalid UniProt IDs: {nodes_without_uniprot}")
            
        # Check 2: Verify protein reference consistency
        ref_proteins = len(self.protein_reference)
        network_proteins = self.network.number_of_nodes()
        
        logger.info(f"Protein reference: {ref_proteins:,} proteins")
        logger.info(f"Network nodes: {network_proteins:,} proteins")
        
        if abs(ref_proteins - network_proteins) / max(ref_proteins, network_proteins) < 0.05:
            logger.info("✓ Protein reference and network size consistent")
        else:
            issues_found.append(f"Size mismatch: reference={ref_proteins}, network={network_proteins}")
            
        # Check 3: Verify edge attributes
        edges_missing_sources = 0
        edges_missing_confidence = 0
        
        for u, v, data in self.network.edges(data=True):
            if 'sources' not in data or not data['sources']:
                edges_missing_sources += 1
            if 'composite_confidence' not in data:
                edges_missing_confidence += 1
                
        if edges_missing_sources == 0:
            logger.info("✓ All edges have source information")
        else:
            issues_found.append(f"Edges missing sources: {edges_missing_sources}")
            
        if edges_missing_confidence == 0:
            logger.info("✓ All edges have confidence scores")
        else:
            issues_found.append(f"Edges missing confidence: {edges_missing_confidence}")
            
        # Check 4: Validate mitochondrial/muscle classifications
        mito_count = sum(1 for n in self.network.nodes() if self.is_mitochondrial(n))
        muscle_count = sum(1 for n in self.network.nodes() if self.is_muscle_expressed(n))
        
        logger.info(f"Mitochondrial proteins in network: {mito_count:,}")
        logger.info(f"Muscle-expressed proteins in network: {muscle_count:,}")
        
        if mito_count > 0 and muscle_count > 0:
            logger.info("✓ Both mitochondrial and muscle proteins present")
        else:
            if mito_count == 0:
                issues_found.append("No mitochondrial proteins found")
            if muscle_count == 0:
                issues_found.append("No muscle-expressed proteins found")
                
        # Summary
        if not issues_found:
            logger.info("\n✓ All data integrity checks passed")
        else:
            logger.info(f"\n⚠ Issues found:")
            for issue in issues_found:
                logger.info(f"  - {issue}")
                
        logger.info("="*60)
        
    def _validate_node_annotations(self):
        """
        Validate quality and completeness of node annotations
        """
        logger.info("\n" + "="*60)
        logger.info("NODE ANNOTATION VALIDATION")
        logger.info("="*60)
        
        total_nodes = self.network.number_of_nodes()
        
        # Count annotations
        annotations_stats = {
            'gene_symbol': 0,
            'gene_description': 0,
            'protein_evidence_level': 0,
            'main_localization': 0,
            'reactome_pathways': 0,
            'protein_complexes': 0,
            'mitocarta_member': 0,
            'muscle_expressed': 0
        }
        
        for node_id in self.network.nodes():
            node_data = self.network.nodes[node_id]
            
            # Check each annotation type
            if node_data.get('gene_symbol'):
                annotations_stats['gene_symbol'] += 1
            if node_data.get('gene_description'):
                annotations_stats['gene_description'] += 1
            if node_data.get('protein_evidence_level'):
                annotations_stats['protein_evidence_level'] += 1
            if node_data.get('main_localization'):
                annotations_stats['main_localization'] += 1
            if node_data.get('reactome_pathways'):
                annotations_stats['reactome_pathways'] += 1
            if node_data.get('protein_complexes'):
                annotations_stats['protein_complexes'] += 1
            if node_data.get('mitocarta_member'):
                annotations_stats['mitocarta_member'] += 1
            if node_data.get('muscle_expressed'):
                annotations_stats['muscle_expressed'] += 1
                
        # Print results
        logger.info("Annotation coverage:")
        for annotation, count in annotations_stats.items():
            percentage = count / total_nodes * 100
            status = "✓" if percentage >= 80 else "⚠" if percentage >= 50 else "✗"
            logger.info(f"  {status} {annotation}: {count:,}/{total_nodes:,} ({percentage:.1f}%)")
            
        # Calculate overall annotation quality
        core_annotations = ['gene_symbol', 'gene_description', 'main_localization']
        core_coverage = sum(annotations_stats[ann] for ann in core_annotations) / (len(core_annotations) * total_nodes) * 100
        
        logger.info(f"\nCore annotation coverage: {core_coverage:.1f}%")
        if core_coverage >= 90:
            logger.info("✓ Excellent annotation quality")
        elif core_coverage >= 70:
            logger.info("⚠ Good annotation quality")
        else:
            logger.info("✗ Poor annotation quality")
            
        logger.info("="*60)
        
    def _compute_network_quality_metrics(self):
        """
        Compute comprehensive network quality metrics
        """
        logger.info("\n" + "="*60)
        logger.info("NETWORK QUALITY METRICS")
        logger.info("="*60)
        
        # Basic metrics
        num_nodes = self.network.number_of_nodes()
        num_edges = self.network.number_of_edges()
        
        # Multi-source validation
        multi_source_edges = sum(1 for _, _, data in self.network.edges(data=True) 
                               if data.get('num_sources', 0) > 1)
        multi_source_percentage = multi_source_edges / num_edges * 100 if num_edges > 0 else 0
        
        # Confidence distribution
        confidences = [data.get('composite_confidence', 0.0) 
                      for _, _, data in self.network.edges(data=True)]
        high_confidence_edges = sum(1 for c in confidences if c >= 0.8)
        medium_confidence_edges = sum(1 for c in confidences if 0.5 <= c < 0.8)
        low_confidence_edges = sum(1 for c in confidences if c < 0.5)
        
        # Network topology quality
        degrees = dict(self.network.degree())
        avg_degree = sum(degrees.values()) / len(degrees) if degrees else 0
        
        # Print metrics
        logger.info(f"Network size: {num_nodes:,} nodes, {num_edges:,} edges")
        logger.info(f"Average degree: {avg_degree:.2f}")
        logger.info(f"Multi-source edges: {multi_source_edges:,} ({multi_source_percentage:.1f}%)")
        
        logger.info(f"\nEdge confidence distribution:")
        logger.info(f"  High (≥0.8): {high_confidence_edges:,} ({high_confidence_edges/num_edges*100:.1f}%)")
        logger.info(f"  Medium (0.5-0.8): {medium_confidence_edges:,} ({medium_confidence_edges/num_edges*100:.1f}%)")
        logger.info(f"  Low (<0.5): {low_confidence_edges:,} ({low_confidence_edges/num_edges*100:.1f}%)")
        
        # Overall quality assessment
        quality_score = 0
        
        # Network size (10 points) - scale for small test network
        if num_nodes >= 3:
            quality_score += 10
        else:
            quality_score += 5
            
        # Multi-source validation (30 points)
        if multi_source_percentage >= 50:  # Lower threshold for test
            quality_score += 30
        elif multi_source_percentage >= 25:
            quality_score += 20
        else:
            quality_score += 10
            
        # High confidence edges (20 points)
        high_conf_percentage = high_confidence_edges / num_edges * 100 if num_edges > 0 else 0
        if high_conf_percentage >= 50:
            quality_score += 20
        elif high_conf_percentage >= 30:
            quality_score += 15
        else:
            quality_score += 10
            
        # Connectivity (20 points)
        components = list(nx.connected_components(self.network))
        largest_comp_size = max(len(comp) for comp in components) if components else 0
        connectivity_percentage = largest_comp_size / num_nodes * 100 if num_nodes > 0 else 0
        
        if connectivity_percentage >= 80:  # Lower threshold for test
            quality_score += 20
        elif connectivity_percentage >= 60:
            quality_score += 15
        else:
            quality_score += 10
            
        # Annotation completeness (20 points)
        annotated_nodes = sum(1 for n in self.network.nodes() 
                            if self.network.nodes[n].get('gene_symbol'))
        annotation_percentage = annotated_nodes / num_nodes * 100 if num_nodes > 0 else 0
        
        if annotation_percentage >= 95:
            quality_score += 20
        elif annotation_percentage >= 80:
            quality_score += 15
        else:
            quality_score += 10
            
        logger.info(f"\nOVERALL QUALITY SCORE: {quality_score}/100")
        
        if quality_score >= 90:
            logger.info("✓ EXCELLENT network quality")
        elif quality_score >= 75:
            logger.info("✓ GOOD network quality")
        elif quality_score >= 60:
            logger.info("⚠ FAIR network quality")
        else:
            logger.info("✗ POOR network quality")
            
        logger.info("="*60)
        
    def phase6_quality_control(self):
        """
        Phase 6: Comprehensive quality control and validation
        """
        logger.info("=== PHASE 6: QUALITY CONTROL AND VALIDATION ===")
        
        # Step 1: Network connectivity analysis
        logger.info("Analyzing network connectivity...")
        self._analyze_network_connectivity()
        
        # Step 2: Validate major mitochondrial proteins
        logger.info("Validating mitochondrial protein presence...")
        self._validate_mitochondrial_proteins()
        
        # Step 3: Check data integrity
        logger.info("Performing data integrity checks...")
        self._check_data_integrity()
        
        # Step 4: Validate node annotations
        logger.info("Validating node annotations...")
        self._validate_node_annotations()
        
        # Step 5: Assess network quality metrics
        logger.info("Computing network quality metrics...")
        self._compute_network_quality_metrics()
        
        logger.info("Phase 6 complete: Quality control finished")
        
    def phase7_generate_exports(self):
        """
        Phase 7: Generate network exports and comprehensive summary report
        """
        logger.info("=== PHASE 7: NETWORK EXPORTS AND REPORTING ===")
        
        # Step 1: Export network to GraphML for Cytoscape
        logger.info("Exporting network to GraphML format...")
        self._export_graphml()
        
        # Step 2: Export network to JSON for web visualization
        logger.info("Exporting network to JSON format...")
        self._export_json()
        
        # Step 3: Generate comprehensive summary report
        logger.info("Generating comprehensive summary report...")
        self._generate_summary_report()
        
        # Step 4: Export node and edge tables
        logger.info("Exporting node and edge data tables...")
        self._export_data_tables()
        
        logger.info("Phase 7 complete: All exports generated successfully")
        
    def _export_graphml(self):
        """Export network to GraphML format for Cytoscape"""
        output_file = "test_mitonet_network.graphml"
        try:
            # Clean and export network
            cleaned_network = self.network.copy()
            
            for node_id in cleaned_network.nodes():
                node_data = cleaned_network.nodes[node_id]
                for key, value in list(node_data.items()):
                    if isinstance(value, list):
                        node_data[key] = "|".join(str(v) for v in value)
                    elif isinstance(value, bool):
                        node_data[key] = str(value).lower()
                    elif value is None:
                        node_data[key] = ""
                        
            # Also clean edge attributes
            for u, v in cleaned_network.edges():
                edge_data = cleaned_network.edges[u, v]
                for key, value in list(edge_data.items()):
                    if isinstance(value, list):
                        edge_data[key] = "|".join(str(v) for v in value)
                    elif isinstance(value, bool):
                        edge_data[key] = str(value).lower()
                    elif value is None:
                        edge_data[key] = ""
                        
            nx.write_graphml(cleaned_network, output_file)
            file_size = Path(output_file).stat().st_size / 1024  # KB
            logger.info(f"✓ GraphML export: {output_file} ({file_size:.1f} KB)")
        except Exception as e:
            logger.error(f"✗ GraphML export failed: {e}")
            
    def _export_json(self):
        """Export network to JSON format for web visualization"""
        output_file = "test_mitonet_network.json"
        try:
            network_data = {
                'metadata': {
                    'version': '1.0',
                    'created': pd.Timestamp.now().isoformat(),
                    'description': 'Test mitochondrial and muscle protein network',
                    'nodes': self.network.number_of_nodes(),
                    'edges': self.network.number_of_edges()
                },
                'nodes': [],
                'edges': []
            }
            
            # Process nodes
            for node_id in self.network.nodes():
                node_data = self.network.nodes[node_id]
                clean_node = {
                    'id': node_id,
                    'label': node_data.get('gene_symbol', node_id),
                    'type': 'Mitochondrial' if self.is_mitochondrial(node_id) else 'Muscle',
                    'attributes': {
                        'gene_symbol': node_data.get('gene_symbol', ''),
                        'mitocarta_member': node_data.get('mitocarta_member', False),
                        'muscle_expressed': node_data.get('muscle_expressed', False)
                    }
                }
                network_data['nodes'].append(clean_node)
                
            # Process edges
            for u, v in self.network.edges():
                edge_data = self.network.edges[u, v]
                clean_edge = {
                    'source': u,
                    'target': v,
                    'confidence': edge_data.get('composite_confidence', 0.0),
                    'sources': edge_data.get('sources', [])
                }
                network_data['edges'].append(clean_edge)
                
            with open(output_file, 'w') as f:
                json.dump(network_data, f, indent=2, default=str)
                
            file_size = Path(output_file).stat().st_size / 1024  # KB
            logger.info(f"✓ JSON export: {output_file} ({file_size:.1f} KB)")
        except Exception as e:
            logger.error(f"✗ JSON export failed: {e}")
            
    def _export_data_tables(self):
        """Export node and edge data as CSV tables"""
        try:
            # Export nodes table
            nodes_data = []
            for node_id in self.network.nodes():
                node_data = self.network.nodes[node_id].copy()
                node_data['uniprot_id'] = node_id
                nodes_data.append(node_data)
                
            nodes_df = pd.DataFrame(nodes_data)
            nodes_file = "test_mitonet_nodes.csv"
            nodes_df.to_csv(nodes_file, index=False)
            logger.info(f"✓ Nodes table: {nodes_file} ({len(nodes_df):,} rows)")
            
            # Export edges table
            edges_data = []
            for u, v in self.network.edges():
                edge_data = self.network.edges[u, v].copy()
                edge_data['source'] = u
                edge_data['target'] = v
                edges_data.append(edge_data)
                
            edges_df = pd.DataFrame(edges_data)
            edges_file = "test_mitonet_edges.csv"
            edges_df.to_csv(edges_file, index=False)
            logger.info(f"✓ Edges table: {edges_file} ({len(edges_df):,} rows)")
            
        except Exception as e:
            logger.error(f"✗ Table export failed: {e}")
            
    def _generate_summary_report(self):
        """Generate comprehensive summary report"""
        report_file = "test_mitonet_summary_report.md"
        try:
            with open(report_file, 'w') as f:
                f.write("# Test Mitochondrial Network - Summary Report\\n\\n")
                f.write(f"**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\\n\\n")
                
                f.write("## Network Overview\\n\\n")
                f.write(f"- **Total nodes:** {self.network.number_of_nodes():,}\\n")
                f.write(f"- **Total edges:** {self.network.number_of_edges():,}\\n")
                
                degrees = dict(self.network.degree())
                avg_degree = sum(degrees.values()) / len(degrees) if degrees else 0
                f.write(f"- **Average degree:** {avg_degree:.2f}\\n")
                f.write(f"- **Network density:** {nx.density(self.network):.6f}\\n\\n")
                
                f.write("## Protein Composition\\n\\n")
                mito_count = sum(1 for n in self.network.nodes() if self.is_mitochondrial(n))
                muscle_count = sum(1 for n in self.network.nodes() if self.is_muscle_expressed(n))
                
                f.write(f"- **Mitochondrial proteins:** {mito_count:,}\\n")
                f.write(f"- **Muscle-expressed proteins:** {muscle_count:,}\\n\\n")
                
                f.write("## Generated Files\\n\\n")
                f.write("- `test_mitonet_network.graphml` - Network file for Cytoscape\\n")
                f.write("- `test_mitonet_network.json` - Network file for web visualization\\n")
                f.write("- `test_mitonet_nodes.csv` - Node attributes table\\n")
                f.write("- `test_mitonet_edges.csv` - Edge attributes table\\n")
                f.write("- `test_mitonet_summary_report.md` - This summary report\\n\\n")
                
            logger.info(f"✓ Summary report: {report_file}")
        except Exception as e:
            logger.error(f"✗ Report generation failed: {e}")

def main():
    """Test Phase 6 and 7 functionality"""
    logger.info("Testing Phase 6: Quality Control and Validation")
    logger.info("Testing Phase 7: Network Exports and Reporting")
    
    # Create mock integrator
    integrator = MockMitoNetIntegrator()
    
    # Execute Phase 6
    integrator.phase6_quality_control()
    
    # Execute Phase 7
    integrator.phase7_generate_exports()
    
    logger.info("Phase 6 and 7 testing complete")

if __name__ == "__main__":
    main()