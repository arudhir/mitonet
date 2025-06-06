#!/usr/bin/env python3
"""
Mitochondrial Network Integration Pipeline
Expert bioinformatics data integration for protein-protein interaction networks
"""

import pandas as pd
import numpy as np
import gzip
import zipfile
import json
import networkx as nx
from pathlib import Path
import logging
from typing import Dict, List, Set, Tuple, Optional
import re
import gc
import psutil
import os

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class MitoNetIntegrator:
    """
    Comprehensive biological network integration pipeline focusing on mitochondrial biology
    """
    
    def __init__(self, data_dir: str = "networks"):
        self.data_dir = Path(data_dir)
        self.inspection_report = []
        self.id_mapping = {}
        self.mitochondrial_proteins = set()
        self.network = nx.Graph()
        
        # Enable automatic memory monitoring
        self._monitor_memory("Initialization")
        
    def phase1_inspect_files(self):
        """
        Phase 1: Systematic inspection of all input files
        """
        logger.info("=== PHASE 1: DATA FILE INSPECTION ===")
        
        # Define expected files and their key columns (updated based on inspection)
        files_to_inspect = {
            'mitocarta': {
                'file': 'mitocarta/Human.MitoCarta3.0.xls',
                'expected_cols': ['Symbol', 'UniProt', 'MitoCarta3.0_SubMitoLocalization', 'MitoCarta3.0_MitoPathways']
            },
            'string_aliases': {
                'file': 'string/9606.protein.aliases.v12.0.txt.gz',
                'expected_cols': ['#string_protein_id', 'alias', 'source']
            },
            'string_links': {
                'file': 'string/9606.protein.links.detailed.v12.0.txt.gz',
                'expected_cols': ['protein1', 'protein2', 'neighborhood', 'fusion', 'cooccurence', 'coexpression', 'experimental', 'database', 'textmining', 'combined_score']
            },
            'string_physical': {
                'file': 'string/9606.protein.physical.links.detailed.v12.0.txt.gz',
                'expected_cols': ['protein1', 'protein2', 'experimental', 'database', 'textmining', 'combined_score']
            },
            'string_info': {
                'file': 'string/9606.protein.info.v12.0.txt.gz',
                'expected_cols': ['#string_protein_id', 'preferred_name', 'protein_size', 'annotation']
            },
            'biogrid': {
                'file': 'biogrid/BIOGRID-ALL-4.4.246.tab3.txt',
                'expected_cols': ['#BioGRID Interaction ID', 'Entrez Gene Interactor A', 'Entrez Gene Interactor B', 'Official Symbol Interactor A', 'Official Symbol Interactor B']
            },
            'reactome_interactions': {
                'file': 'reactome/reactome.homo_sapiens.interactions.tab-delimited.txt',
                'expected_cols': ['# Interactor 1 uniprot id', 'Interactor 2 uniprot id', 'Interaction type', 'Interaction context']
            },
            'reactome_uniprot': {
                'file': 'reactome/UniProt2Reactome_All_Levels.txt',
                'expected_cols': ['UniProt', 'Reactome_Pathway_ID', 'URL', 'Event_Name', 'Evidence_Code', 'Species']
            },
            'corum_complexes': {
                'file': 'corum/corum_humanComplexes.txt',
                'expected_cols': ['complex_id', 'complex_name', 'organism']  # subunits come from mapping file
            },
            'corum_mapping': {
                'file': 'corum/corum_uniprotCorumMapping.txt',
                'expected_cols': ['UniProtKB_accession_number', 'corum_id']
            },
            'hpa_skeletal_muscle': {
                'file': 'hpa/hpa_skm.tsv',
                'expected_cols': ['Gene', 'Gene description', 'Evidence', 'Tissue RNA - skeletal muscle [nTPM]', 'Subcellular main location']
            }
        }
        
        for source_name, file_info in files_to_inspect.items():
            self._inspect_single_file(source_name, file_info['file'], file_info['expected_cols'])
            
        # Print comprehensive inspection report
        self._print_inspection_report()
        self._monitor_memory("Phase 1 Complete")
        
    def _inspect_single_file(self, source_name: str, relative_path: str, expected_cols: List[str]):
        """
        Inspect a single file and generate detailed report
        """
        file_path = self.data_dir / relative_path
        
        report_entry = {
            'source': source_name,
            'file': str(file_path),
            'status': 'unknown',
            'actual_columns': [],
            'expected_columns': expected_cols,
            'sample_rows': [],
            'issues': [],
            'file_info': {}
        }
        
        try:
            # Check if file exists
            if not file_path.exists():
                report_entry['status'] = 'FILE_NOT_FOUND'
                report_entry['issues'].append(f"File not found: {file_path}")
                self.inspection_report.append(report_entry)
                return
                
            # Get file info
            report_entry['file_info'] = {
                'size_mb': round(file_path.stat().st_size / (1024*1024), 2),
                'extension': file_path.suffix
            }
            
            # Load file based on extension
            df = self._load_file_safely(file_path)
            
            if df is not None:
                report_entry['status'] = 'SUCCESS'
                report_entry['actual_columns'] = list(df.columns)
                report_entry['file_info']['rows'] = len(df)
                report_entry['file_info']['columns'] = len(df.columns)
                
                # Get sample rows (first 2 non-empty rows)
                sample_df = df.head(3).fillna('NaN')
                report_entry['sample_rows'] = sample_df.to_dict('records')
                
                # Check for expected columns
                missing_cols = set(expected_cols) - set(df.columns)
                if missing_cols:
                    report_entry['issues'].append(f"Missing expected columns: {missing_cols}")
                    
                # Check for completely empty columns
                empty_cols = [col for col in df.columns if df[col].isna().all()]
                if empty_cols:
                    report_entry['issues'].append(f"Completely empty columns: {empty_cols}")
                    
            else:
                report_entry['status'] = 'LOAD_FAILED'
                report_entry['issues'].append("Failed to load file - format not recognized")
                
        except Exception as e:
            report_entry['status'] = 'ERROR'
            report_entry['issues'].append(f"Error during inspection: {str(e)}")
            
        self.inspection_report.append(report_entry)
        
    def _load_file_safely(self, file_path: Path) -> Optional[pd.DataFrame]:
        """
        Safely load different file formats with appropriate handling
        """
        try:
            if file_path.suffix == '.gz':
                # Handle STRING files which are space-separated despite .txt extension
                if 'protein.links' in file_path.name or 'protein.physical' in file_path.name:
                    # STRING interaction files are space-separated
                    return pd.read_csv(file_path, sep=' ', low_memory=False, nrows=1000)
                elif 'tab' in file_path.stem or 'txt' in file_path.stem:
                    return pd.read_csv(file_path, sep='\t', low_memory=False, nrows=1000)
                else:
                    return pd.read_csv(file_path, low_memory=False, nrows=1000)
                    
            elif file_path.suffix in ['.txt', '.tsv']:
                # Special handling for specific files
                if 'UniProt2Reactome' in file_path.name:
                    # This file has no headers, assign column names
                    df = pd.read_csv(file_path, sep='\t', low_memory=False, nrows=1000, header=None)
                    df.columns = ['UniProt', 'Reactome_Pathway_ID', 'URL', 'Event_Name', 'Evidence_Code', 'Species']
                    return df
                else:
                    # Try tab-delimited first, then comma
                    try:
                        return pd.read_csv(file_path, sep='\t', low_memory=False, nrows=1000)
                    except:
                        return pd.read_csv(file_path, sep=',', low_memory=False, nrows=1000)
                    
            elif file_path.suffix == '.xls' or file_path.suffix == '.xlsx':
                # Special handling for MitoCarta Excel file
                if 'MitoCarta3.0' in file_path.name:
                    # Load specific sheet with actual data (sheet 'A Human MitoCarta3.0')
                    return pd.read_excel(file_path, sheet_name='A Human MitoCarta3.0', nrows=1000, engine='xlrd' if file_path.suffix == '.xls' else None)
                else:
                    return pd.read_excel(file_path, nrows=1000, engine='xlrd' if file_path.suffix == '.xls' else None)
                
            elif file_path.suffix == '.zip':
                # Handle BioGRID zip file
                with zipfile.ZipFile(file_path, 'r') as zip_ref:
                    # Find the main data file
                    txt_files = [f for f in zip_ref.namelist() if f.endswith('.txt') and 'tab3' in f]
                    if txt_files:
                        with zip_ref.open(txt_files[0]) as f:
                            return pd.read_csv(f, sep='\t', low_memory=False, nrows=1000)
                            
            return None
            
        except Exception as e:
            logger.error(f"Error loading {file_path}: {e}")
            return None
            
    def _print_inspection_report(self):
        """
        Print comprehensive inspection report
        """
        logger.info("\n" + "="*80)
        logger.info("DATA FILE INSPECTION REPORT")
        logger.info("="*80)
        
        for report in self.inspection_report:
            logger.info(f"\nDATA SOURCE: {report['source'].upper()}")
            logger.info(f"FILE: {report['file']}")
            logger.info(f"STATUS: {report['status']}")
            
            if report['file_info']:
                info = report['file_info']
                logger.info(f"SIZE: {info.get('size_mb', 'unknown')} MB")
                logger.info(f"DIMENSIONS: {info.get('rows', 'unknown')} rows × {info.get('columns', 'unknown')} columns")
                
            logger.info(f"COLUMNS FOUND: {report['actual_columns'][:10]}{'...' if len(report['actual_columns']) > 10 else ''}")
            logger.info(f"EXPECTED COLUMNS: {report['expected_columns']}")
            
            if report['issues']:
                logger.info(f"ISSUES: {'; '.join(report['issues'])}")
            else:
                logger.info("ISSUES: None")
                
            if report['sample_rows'] and len(report['sample_rows']) > 0:
                logger.info("SAMPLE ROWS:")
                for i, row in enumerate(report['sample_rows'][:2]):
                    logger.info(f"  Row {i+1}: {dict(list(row.items())[:3])}...")  # Show first 3 columns
                    
        logger.info("\n" + "="*80)
        
        # Summary
        success_count = sum(1 for r in self.inspection_report if r['status'] == 'SUCCESS')
        total_count = len(self.inspection_report)
        logger.info(f"INSPECTION SUMMARY: {success_count}/{total_count} files successfully loaded")
        
        failed_files = [r['source'] for r in self.inspection_report if r['status'] != 'SUCCESS']
        if failed_files:
            logger.info(f"FAILED FILES: {failed_files}")
            
        logger.info("="*80)
        
    def phase2_build_id_mapping(self):
        """
        Phase 2: Build comprehensive ID mapping system using UniProt as canonical identifier
        """
        logger.info("=== PHASE 2: IDENTIFIER STANDARDIZATION ===")
        
        # Initialize mapping dictionaries
        self.id_mapping = {
            'string_to_uniprot': {},
            'ensembl_to_uniprot': {},
            'entrez_to_uniprot': {},
            'symbol_to_uniprot': {},
            'uniprot_to_symbol': {},
            'uniprot_info': {}
        }
        
        # Step 1: Load STRING aliases for comprehensive mapping (chunked)
        logger.info("Loading STRING protein aliases (chunked)...")
        
        chunk_count = 0
        total_processed = 0
        for chunk in self._load_file_chunked('string/9606.protein.aliases.v12.0.txt.gz', chunk_size=100000):
            if chunk.empty:
                continue
            chunk_count += 1
            total_processed += len(chunk)
            logger.info(f"  Processing aliases chunk {chunk_count} ({len(chunk):,} rows, {total_processed:,} total)...")
            
            # Process STRING aliases to build mapping
            uniprot_aliases = chunk[chunk['source'] == 'UniProt_AC'].copy()
            symbol_aliases = chunk[chunk['source'].str.contains('BLAST_UniProt_GN|UniProt_GN', na=False)].copy()
            
            # Build STRING to UniProt mapping
            for _, row in uniprot_aliases.iterrows():
                string_id = row['#string_protein_id']
                uniprot_id = row['alias']
                self.id_mapping['string_to_uniprot'][string_id] = uniprot_id
                
            # Build Symbol mappings  
            for _, row in symbol_aliases.iterrows():
                string_id = row['#string_protein_id']
                symbol = row['alias']
                if string_id in self.id_mapping['string_to_uniprot']:
                    uniprot_id = self.id_mapping['string_to_uniprot'][string_id]
                    self.id_mapping['symbol_to_uniprot'][symbol] = uniprot_id
                    self.id_mapping['uniprot_to_symbol'][uniprot_id] = symbol
            
            # Force garbage collection after each chunk
            import gc
            gc.collect()
        
        logger.info(f"Built STRING→UniProt mapping: {len(self.id_mapping['string_to_uniprot'])} entries")
        logger.info(f"Built Symbol↔UniProt mapping: {len(self.id_mapping['symbol_to_uniprot'])} entries")
        
        # Step 2: Load STRING protein info for additional annotations (chunked)
        logger.info("Loading STRING protein info (chunked)...")
        
        chunk_count = 0
        for chunk in self._load_file_chunked('string/9606.protein.info.v12.0.txt.gz', chunk_size=50000):
            if chunk.empty:
                continue
            chunk_count += 1
            logger.info(f"  Processing info chunk {chunk_count} ({len(chunk):,} rows)...")
            
            for _, row in chunk.iterrows():
                string_id = row['#string_protein_id']
                if string_id in self.id_mapping['string_to_uniprot']:
                    uniprot_id = self.id_mapping['string_to_uniprot'][string_id]
                    self.id_mapping['uniprot_info'][uniprot_id] = {
                        'preferred_name': row['preferred_name'],
                        'protein_size': row['protein_size'],
                        'annotation': row['annotation']
                    }
            
            # Force garbage collection after each chunk
            import gc
            gc.collect()
        
        logger.info(f"Added protein info for {len(self.id_mapping['uniprot_info'])} UniProt entries")
        
        # Step 3: Load Reactome UniProt mappings for human only
        logger.info("Loading Reactome UniProt mappings...")
        reactome_df = self._load_file_completely('reactome/UniProt2Reactome_All_Levels.txt')
        
        # Filter for human (Homo sapiens) only
        human_reactome = reactome_df[reactome_df['Species'] == 'Homo sapiens'].copy()
        logger.info(f"Found {len(human_reactome)} human Reactome pathway mappings")
        
        # Step 4: Build comprehensive mapping statistics
        self._print_mapping_statistics()
        
        # Step 5: Create reverse lookup functions
        self._create_mapping_functions()
        
        logger.info("Phase 2 complete: ID mapping system built")
        self._monitor_memory("Phase 2 Complete")
        
    def phase3_create_mitochondrial_reference(self):
        """
        Phase 3: Create comprehensive protein reference set from MitoCarta + HPA skeletal muscle
        """
        logger.info("=== PHASE 3: COMPREHENSIVE PROTEIN REFERENCE ===")
        
        # Step 1: Load MitoCarta data
        logger.info("Loading MitoCarta 3.0 data...")
        mitocarta_df = self._load_file_completely('mitocarta/Human.MitoCarta3.0.xls')
        
        if mitocarta_df.empty:
            logger.error("Failed to load MitoCarta data")
            return
            
        logger.info(f"Loaded {len(mitocarta_df)} entries from MitoCarta 3.0")
        
        # Step 2: Load HPA skeletal muscle data
        logger.info("Loading HPA skeletal muscle expression data...")
        hpa_df = self._load_file_completely('hpa/hpa_skm.tsv')
        
        if hpa_df.empty:
            logger.error("Failed to load HPA data")
            return
            
        logger.info(f"Loaded {len(hpa_df)} entries from HPA skeletal muscle")
        
        # Step 3: Load MitoCarta pathway mappings
        logger.info("Loading MitoCarta pathway mappings...")
        pathway_mappings = self._load_mitocarta_pathways()
        
        # Initialize comprehensive protein data structure
        self.protein_reference = {}
        self.mitochondrial_uniprot_ids = set()
        self.muscle_expressed_uniprot_ids = set()
        
        # Step 4: Process MitoCarta entries
        mitocarta_mapped = 0
        mitocarta_unmapped = []
        
        for _, row in mitocarta_df.iterrows():
            symbol = row['Symbol']
            uniprot_id = self.map_to_uniprot(symbol, 'symbol')
            
            if uniprot_id:
                self.protein_reference[uniprot_id] = {
                    'gene_symbol': symbol,
                    'gene_description': row.get('Description', ''),
                    'mitocarta_member': True,
                    'mitocarta_list': row.get('MitoCarta3.0_List', ''),
                    'mitocarta_evidence': row.get('MitoCarta3.0_Evidence', ''),
                    'mitocarta_sub_localization': row.get('MitoCarta3.0_SubMitoLocalization', ''),
                    'mitocarta_pathways': row.get('MitoCarta3.0_MitoPathways', ''),
                    'mitocarta_pathways_detailed': pathway_mappings.get(symbol, []),
                    'human_gene_id': row.get('HumanGeneID', ''),
                    'synonyms': row.get('Synonyms', ''),
                    # Initialize HPA fields as None (will be filled if available)
                    'protein_evidence_level': None,
                    'muscle_TPM': 0.0,
                    'muscle_specificity_label': None,
                    'muscle_specificity_score': None,
                    'main_localization': None,
                    'additional_localization': None,
                    'hpa_interactions': None
                }
                self.mitochondrial_uniprot_ids.add(uniprot_id)
                mitocarta_mapped += 1
            else:
                mitocarta_unmapped.append(symbol)
        
        logger.info(f"Mapped {mitocarta_mapped}/{len(mitocarta_df)} MitoCarta proteins ({mitocarta_mapped/len(mitocarta_df)*100:.1f}%)")
        
        # Step 5: Process HPA skeletal muscle entries
        hpa_mapped = 0
        hpa_muscle_expressed = 0
        hpa_unmapped = []
        
        for _, row in hpa_df.iterrows():
            symbol = row['Gene']
            muscle_tpm = float(row.get('Tissue RNA - skeletal muscle [nTPM]', 0.0))
            
            # Only include if expressed in skeletal muscle (TPM > 0)
            if muscle_tpm > 0:
                uniprot_id = self.map_to_uniprot(symbol, 'symbol')
                
                if uniprot_id:
                    # If protein already exists (from MitoCarta), update with HPA data
                    if uniprot_id in self.protein_reference:
                        existing = self.protein_reference[uniprot_id]
                    else:
                        # Create new entry for muscle-expressed non-mitochondrial protein
                        existing = {
                            'gene_symbol': symbol,
                            'gene_description': '',
                            'mitocarta_member': False,
                            'mitocarta_list': '',
                            'mitocarta_evidence': '',
                            'mitocarta_sub_localization': '',
                            'mitocarta_pathways': '',
                            'mitocarta_pathways_detailed': [],
                            'human_gene_id': '',
                            'synonyms': ''
                        }
                        self.protein_reference[uniprot_id] = existing
                        
                    # Add/update HPA data
                    existing.update({
                        'gene_description': existing['gene_description'] or row.get('Gene description', ''),
                        'protein_evidence_level': row.get('Evidence', ''),
                        'muscle_TPM': muscle_tpm,
                        'muscle_specificity_label': row.get('RNA tissue specificity', ''),
                        'muscle_specificity_score': self._safe_float(row.get('RNA tissue specificity score', '')),
                        'main_localization': row.get('Subcellular main location', ''),
                        'additional_localization': row.get('Subcellular additional location', ''),
                        'hpa_interactions': row.get('Interactions', '')
                    })
                    
                    self.muscle_expressed_uniprot_ids.add(uniprot_id)
                    hpa_muscle_expressed += 1
                    hpa_mapped += 1
                else:
                    hpa_unmapped.append(symbol)
        
        logger.info(f"Mapped {hpa_mapped} HPA proteins, {hpa_muscle_expressed} muscle-expressed")
        
        # Step 6: Calculate priority scores
        for uniprot_id, data in self.protein_reference.items():
            mitocarta_score = 1.0 if data['mitocarta_member'] else 0.0
            muscle_score = np.log1p(data['muscle_TPM'])
            data['priority_score'] = mitocarta_score + muscle_score
        
        # Step 7: Create reference sets
        self.included_uniprot_ids = set(self.protein_reference.keys())
        
        # Print comprehensive statistics
        self._print_reference_statistics()
        
        logger.info("Phase 3 complete: Comprehensive protein reference created")
        self._monitor_memory("Phase 3 Complete")
        
    def phase4_integrate_networks(self):
        """
        Phase 4: Process and integrate multi-source protein-protein interaction networks
        """
        logger.info("=== PHASE 4: MULTI-SOURCE NETWORK INTEGRATION ===")
        
        # Initialize network
        self.network = nx.Graph()
        self.edge_sources = {}  # Track which sources contributed each edge
        
        # Step 1: Process STRING networks
        logger.info("Processing STRING networks...")
        string_full_edges = self._process_string_network('string/9606.protein.links.detailed.v12.0.txt.gz', 'STRING_full')
        string_physical_edges = self._process_string_network('string/9606.protein.physical.links.detailed.v12.0.txt.gz', 'STRING_physical')
        
        # Step 2: Process BioGRID interactions
        logger.info("Processing BioGRID interactions...")
        biogrid_edges = self._process_biogrid_network()
        
        # Step 3: Process Reactome interactions
        logger.info("Processing Reactome interactions...")
        reactome_edges = self._process_reactome_network()
        
        # Step 4: Process CORUM complex interactions
        logger.info("Processing CORUM complex interactions...")
        corum_edges = self._process_corum_network()
        
        # Step 5: Merge all edges with comprehensive attributes
        logger.info("Merging and attributing network edges...")
        self._merge_network_sources()
        
        # Step 6: Filter network to include only relevant proteins
        logger.info("Filtering network for mitochondrial/muscle proteins...")
        self._filter_network()
        
        self._monitor_memory("Phase 4 Complete")
        
        # Step 7: Generate network statistics
        self._print_network_statistics()
        
        logger.info("Phase 4 complete: Integrated network ready")
        
    def phase5_annotate_nodes(self):
        """
        Phase 5: Add comprehensive node annotations and enrichment
        """
        logger.info("=== PHASE 5: NODE ANNOTATION AND ENRICHMENT ===")
        
        # Load Reactome pathway data for pathway enrichment
        logger.info("Loading Reactome pathway data...")
        reactome_pathways = self._build_reactome_pathway_mappings()
        
        # Annotate each node in the network
        logger.info("Annotating network nodes...")
        nodes_annotated = 0
        
        for node_id in self.network.nodes():
            # Get comprehensive protein information
            protein_info = self.get_protein_info(node_id)
            string_info = self.get_protein_info(node_id)  # From STRING info
            
            # Build comprehensive node attributes
            node_attrs = {
                # Core identifiers
                'uniprot_id': node_id,
                'gene_symbol': protein_info.get('gene_symbol', ''),
                'gene_description': protein_info.get('gene_description', ''),
                'alternative_ids': self._get_alternative_ids(node_id),
                
                # Protein evidence and quality
                'protein_evidence_level': protein_info.get('protein_evidence_level', ''),
                'protein_size': string_info.get('protein_size', ''),
                'annotation_score': self._calculate_annotation_score(protein_info),
                
                # Mitochondrial annotations (from MitoCarta)
                'mitocarta_member': protein_info.get('mitocarta_member', False),
                'mitocarta_confidence': protein_info.get('mitocarta_evidence', ''),
                'mitocarta_sub_localization': protein_info.get('mitocarta_sub_localization', ''),
                'mitocarta_pathways': protein_info.get('mitocarta_pathways', ''),
                'mitocarta_pathways_detailed': protein_info.get('mitocarta_pathways_detailed', []),
                
                # Muscle expression (from HPA)
                'muscle_expressed': self.is_muscle_expressed(node_id),
                'muscle_TPM': protein_info.get('muscle_TPM', 0.0),
                'muscle_specificity_label': protein_info.get('muscle_specificity_label', ''),
                'muscle_specificity_score': protein_info.get('muscle_specificity_score', 0.0),
                
                # Subcellular localization
                'main_localization': protein_info.get('main_localization', ''),
                'additional_localization': protein_info.get('additional_localization', ''),
                'localization_sources': self._combine_localization_sources(protein_info),
                
                # Pathway memberships
                'reactome_pathways': reactome_pathways.get(node_id, []),
                'top_level_pathways': self._get_top_level_pathways(reactome_pathways.get(node_id, [])),
                'pathway_count': len(reactome_pathways.get(node_id, [])),
                
                # Complex memberships (from CORUM)
                'protein_complexes': self._get_corum_complexes(node_id),
                'complex_count': len(self._get_corum_complexes(node_id)),
                
                # Network topology features
                'degree': self.network.degree(node_id),
                'degree_centrality': self._calculate_degree_centrality(node_id),
                'betweenness_centrality': 0.0,  # Will be calculated later if needed
                'clustering_coefficient': 0.0,  # Will be calculated later if needed
                
                # Priority and classification
                'priority_score': protein_info.get('priority_score', 0.0),
                'protein_class': self._classify_protein(protein_info),
                'functional_category': self._get_functional_category(protein_info),
                
                # Cross-references and interactions
                'hpa_interactions': protein_info.get('hpa_interactions', ''),
                'num_sources': self._count_interaction_sources(node_id),
                'max_edge_confidence': self._get_max_edge_confidence(node_id)
            }
            
            # Add all attributes to the node
            self.network.nodes[node_id].update(node_attrs)
            nodes_annotated += 1
            
        logger.info(f"Annotated {nodes_annotated:,} nodes with comprehensive metadata")
        
        # Generate annotation statistics
        self._print_annotation_statistics()
        
        logger.info("Phase 5 complete: Node annotation finished")
        
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
        """
        Export network to GraphML format for Cytoscape
        """
        output_file = "outputs/mitonet_network.graphml"
        
        try:
            # Clean node attributes for GraphML compatibility
            cleaned_network = self.network.copy()
            
            for node_id in cleaned_network.nodes():
                node_data = cleaned_network.nodes[node_id]
                # Convert lists to comma-separated strings for GraphML
                for key, value in list(node_data.items()):
                    if isinstance(value, list):
                        node_data[key] = "|".join(str(v) for v in value)
                    elif isinstance(value, bool):
                        node_data[key] = str(value).lower()
                    elif value is None:
                        node_data[key] = ""
                        
            # Clean edge attributes
            for u, v in cleaned_network.edges():
                edge_data = cleaned_network.edges[u, v]
                for key, value in list(edge_data.items()):
                    if isinstance(value, list):
                        edge_data[key] = "|".join(str(v) for v in value)
                    elif isinstance(value, bool):
                        edge_data[key] = str(value).lower()
                    elif value is None:
                        edge_data[key] = ""
                        
            # Export to GraphML
            nx.write_graphml(cleaned_network, output_file)
            
            file_size = Path(output_file).stat().st_size / (1024*1024)  # MB
            logger.info(f"✓ GraphML export: {output_file} ({file_size:.1f} MB)")
            
        except Exception as e:
            logger.error(f"✗ GraphML export failed: {e}")
            
    def _export_json(self):
        """
        Export network to JSON format for web visualization
        """
        output_file = "outputs/mitonet_network.json"
        
        try:
            # Create JSON-compatible data structure
            network_data = {
                'metadata': {
                    'version': '1.0',
                    'created': pd.Timestamp.now().isoformat(),
                    'description': 'Mitochondrial and muscle protein interaction network',
                    'nodes': self.network.number_of_nodes(),
                    'edges': self.network.number_of_edges(),
                    'sources': ['STRING', 'BioGRID', 'Reactome', 'CORUM'],
                    'databases': ['MitoCarta3.0', 'HPA_skeletal_muscle']
                },
                'nodes': [],
                'edges': []
            }
            
            # Process nodes
            for node_id in self.network.nodes():
                node_data = self.network.nodes[node_id]
                
                # Create clean node record
                clean_node = {
                    'id': node_id,
                    'label': node_data.get('gene_symbol', node_id),
                    'type': node_data.get('protein_class', 'Unknown'),
                    'score': node_data.get('priority_score', 0.0),
                    'attributes': {}
                }
                
                # Add selected attributes
                key_attributes = [
                    'gene_symbol', 'gene_description', 'mitocarta_member', 
                    'muscle_expressed', 'muscle_TPM', 'main_localization',
                    'functional_category', 'degree', 'annotation_score'
                ]
                
                for attr in key_attributes:
                    value = node_data.get(attr)
                    if isinstance(value, list):
                        clean_node['attributes'][attr] = value[:5]  # Limit list size
                    elif value is not None:
                        clean_node['attributes'][attr] = value
                        
                network_data['nodes'].append(clean_node)
                
            # Process edges
            for u, v in self.network.edges():
                edge_data = self.network.edges[u, v]
                
                clean_edge = {
                    'source': u,
                    'target': v,
                    'confidence': edge_data.get('composite_confidence', 0.0),
                    'sources': edge_data.get('sources', []),
                    'type': edge_data.get('interaction_types', ['unknown'])[0] if edge_data.get('interaction_types') else 'unknown'
                }
                
                network_data['edges'].append(clean_edge)
                
            # Write JSON file
            with open(output_file, 'w') as f:
                json.dump(network_data, f, indent=2, default=str)
                
            file_size = Path(output_file).stat().st_size / (1024*1024)  # MB
            logger.info(f"✓ JSON export: {output_file} ({file_size:.1f} MB)")
            
        except Exception as e:
            logger.error(f"✗ JSON export failed: {e}")
            
    def _export_data_tables(self):
        """
        Export node and edge data as CSV tables
        """
        try:
            # Export nodes table
            nodes_data = []
            for node_id in self.network.nodes():
                node_data = self.network.nodes[node_id].copy()
                node_data['uniprot_id'] = node_id
                
                # Convert lists to strings
                for key, value in node_data.items():
                    if isinstance(value, list):
                        node_data[key] = "|".join(str(v) for v in value)
                        
                nodes_data.append(node_data)
                
            nodes_df = pd.DataFrame(nodes_data)
            nodes_file = "outputs/mitonet_nodes.csv"
            nodes_df.to_csv(nodes_file, index=False)
            logger.info(f"✓ Nodes table: {nodes_file} ({len(nodes_df):,} rows)")
            
            # Export edges table
            edges_data = []
            for u, v in self.network.edges():
                edge_data = self.network.edges[u, v].copy()
                edge_data['source'] = u
                edge_data['target'] = v
                
                # Convert lists to strings
                for key, value in edge_data.items():
                    if isinstance(value, list):
                        edge_data[key] = "|".join(str(v) for v in value)
                        
                edges_data.append(edge_data)
                
            edges_df = pd.DataFrame(edges_data)
            edges_file = "outputs/mitonet_edges.csv"
            edges_df.to_csv(edges_file, index=False)
            logger.info(f"✓ Edges table: {edges_file} ({len(edges_df):,} rows)")
            
        except Exception as e:
            logger.error(f"✗ Table export failed: {e}")
            
    def _generate_summary_report(self):
        """
        Generate comprehensive summary report
        """
        report_file = "outputs/mitonet_summary_report.md"
        
        try:
            with open(report_file, 'w') as f:
                # Header
                f.write("# Mitochondrial Network Integration Pipeline - Summary Report\\n\\n")
                f.write(f"**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\\n\\n")
                
                # Network Overview
                f.write("## Network Overview\\n\\n")
                f.write(f"- **Total nodes:** {self.network.number_of_nodes():,}\\n")
                f.write(f"- **Total edges:** {self.network.number_of_edges():,}\\n")
                
                # Calculate basic statistics
                degrees = dict(self.network.degree())
                avg_degree = sum(degrees.values()) / len(degrees) if degrees else 0
                f.write(f"- **Average degree:** {avg_degree:.2f}\\n")
                f.write(f"- **Network density:** {nx.density(self.network):.6f}\\n")
                
                # Connected components
                components = list(nx.connected_components(self.network))
                largest_comp = max(len(comp) for comp in components) if components else 0
                f.write(f"- **Connected components:** {len(components)}\\n")
                f.write(f"- **Largest component:** {largest_comp:,} nodes ({largest_comp/self.network.number_of_nodes()*100:.1f}%)\\n\\n")
                
                # Protein Composition
                f.write("## Protein Composition\\n\\n")
                mito_count = sum(1 for n in self.network.nodes() if self.is_mitochondrial(n))
                muscle_count = sum(1 for n in self.network.nodes() if self.is_muscle_expressed(n))
                overlap_count = sum(1 for n in self.network.nodes() if self.is_mitochondrial(n) and self.is_muscle_expressed(n))
                
                f.write(f"- **Mitochondrial proteins:** {mito_count:,} ({mito_count/self.network.number_of_nodes()*100:.1f}%)\\n")
                f.write(f"- **Muscle-expressed proteins:** {muscle_count:,} ({muscle_count/self.network.number_of_nodes()*100:.1f}%)\\n")
                f.write(f"- **Mitochondrial + Muscle overlap:** {overlap_count:,} ({overlap_count/self.network.number_of_nodes()*100:.1f}%)\\n\\n")
                
                # Data Sources
                f.write("## Data Sources\\n\\n")
                source_counts = {}
                for _, _, data in self.network.edges(data=True):
                    for source in data.get('sources', []):
                        source_counts[source] = source_counts.get(source, 0) + 1
                        
                f.write("### Edge Sources\\n")
                for source, count in sorted(source_counts.items(), key=lambda x: x[1], reverse=True):
                    f.write(f"- **{source}:** {count:,} edges\\n")
                    
                # Multi-source validation
                multi_source_edges = sum(1 for _, _, data in self.network.edges(data=True) 
                                       if data.get('num_sources', 0) > 1)
                f.write(f"\\n### Validation\\n")
                f.write(f"- **Multi-source edges:** {multi_source_edges:,} ({multi_source_edges/self.network.number_of_edges()*100:.1f}%)\\n\\n")
                
                # Confidence Distribution
                f.write("## Edge Confidence Distribution\\n\\n")
                confidences = [data.get('composite_confidence', 0.0) 
                             for _, _, data in self.network.edges(data=True)]
                
                high_conf = sum(1 for c in confidences if c >= 0.8)
                med_conf = sum(1 for c in confidences if 0.5 <= c < 0.8)
                low_conf = sum(1 for c in confidences if c < 0.5)
                
                f.write(f"- **High confidence (≥0.8):** {high_conf:,} ({high_conf/len(confidences)*100:.1f}%)\\n")
                f.write(f"- **Medium confidence (0.5-0.8):** {med_conf:,} ({med_conf/len(confidences)*100:.1f}%)\\n")
                f.write(f"- **Low confidence (<0.5):** {low_conf:,} ({low_conf/len(confidences)*100:.1f}%)\\n\\n")
                
                # Top proteins by degree
                f.write("## Top Hub Proteins\\n\\n")
                top_hubs = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:10]
                
                f.write("| UniProt ID | Gene Symbol | Degree | Type |\\n")
                f.write("|------------|-------------|--------|------|\\n")
                
                for uniprot_id, degree in top_hubs:
                    symbol = self.network.nodes[uniprot_id].get('gene_symbol', uniprot_id)
                    protein_type = self.network.nodes[uniprot_id].get('protein_class', 'Unknown')
                    f.write(f"| {uniprot_id} | {symbol} | {degree} | {protein_type} |\\n")
                    
                f.write("\\n")
                
                # File Outputs
                f.write("## Generated Files\\n\\n")
                f.write("- `outputs/mitonet_network.graphml` - Network file for Cytoscape visualization\\n")
                f.write("- `outputs/mitonet_network.json` - Network file for web-based visualization\\n")
                f.write("- `outputs/mitonet_nodes.csv` - Node attributes table\\n")
                f.write("- `outputs/mitonet_edges.csv` - Edge attributes table\\n")
                f.write("- `outputs/mitonet_summary_report.md` - This summary report\\n\\n")
                
                # Usage Instructions
                f.write("## Usage Instructions\\n\\n")
                f.write("### Cytoscape Visualization\\n")
                f.write("1. Open Cytoscape\\n")
                f.write("2. Import → Network from File → Select `outputs/mitonet_network.graphml`\\n")
                f.write("3. Apply layout (e.g., Prefuse Force Directed)\\n")
                f.write("4. Style nodes by `protein_class` and `priority_score`\\n")
                f.write("5. Style edges by `composite_confidence`\\n\\n")
                
                f.write("### Web Visualization\\n")
                f.write("- Use `outputs/mitonet_network.json` with D3.js, Cytoscape.js, or similar libraries\\n")
                f.write("- Node attributes include classification, scores, and annotations\\n")
                f.write("- Edge attributes include confidence scores and source information\\n\\n")
                
                # Analysis Recommendations
                f.write("## Analysis Recommendations\\n\\n")
                f.write("1. **Network Topology Analysis**\\n")
                f.write("   - Identify highly connected hub proteins\\n")
                f.write("   - Analyze community structure and clustering\\n")
                f.write("   - Examine shortest paths between protein classes\\n\\n")
                
                f.write("2. **Functional Analysis**\\n")
                f.write("   - Pathway enrichment of network communities\\n")
                f.write("   - Gene Ontology analysis of protein clusters\\n")
                f.write("   - Disease association analysis\\n\\n")
                
                f.write("3. **Integration Analysis**\\n")
                f.write("   - Compare mitochondrial vs muscle expression patterns\\n")
                f.write("   - Analyze tissue-specific interaction modules\\n")
                f.write("   - Identify candidate therapeutic targets\\n\\n")
                
            logger.info(f"✓ Summary report: {report_file}")
            
        except Exception as e:
            logger.error(f"✗ Report generation failed: {e}")
        
    def _build_reactome_pathway_mappings(self) -> Dict[str, List[str]]:
        """
        Build UniProt → Reactome pathway mappings
        """
        reactome_df = self._load_file_completely('reactome/UniProt2Reactome_All_Levels.txt')
        
        # Filter for human pathways only
        human_reactome = reactome_df[reactome_df['Species'] == 'Homo sapiens'].copy()
        
        pathway_mappings = {}
        for _, row in human_reactome.iterrows():
            uniprot_id = row['UniProt']
            pathway_name = row['Event_Name']
            
            if uniprot_id not in pathway_mappings:
                pathway_mappings[uniprot_id] = []
            pathway_mappings[uniprot_id].append(pathway_name)
            
        logger.info(f"Built Reactome mappings for {len(pathway_mappings):,} proteins")
        return pathway_mappings
        
    def _get_alternative_ids(self, uniprot_id: str) -> List[str]:
        """
        Get alternative identifiers for a protein
        """
        alt_ids = []
        
        # Add gene symbol
        symbol = self.get_symbol(uniprot_id)
        if symbol:
            alt_ids.append(f"SYMBOL:{symbol}")
            
        # Add STRING ID if available
        for string_id, mapped_uniprot in self.id_mapping['string_to_uniprot'].items():
            if mapped_uniprot == uniprot_id:
                alt_ids.append(f"STRING:{string_id}")
                break
                
        return alt_ids
        
    def _calculate_annotation_score(self, protein_info: Dict) -> float:
        """
        Calculate comprehensive annotation quality score
        """
        score = 0.0
        
        # MitoCarta evidence
        if protein_info.get('mitocarta_member', False):
            score += 0.3
            
        # HPA evidence level
        evidence = protein_info.get('protein_evidence_level', '')
        if 'protein level' in evidence:
            score += 0.3
        elif 'transcript level' in evidence:
            score += 0.2
            
        # Localization data
        if protein_info.get('main_localization'):
            score += 0.2
            
        # Pathway data
        if protein_info.get('mitocarta_pathways_detailed'):
            score += 0.2
            
        return min(score, 1.0)
        
    def _combine_localization_sources(self, protein_info: Dict) -> List[str]:
        """
        Combine localization information from multiple sources
        """
        sources = []
        
        # MitoCarta localization
        mito_loc = protein_info.get('mitocarta_sub_localization', '')
        if mito_loc:
            sources.append(f"MitoCarta:{mito_loc}")
            
        # HPA localization
        hpa_main = protein_info.get('main_localization', '')
        if hpa_main:
            sources.append(f"HPA_main:{hpa_main}")
            
        hpa_additional = protein_info.get('additional_localization', '')
        if hpa_additional:
            sources.append(f"HPA_additional:{hpa_additional}")
            
        return sources
        
    def _get_top_level_pathways(self, pathways: List[str]) -> List[str]:
        """
        Extract top-level pathway categories
        """
        top_level = set()
        
        for pathway in pathways:
            # Extract high-level categories from pathway names
            if 'metabolism' in pathway.lower():
                top_level.add('Metabolism')
            elif 'signal' in pathway.lower():
                top_level.add('Signaling')
            elif 'transport' in pathway.lower():
                top_level.add('Transport')
            elif 'immune' in pathway.lower():
                top_level.add('Immune Response')
            elif 'cell cycle' in pathway.lower():
                top_level.add('Cell Cycle')
            elif 'apoptosis' in pathway.lower():
                top_level.add('Cell Death')
            else:
                top_level.add('Other')
                
        return list(top_level)
        
    def _get_corum_complexes(self, uniprot_id: str) -> List[str]:
        """
        Get CORUM complex memberships for a protein
        """
        complexes = []
        
        # Load CORUM mapping if not already done
        mapping_df = self._load_file_completely('corum/corum_uniprotCorumMapping.txt')
        complexes_df = self._load_file_completely('corum/corum_humanComplexes.txt')
        
        # Find complexes for this protein
        protein_complexes = mapping_df[mapping_df['UniProtKB_accession_number'] == uniprot_id]
        
        for _, row in protein_complexes.iterrows():
            complex_id = row['corum_id']
            complex_info = complexes_df[complexes_df['complex_id'] == complex_id]
            
            if not complex_info.empty:
                complex_name = complex_info.iloc[0]['complex_name']
                complexes.append(complex_name)
                
        return complexes
        
    def _calculate_degree_centrality(self, node_id: str) -> float:
        """
        Calculate degree centrality for a node
        """
        if self.network.number_of_nodes() <= 1:
            return 0.0
            
        degree = self.network.degree(node_id)
        max_possible = self.network.number_of_nodes() - 1
        return degree / max_possible if max_possible > 0 else 0.0
        
    def _classify_protein(self, protein_info: Dict) -> str:
        """
        Classify protein based on its characteristics
        """
        uniprot_id = protein_info.get('gene_symbol', '')  # Use gene symbol as fallback for lookup
        
        # Find the actual UniProt ID for this protein
        for uid, pinfo in self.protein_reference.items():
            if pinfo.get('gene_symbol') == uniprot_id or uid == uniprot_id:
                if pinfo.get('mitocarta_member', False) and self.is_muscle_expressed(uid):
                    return 'Mitochondrial_Muscle'
                elif pinfo.get('mitocarta_member', False):
                    return 'Mitochondrial_Only'
                elif self.is_muscle_expressed(uid):
                    return 'Muscle_Only'
                else:
                    return 'Associated'
        return 'Unknown'
            
    def _get_functional_category(self, protein_info: Dict) -> str:
        """
        Determine functional category based on pathways and localization
        """
        pathways = protein_info.get('mitocarta_pathways', '').lower()
        
        if 'oxphos' in pathways:
            return 'OXPHOS'
        elif 'metabolism' in pathways:
            return 'Metabolism'
        elif 'transport' in pathways:
            return 'Transport'
        elif 'ribosome' in pathways:
            return 'Translation'
        elif 'import' in pathways:
            return 'Protein Import'
        else:
            return 'Other'
            
    def _count_interaction_sources(self, node_id: str) -> int:
        """
        Count how many data sources contribute interactions for this node
        """
        sources = set()
        
        for neighbor in self.network.neighbors(node_id):
            edge_data = self.network.edges[node_id, neighbor]
            edge_sources = edge_data.get('sources', [])
            sources.update(edge_sources)
            
        return len(sources)
        
    def _get_max_edge_confidence(self, node_id: str) -> float:
        """
        Get maximum edge confidence for this node's interactions
        """
        max_confidence = 0.0
        
        for neighbor in self.network.neighbors(node_id):
            edge_data = self.network.edges[node_id, neighbor]
            confidence = edge_data.get('composite_confidence', 0.0)
            max_confidence = max(max_confidence, confidence)
            
        return max_confidence
        
    def _print_annotation_statistics(self):
        """
        Print comprehensive annotation statistics
        """
        logger.info("\n" + "="*70)
        logger.info("NODE ANNOTATION STATISTICS")
        logger.info("="*70)
        
        # Basic counts
        total_nodes = self.network.number_of_nodes()
        
        # Protein classification distribution
        class_counts = {}
        for node_id in self.network.nodes():
            node_data = self.network.nodes[node_id]
            protein_class = node_data.get('protein_class', 'Unknown')
            class_counts[protein_class] = class_counts.get(protein_class, 0) + 1
            
        logger.info("Protein classification:")
        for pclass, count in sorted(class_counts.items()):
            percentage = count / total_nodes * 100
            logger.info(f"  {pclass}: {count:,} ({percentage:.1f}%)")
            
        # Functional category distribution
        func_counts = {}
        for node_id in self.network.nodes():
            node_data = self.network.nodes[node_id]
            func_cat = node_data.get('functional_category', 'Unknown')
            func_counts[func_cat] = func_counts.get(func_cat, 0) + 1
            
        logger.info("\nFunctional categories:")
        for func, count in sorted(func_counts.items(), key=lambda x: x[1], reverse=True):
            percentage = count / total_nodes * 100
            logger.info(f"  {func}: {count:,} ({percentage:.1f}%)")
            
        # Annotation quality
        annotation_scores = [self.network.nodes[n].get('annotation_score', 0.0) for n in self.network.nodes()]
        if annotation_scores:
            logger.info(f"\nAnnotation quality:")
            logger.info(f"  Mean score: {np.mean(annotation_scores):.3f}")
            logger.info(f"  High quality (>0.8): {sum(1 for s in annotation_scores if s > 0.8):,}")
            logger.info(f"  Medium quality (0.5-0.8): {sum(1 for s in annotation_scores if 0.5 <= s <= 0.8):,}")
            logger.info(f"  Low quality (<0.5): {sum(1 for s in annotation_scores if s < 0.5):,}")
            
        # Pathway coverage
        pathway_counts = [len(self.network.nodes[n].get('reactome_pathways', [])) for n in self.network.nodes()]
        if pathway_counts:
            logger.info(f"\nPathway annotation coverage:")
            logger.info(f"  Proteins with pathways: {sum(1 for c in pathway_counts if c > 0):,}")
            if any(c > 0 for c in pathway_counts):
                logger.info(f"  Mean pathways per protein: {np.mean([c for c in pathway_counts if c > 0]):.1f}")
            
        logger.info("="*70)
        
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
        
        # Giant component analysis
        if component_sizes[0] > 1:
            giant_component = max(components, key=len)
            giant_graph = self.network.subgraph(giant_component)
            
            # Average path length (on a sample for efficiency)
            if len(giant_component) <= 1000:
                try:
                    avg_path_length = nx.average_shortest_path_length(giant_graph)
                    logger.info(f"Average path length (giant component): {avg_path_length:.2f}")
                except:
                    logger.info("Average path length: Could not compute")
            else:
                logger.info("Average path length: Skipped (large network)")
                
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
            'OXPHOS': ['NDUFA1', 'NDUFB1', 'COX1', 'ATP5F1A', 'CYCS'],
            'TCA_cycle': ['CS', 'IDH2', 'OGDH', 'SUCLA2', 'MDH2'],
            'Import': ['TOMM40', 'TIMM23', 'TIMM44', 'HSPA9', 'DNAJC19'],
            'Dynamics': ['DNM1L', 'MFN1', 'MFN2', 'OPA1'],
            'Metabolism': ['GLUD1', 'PCCA', 'ACADM', 'CPT1A']
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
        
        # Network size (10 points)
        if num_nodes >= 1000:
            quality_score += 10
        elif num_nodes >= 500:
            quality_score += 7
        elif num_nodes >= 100:
            quality_score += 5
        else:
            quality_score += 2
            
        # Multi-source validation (30 points)
        if multi_source_percentage >= 80:
            quality_score += 30
        elif multi_source_percentage >= 60:
            quality_score += 25
        elif multi_source_percentage >= 40:
            quality_score += 20
        else:
            quality_score += 10
            
        # High confidence edges (20 points)
        high_conf_percentage = high_confidence_edges / num_edges * 100 if num_edges > 0 else 0
        if high_conf_percentage >= 50:
            quality_score += 20
        elif high_conf_percentage >= 30:
            quality_score += 15
        elif high_conf_percentage >= 10:
            quality_score += 10
        else:
            quality_score += 5
            
        # Connectivity (20 points)
        components = list(nx.connected_components(self.network))
        largest_comp_size = max(len(comp) for comp in components) if components else 0
        connectivity_percentage = largest_comp_size / num_nodes * 100 if num_nodes > 0 else 0
        
        if connectivity_percentage >= 90:
            quality_score += 20
        elif connectivity_percentage >= 70:
            quality_score += 15
        elif connectivity_percentage >= 50:
            quality_score += 10
        else:
            quality_score += 5
            
        # Annotation completeness (20 points)
        annotated_nodes = sum(1 for n in self.network.nodes() 
                            if self.network.nodes[n].get('gene_symbol'))
        annotation_percentage = annotated_nodes / num_nodes * 100 if num_nodes > 0 else 0
        
        if annotation_percentage >= 95:
            quality_score += 20
        elif annotation_percentage >= 85:
            quality_score += 15
        elif annotation_percentage >= 70:
            quality_score += 10
        else:
            quality_score += 5
            
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
        
    def _process_string_network(self, filename: str, source_label: str) -> int:
        """
        Process STRING interaction data using chunked loading for memory efficiency
        """
        logger.info(f"  Loading {source_label} (chunked)...")
        
        edges_added = 0
        edges_filtered = 0
        chunk_count = 0
        
        # Process file in chunks to avoid memory overflow
        for chunk in self._load_file_chunked(filename, chunk_size=50000):
            if chunk.empty:
                continue
                
            chunk_count += 1
            logger.info(f"    Processing chunk {chunk_count} ({len(chunk):,} rows)...")
            
            for _, row in chunk.iterrows():
                # Map STRING IDs to UniProt
                uniprot1 = self.map_to_uniprot(row['protein1'], 'string')
                uniprot2 = self.map_to_uniprot(row['protein2'], 'string')
                
                if not uniprot1 or not uniprot2:
                    continue
                    
                # Apply inclusion filter: at least one protein must be in our reference set
                if not (self.is_included(uniprot1) or self.is_included(uniprot2)):
                    edges_filtered += 1
                    continue
                    
                # Create edge key
                edge_key = tuple(sorted([uniprot1, uniprot2]))
                
                # Extract edge attributes
                edge_data = {
                    'source': source_label,
                    'confidence_score': row.get('combined_score', 0),
                    'experimental_score': row.get('experimental', 0),
                    'database_score': row.get('database', 0),
                    'textmining_score': row.get('textmining', 0),
                    'evidence_type': self._classify_string_evidence(row),
                    'interaction_type': 'physical' if 'physical' in source_label else 'functional',
                    'source_specific_id': f"{row['protein1']}___{row['protein2']}"
                }
                
                # Add additional scores if available
                if 'neighborhood' in row:
                    edge_data['neighborhood_score'] = row.get('neighborhood', 0)
                if 'fusion' in row:
                    edge_data['fusion_score'] = row.get('fusion', 0)
                if 'cooccurence' in row:
                    edge_data['cooccurrence_score'] = row.get('cooccurence', 0)
                if 'coexpression' in row:
                    edge_data['coexpression_score'] = row.get('coexpression', 0)
                    
                # Store edge data
                if edge_key not in self.edge_sources:
                    self.edge_sources[edge_key] = []
                self.edge_sources[edge_key].append(edge_data)
                
                edges_added += 1
            
            # Force garbage collection after each chunk
            import gc
            gc.collect()
            
        logger.info(f"  {source_label}: {edges_added:,} edges added, {edges_filtered:,} filtered from {chunk_count} chunks")
        return edges_added
        
    def _classify_string_evidence(self, row) -> str:
        """
        Classify STRING evidence type based on score distribution
        """
        experimental = row.get('experimental', 0)
        database = row.get('database', 0)
        textmining = row.get('textmining', 0)
        
        if experimental > max(database, textmining):
            return 'experimental'
        elif database > textmining:
            return 'database'
        else:
            return 'text_mining'
            
    def _process_biogrid_network(self) -> int:
        """
        Process BioGRID interaction data using chunked loading for memory efficiency
        """
        logger.info("  Loading BioGRID (chunked)...")
        
        edges_added = 0
        edges_filtered = 0
        chunk_count = 0
        
        # Process file in chunks
        for chunk in self._load_file_chunked('biogrid/BIOGRID-ALL-4.4.246.tab3.txt', chunk_size=25000):
            if chunk.empty:
                continue
                
            chunk_count += 1
            logger.info(f"    Processing BioGRID chunk {chunk_count} ({len(chunk):,} rows)...")
            
            # Filter for human interactions if organism columns exist
            if 'Organism Interactor A' in chunk.columns:
                human_chunk = chunk[chunk['Organism Interactor A'] == 'Homo sapiens'].copy()
                human_chunk = human_chunk[human_chunk['Organism Interactor B'] == 'Homo sapiens'].copy()
            else:
                # If no organism columns, assume all are human
                human_chunk = chunk.copy()
                if chunk_count == 1:
                    logger.info("    No organism columns found, assuming all interactions are human")
            
            for _, row in human_chunk.iterrows():
                # Map gene symbols to UniProt
                symbol1 = row['Official Symbol Interactor A']
                symbol2 = row['Official Symbol Interactor B']
                
                uniprot1 = self.map_to_uniprot(symbol1, 'symbol')
                uniprot2 = self.map_to_uniprot(symbol2, 'symbol')
                
                if not uniprot1 or not uniprot2:
                    continue
                    
                # Apply inclusion filter
                if not (self.is_included(uniprot1) or self.is_included(uniprot2)):
                    edges_filtered += 1
                    continue
                    
                # Create edge key
                edge_key = tuple(sorted([uniprot1, uniprot2]))
                
                # Extract edge attributes
                edge_data = {
                    'source': 'BioGRID',
                    'confidence_score': 1.0,  # BioGRID doesn't provide scores, use 1.0 for curated
                    'experimental_score': 1.0,
                    'evidence_type': 'experimental',
                    'interaction_type': self._classify_biogrid_interaction_type(row),
                    'experimental_system': row.get('Experimental System', ''),
                    'publication': row.get('Publication Source', ''),
                    'source_specific_id': str(row.get('#BioGRID Interaction ID', ''))
                }
                
                # Store edge data
                if edge_key not in self.edge_sources:
                    self.edge_sources[edge_key] = []
                self.edge_sources[edge_key].append(edge_data)
                
                edges_added += 1
            
            # Force garbage collection after each chunk
            import gc
            gc.collect()
            
        logger.info(f"  BioGRID: {edges_added:,} edges added, {edges_filtered:,} filtered from {chunk_count} chunks")
        return edges_added
        
    def _classify_biogrid_interaction_type(self, row) -> str:
        """
        Classify BioGRID interaction type
        """
        exp_system = str(row.get('Experimental System', '')).lower()
        
        if any(term in exp_system for term in ['physical', 'binding', 'immunoprecipitation', 'pull down']):
            return 'physical'
        elif any(term in exp_system for term in ['genetic', 'epistatic', 'synthetic']):
            return 'genetic'
        else:
            return 'physical'  # Default to physical for BioGRID
            
    def _process_reactome_network(self) -> int:
        """
        Process Reactome interaction data
        """
        logger.info("  Loading Reactome...")
        df = self._load_file_completely('reactome/reactome.homo_sapiens.interactions.tab-delimited.txt')
        
        if df.empty:
            logger.warning("  No Reactome data loaded")
            return 0
            
        edges_added = 0
        edges_filtered = 0
        
        for _, row in df.iterrows():
            # Extract UniProt IDs
            uniprot1_raw = row['# Interactor 1 uniprot id']
            uniprot2_raw = row['Interactor 2 uniprot id']
            
            # Clean UniProt IDs (remove uniprotkb: prefix)
            uniprot1 = uniprot1_raw.replace('uniprotkb:', '') if pd.notna(uniprot1_raw) else None
            uniprot2 = uniprot2_raw.replace('uniprotkb:', '') if pd.notna(uniprot2_raw) else None
            
            if not uniprot1 or not uniprot2:
                continue
                
            # Apply inclusion filter
            if not (self.is_included(uniprot1) or self.is_included(uniprot2)):
                edges_filtered += 1
                continue
                
            # Create edge key
            edge_key = tuple(sorted([uniprot1, uniprot2]))
            
            # Extract edge attributes
            edge_data = {
                'source': 'Reactome',
                'confidence_score': 0.8,  # Reactome is curated, assign high confidence
                'experimental_score': 0.8,
                'evidence_type': 'database',
                'interaction_type': 'pathway_derived',
                'reactome_interaction_type': row.get('Interaction type', ''),
                'reactome_context': row.get('Interaction context', ''),
                'pubmed_refs': row.get('Pubmed references', ''),
                'source_specific_id': f"{uniprot1}___{uniprot2}"
            }
            
            # Check if proteins are in same pathway
            edge_data['same_pathway_flag'] = self._check_same_pathway(uniprot1, uniprot2)
            
            # Store edge data
            if edge_key not in self.edge_sources:
                self.edge_sources[edge_key] = []
            self.edge_sources[edge_key].append(edge_data)
            
            edges_added += 1
            
        logger.info(f"  Reactome: {edges_added:,} edges added, {edges_filtered:,} filtered")
        return edges_added
        
    def _check_same_pathway(self, uniprot1: str, uniprot2: str) -> bool:
        """
        Check if two proteins are in the same Reactome pathway
        """
        # This would require loading the full Reactome pathway data
        # For now, return False as placeholder
        return False
        
    def _process_corum_network(self) -> int:
        """
        Process CORUM complex data to create co-complex interactions
        """
        logger.info("  Loading CORUM...")
        
        # Load CORUM complexes and mapping
        complexes_df = self._load_file_completely('corum/corum_humanComplexes.txt')
        mapping_df = self._load_file_completely('corum/corum_uniprotCorumMapping.txt')
        
        if complexes_df.empty or mapping_df.empty:
            logger.warning("  No CORUM data loaded")
            return 0
            
        # Filter for human complexes
        human_complexes = complexes_df[complexes_df['organism'] == 'Human'].copy()
        
        # Build complex membership mapping
        complex_members = {}
        for _, row in mapping_df.iterrows():
            uniprot_id = row['UniProtKB_accession_number']
            complex_id = row['corum_id']
            
            if complex_id not in complex_members:
                complex_members[complex_id] = []
            complex_members[complex_id].append(uniprot_id)
            
        edges_added = 0
        edges_filtered = 0
        
        # Generate pairwise interactions within each complex
        for _, complex_row in human_complexes.iterrows():
            complex_id = complex_row['complex_id']
            complex_name = complex_row['complex_name']
            
            if complex_id not in complex_members:
                continue
                
            members = complex_members[complex_id]
            
            # Filter members to include only those in our reference set
            included_members = [m for m in members if self.is_included(m)]
            
            if len(included_members) < 2:
                continue
                
            # Create all pairwise combinations
            for i, uniprot1 in enumerate(included_members):
                for uniprot2 in included_members[i+1:]:
                    edge_key = tuple(sorted([uniprot1, uniprot2]))
                    
                    # Extract edge attributes
                    edge_data = {
                        'source': 'CORUM',
                        'confidence_score': 0.9,  # High confidence for co-complex
                        'experimental_score': 0.9,
                        'evidence_type': 'database',
                        'interaction_type': 'co_complex',
                        'complex_name': complex_name,
                        'complex_id': complex_id,
                        'source_specific_id': f"CORUM_{complex_id}"
                    }
                    
                    # Store edge data
                    if edge_key not in self.edge_sources:
                        self.edge_sources[edge_key] = []
                    self.edge_sources[edge_key].append(edge_data)
                    
                    edges_added += 1
                    
        logger.info(f"  CORUM: {edges_added:,} edges added from complexes")
        return edges_added
        
    def _merge_network_sources(self):
        """
        Merge edges from all sources into the network with comprehensive attributes
        """
        logger.info("  Merging edges from all sources...")
        
        for edge_key, source_list in self.edge_sources.items():
            uniprot1, uniprot2 = edge_key
            
            # Combine attributes from all sources
            merged_attrs = {
                'sources': [s['source'] for s in source_list],
                'num_sources': len(source_list),
                'max_confidence': max(s['confidence_score'] for s in source_list),
                'mean_confidence': sum(s['confidence_score'] for s in source_list) / len(source_list),
                'evidence_types': list(set(s['evidence_type'] for s in source_list)),
                'interaction_types': list(set(s['interaction_type'] for s in source_list))
            }
            
            # Add source-specific attributes
            for source_data in source_list:
                source_name = source_data['source']
                for key, value in source_data.items():
                    if key != 'source':
                        merged_attrs[f"{source_name}_{key}"] = value
                        
            # Calculate composite confidence score
            merged_attrs['composite_confidence'] = self._calculate_composite_confidence(source_list)
            
            # Add edge to network
            self.network.add_edge(uniprot1, uniprot2, **merged_attrs)
            
    def _calculate_composite_confidence(self, source_list: List[Dict]) -> float:
        """
        Calculate composite confidence score from multiple sources
        """
        # Weight sources differently
        source_weights = {
            'STRING_physical': 0.9,
            'STRING_full': 0.7,
            'BioGRID': 1.0,
            'Reactome': 0.8,
            'CORUM': 0.9
        }
        
        weighted_sum = 0
        total_weight = 0
        
        for source_data in source_list:
            source = source_data['source']
            confidence = source_data['confidence_score']
            weight = source_weights.get(source, 0.5)
            
            weighted_sum += confidence * weight
            total_weight += weight
            
        return weighted_sum / total_weight if total_weight > 0 else 0.5
        
    def _filter_network(self):
        """
        Filter network to keep only edges involving included proteins
        """
        # Remove nodes not in our reference set
        nodes_to_remove = [n for n in self.network.nodes() if not self.is_included(n)]
        self.network.remove_nodes_from(nodes_to_remove)
        
        logger.info(f"  Removed {len(nodes_to_remove):,} nodes not in reference set")
        
    def _print_network_statistics(self):
        """
        Print comprehensive network statistics
        """
        logger.info("\n" + "="*70)
        logger.info("INTEGRATED NETWORK STATISTICS")
        logger.info("="*70)
        
        num_nodes = self.network.number_of_nodes()
        num_edges = self.network.number_of_edges()
        
        logger.info(f"Total nodes: {num_nodes:,}")
        logger.info(f"Total edges: {num_edges:,}")
        
        if num_nodes > 0:
            density = nx.density(self.network)
            logger.info(f"Network density: {density:.4f}")
            
            # Calculate degree statistics
            degrees = dict(self.network.degree())
            if degrees:
                avg_degree = sum(degrees.values()) / len(degrees)
                max_degree = max(degrees.values())
                logger.info(f"Average degree: {avg_degree:.1f}")
                logger.info(f"Maximum degree: {max_degree}")
                
        # Source statistics
        logger.info("\nSource contributions:")
        source_counts = {}
        for _, _, data in self.network.edges(data=True):
            for source in data.get('sources', []):
                source_counts[source] = source_counts.get(source, 0) + 1
                
        for source, count in sorted(source_counts.items(), key=lambda x: x[1], reverse=True):
            logger.info(f"  {source}: {count:,} edges")
            
        # Multi-source edges
        multi_source_edges = sum(1 for _, _, data in self.network.edges(data=True) 
                                if data.get('num_sources', 0) > 1)
        if num_edges > 0:
            logger.info(f"\nMulti-source edges: {multi_source_edges:,} ({multi_source_edges/num_edges*100:.1f}%)")
        else:
            logger.info(f"\nMulti-source edges: {multi_source_edges:,} (0.0%)")
        
        # Protein type distribution
        mito_nodes = sum(1 for n in self.network.nodes() if self.is_mitochondrial(n))
        muscle_nodes = sum(1 for n in self.network.nodes() if self.is_muscle_expressed(n))
        logger.info(f"\nNode types:")
        logger.info(f"  Mitochondrial: {mito_nodes:,}")
        logger.info(f"  Muscle-expressed: {muscle_nodes:,}")
        
        logger.info("="*70)
        
    def _load_mitocarta_pathways(self) -> Dict[str, List[str]]:
        """
        Load MitoCarta pathway mappings from GMX file
        """
        try:
            gmx_path = self.data_dir / 'mitocarta/Human.MitoPathways3.0.gmx'
            
            # Read the GMX file
            with open(gmx_path, 'r') as f:
                lines = f.readlines()
            
            if len(lines) < 3:
                return {}
                
            # Parse pathway hierarchy from first two lines
            pathway_names = lines[0].strip().split('\t')
            pathway_full_names = lines[1].strip().split('\t')
            
            # Parse gene memberships from subsequent lines
            gene_pathways = {}
            
            for line_idx in range(2, len(lines)):
                genes = lines[line_idx].strip().split('\t')
                
                for col_idx, gene in enumerate(genes):
                    if gene and col_idx < len(pathway_names):
                        pathway_name = pathway_names[col_idx]
                        pathway_full = pathway_full_names[col_idx]
                        
                        if gene not in gene_pathways:
                            gene_pathways[gene] = []
                        gene_pathways[gene].append(pathway_full)
            
            logger.info(f"Loaded pathway mappings for {len(gene_pathways)} genes")
            return gene_pathways
            
        except Exception as e:
            logger.warning(f"Could not load MitoCarta pathways: {e}")
            return {}
            
    def _safe_float(self, value) -> Optional[float]:
        """
        Safely convert value to float, return None if not possible
        """
        try:
            if pd.isna(value) or value == '' or value == 'nan':
                return None
            return float(value)
        except:
            return None
            
    def _print_reference_statistics(self):
        """
        Print comprehensive reference set statistics
        """
        logger.info("\n" + "="*70)
        logger.info("COMPREHENSIVE PROTEIN REFERENCE STATISTICS")
        logger.info("="*70)
        
        total_proteins = len(self.protein_reference)
        mitochondrial_count = len(self.mitochondrial_uniprot_ids)
        muscle_expressed_count = len(self.muscle_expressed_uniprot_ids)
        
        # Calculate overlaps
        mito_and_muscle = len(self.mitochondrial_uniprot_ids & self.muscle_expressed_uniprot_ids)
        mito_only = mitochondrial_count - mito_and_muscle
        muscle_only = muscle_expressed_count - mito_and_muscle
        
        logger.info(f"Total proteins in reference: {total_proteins:,}")
        logger.info(f"Mitochondrial proteins: {mitochondrial_count:,}")
        logger.info(f"Muscle-expressed proteins: {muscle_expressed_count:,}")
        logger.info(f"Mitochondrial + Muscle overlap: {mito_and_muscle:,}")
        logger.info(f"Mitochondrial only: {mito_only:,}")
        logger.info(f"Muscle only: {muscle_only:,}")
        
        # Analyze muscle expression levels
        muscle_tpms = [data['muscle_TPM'] for data in self.protein_reference.values() if data['muscle_TPM'] > 0]
        if muscle_tpms:
            logger.info(f"\nMuscle expression statistics:")
            logger.info(f"  Mean TPM: {np.mean(muscle_tpms):.1f}")
            logger.info(f"  Median TPM: {np.median(muscle_tpms):.1f}")
            logger.info(f"  Max TPM: {np.max(muscle_tpms):.1f}")
            
        # Analyze subcellular localizations
        self._analyze_comprehensive_localization()
        
        # Analyze evidence levels
        self._analyze_evidence_levels()
        
        logger.info("="*70)
        
    def _analyze_comprehensive_localization(self):
        """
        Analyze subcellular localization patterns across all proteins
        """
        logger.info("\nSubcellular localization analysis:")
        
        # MitoCarta localizations
        mito_locs = {}
        for uniprot_id in self.mitochondrial_uniprot_ids:
            data = self.protein_reference[uniprot_id]
            loc = data['mitocarta_sub_localization']
            if loc and loc != '':
                for l in str(loc).split('|'):
                    l = l.strip()
                    mito_locs[l] = mito_locs.get(l, 0) + 1
                    
        logger.info("Top MitoCarta localizations:")
        for loc, count in sorted(mito_locs.items(), key=lambda x: x[1], reverse=True)[:5]:
            logger.info(f"  {loc}: {count} proteins")
            
        # HPA main localizations
        hpa_locs = {}
        for data in self.protein_reference.values():
            loc = data['main_localization']
            if loc and loc != '':
                for l in str(loc).split(','):
                    l = l.strip()
                    hpa_locs[l] = hpa_locs.get(l, 0) + 1
                    
        logger.info("Top HPA main localizations:")
        for loc, count in sorted(hpa_locs.items(), key=lambda x: x[1], reverse=True)[:5]:
            logger.info(f"  {loc}: {count} proteins")
            
    def _analyze_evidence_levels(self):
        """
        Analyze protein evidence levels from HPA
        """
        evidence_counts = {}
        for data in self.protein_reference.values():
            evidence = data['protein_evidence_level']
            if evidence:
                evidence_counts[evidence] = evidence_counts.get(evidence, 0) + 1
                
        logger.info("\nProtein evidence levels:")
        for evidence, count in sorted(evidence_counts.items(), key=lambda x: x[1], reverse=True):
            logger.info(f"  {evidence}: {count} proteins")
        
    def _analyze_mitochondrial_localization(self):
        """
        Analyze subcellular localization patterns in mitochondrial proteins
        """
        logger.info("\n" + "="*50)
        logger.info("MITOCHONDRIAL LOCALIZATION ANALYSIS")
        logger.info("="*50)
        
        # Count localization patterns
        localization_counts = {}
        
        for uniprot_id, data in self.mitochondrial_proteins.items():
            loc = data['sub_localization']
            if pd.isna(loc) or loc == '':
                loc = 'Unknown'
            
            # Handle multiple localizations
            if '|' in str(loc):
                locs = str(loc).split('|')
                for l in locs:
                    l = l.strip()
                    localization_counts[l] = localization_counts.get(l, 0) + 1
            else:
                localization_counts[str(loc)] = localization_counts.get(str(loc), 0) + 1
                
        # Sort by frequency
        sorted_locs = sorted(localization_counts.items(), key=lambda x: x[1], reverse=True)
        
        logger.info("Top subcellular localizations:")
        for loc, count in sorted_locs[:10]:
            percentage = count / len(self.mitochondrial_proteins) * 100
            logger.info(f"  {loc}: {count} proteins ({percentage:.1f}%)")
            
        logger.info("="*50)
        
    def _analyze_mitochondrial_pathways(self):
        """
        Analyze mitochondrial pathway patterns
        """
        logger.info("\n" + "="*50)
        logger.info("MITOCHONDRIAL PATHWAY ANALYSIS")
        logger.info("="*50)
        
        # Count pathway patterns
        pathway_counts = {}
        
        for uniprot_id, data in self.mitochondrial_proteins.items():
            pathways = data['pathways']
            if pd.isna(pathways) or pathways == '':
                pathways = 'Unknown'
                
            # Handle multiple pathways
            if '|' in str(pathways):
                paths = str(pathways).split('|')
                for p in paths:
                    p = p.strip()
                    pathway_counts[p] = pathway_counts.get(p, 0) + 1
            else:
                pathway_counts[str(pathways)] = pathway_counts.get(str(pathways), 0) + 1
                
        # Sort by frequency
        sorted_pathways = sorted(pathway_counts.items(), key=lambda x: x[1], reverse=True)
        
        logger.info("Top mitochondrial pathways:")
        for pathway, count in sorted_pathways[:10]:
            percentage = count / len(self.mitochondrial_proteins) * 100
            logger.info(f"  {pathway}: {count} proteins ({percentage:.1f}%)")
            
        logger.info("="*50)
        
    def is_mitochondrial(self, uniprot_id: str) -> bool:
        """
        Check if a UniProt ID corresponds to a mitochondrial protein
        """
        return uniprot_id in self.mitochondrial_uniprot_ids
        
    def is_muscle_expressed(self, uniprot_id: str) -> bool:
        """
        Check if a UniProt ID corresponds to a muscle-expressed protein
        """
        return uniprot_id in self.muscle_expressed_uniprot_ids
        
    def is_included(self, uniprot_id: str) -> bool:
        """
        Check if a UniProt ID should be included in the network
        (mitochondrial OR muscle-expressed)
        """
        return uniprot_id in self.included_uniprot_ids
        
    def get_protein_info(self, uniprot_id: str) -> Dict:
        """
        Get comprehensive protein annotation for a UniProt ID
        """
        return self.protein_reference.get(uniprot_id, {})
        
    def _load_file_completely(self, relative_path: str) -> pd.DataFrame:
        """
        Load complete file without row limit for processing
        """
        file_path = self.data_dir / relative_path
        
        try:
            if file_path.suffix == '.gz':
                if 'protein.links' in file_path.name or 'protein.physical' in file_path.name:
                    return pd.read_csv(file_path, sep=' ', low_memory=False)
                elif 'tab' in file_path.stem or 'txt' in file_path.stem:
                    return pd.read_csv(file_path, sep='\t', low_memory=False)
                else:
                    return pd.read_csv(file_path, low_memory=False)
                    
            elif file_path.suffix in ['.txt', '.tsv']:
                if 'UniProt2Reactome' in file_path.name:
                    df = pd.read_csv(file_path, sep='\t', low_memory=False, header=None)
                    df.columns = ['UniProt', 'Reactome_Pathway_ID', 'URL', 'Event_Name', 'Evidence_Code', 'Species']
                    return df
                else:
                    return pd.read_csv(file_path, sep='\t', low_memory=False)
                    
            elif file_path.suffix == '.xls' or file_path.suffix == '.xlsx':
                if 'MitoCarta3.0' in file_path.name:
                    return pd.read_excel(file_path, sheet_name='A Human MitoCarta3.0', engine='xlrd' if file_path.suffix == '.xls' else None)
                else:
                    return pd.read_excel(file_path, engine='xlrd' if file_path.suffix == '.xls' else None)
                    
            elif file_path.suffix == '.zip':
                with zipfile.ZipFile(file_path, 'r') as zip_ref:
                    txt_files = [f for f in zip_ref.namelist() if f.endswith('.txt') and 'tab3' in f]
                    if txt_files:
                        with zip_ref.open(txt_files[0]) as f:
                            return pd.read_csv(f, sep='\t', low_memory=False)
                            
        except Exception as e:
            logger.error(f"Error loading complete file {file_path}: {e}")
            return pd.DataFrame()
            
        return pd.DataFrame()

    def _load_file_chunked(self, relative_path: str, chunk_size: int = 50000):
        """
        Load file in chunks to reduce memory usage - returns an iterator
        """
        file_path = self.data_dir / relative_path
        logger.debug(f"Loading chunked file: {file_path}")
        
        try:
            if file_path.suffix == '.gz':
                if 'protein.links' in file_path.name or 'protein.physical' in file_path.name:
                    logger.debug(f"Using space separation for {file_path}")
                    return pd.read_csv(file_path, sep=' ', low_memory=False, chunksize=chunk_size)
                else:
                    # Default to tab separation for .gz files (most STRING files are tab-delimited)
                    logger.debug(f"Using tab separation for {file_path}")
                    return pd.read_csv(file_path, sep='\t', low_memory=False, chunksize=chunk_size)
                    
            elif file_path.suffix in ['.txt', '.tsv']:
                if 'UniProt2Reactome' in file_path.name:
                    chunk_iter = pd.read_csv(file_path, sep='\t', low_memory=False, header=None, chunksize=chunk_size)
                    for chunk in chunk_iter:
                        chunk.columns = ['UniProt', 'Reactome_Pathway_ID', 'URL', 'Event_Name', 'Evidence_Code', 'Species']
                        yield chunk
                    return
                else:
                    return pd.read_csv(file_path, sep='\t', low_memory=False, chunksize=chunk_size)
                    
            elif file_path.suffix == '.zip':
                with zipfile.ZipFile(file_path, 'r') as zip_ref:
                    txt_files = [f for f in zip_ref.namelist() if f.endswith('.txt') and 'tab3' in f]
                    if txt_files:
                        with zip_ref.open(txt_files[0]) as f:
                            return pd.read_csv(f, sep='\t', low_memory=False, chunksize=chunk_size)
                            
        except Exception as e:
            logger.error(f"Error loading chunked file {file_path}: {e}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            return iter([pd.DataFrame()])
            
        logger.warning(f"No matching file type for {file_path}")
        return iter([pd.DataFrame()])
        
    def _monitor_memory(self, phase_name: str):
        """
        Monitor memory usage and log current stats
        """
        try:
            process = psutil.Process(os.getpid())
            memory_info = process.memory_info()
            memory_mb = memory_info.rss / (1024 * 1024)
            
            # Get system memory info
            system_memory = psutil.virtual_memory()
            available_mb = system_memory.available / (1024 * 1024)
            
            logger.info(f"📊 Memory Monitor [{phase_name}]: "
                       f"Process={memory_mb:.1f}MB, "
                       f"Available={available_mb:.1f}MB, "
                       f"Usage={system_memory.percent:.1f}%")
            
            # Warning if memory usage is high
            if memory_mb > 2000:  # More than 2GB
                logger.warning(f"⚠️  High memory usage detected: {memory_mb:.1f}MB")
                gc.collect()  # Force garbage collection
                
        except Exception as e:
            logger.debug(f"Memory monitoring failed: {e}")
    
    def _optimize_memory(self):
        """
        Force garbage collection and memory optimization
        """
        gc.collect()
        logger.debug("🧹 Memory optimization: garbage collection completed")
        
    def _print_mapping_statistics(self):
        """
        Print comprehensive mapping statistics
        """
        logger.info("\n" + "="*60)
        logger.info("ID MAPPING STATISTICS")
        logger.info("="*60)
        
        stats = {
            'STRING→UniProt': len(self.id_mapping['string_to_uniprot']),
            'Symbol→UniProt': len(self.id_mapping['symbol_to_uniprot']),
            'UniProt→Symbol': len(self.id_mapping['uniprot_to_symbol']),
            'UniProt annotations': len(self.id_mapping['uniprot_info'])
        }
        
        for mapping_type, count in stats.items():
            logger.info(f"{mapping_type}: {count:,} entries")
            
        # Sample mappings for verification
        logger.info("\nSample mappings:")
        for i, (string_id, uniprot_id) in enumerate(list(self.id_mapping['string_to_uniprot'].items())[:3]):
            symbol = self.id_mapping['uniprot_to_symbol'].get(uniprot_id, 'N/A')
            logger.info(f"  {string_id} → {uniprot_id} ({symbol})")
            
        logger.info("="*60)
        
    def _create_mapping_functions(self):
        """
        Create helper functions for ID conversion
        """
        def map_to_uniprot(identifier: str, id_type: str = 'auto') -> Optional[str]:
            """Map any identifier to UniProt"""
            if id_type == 'auto':
                # Auto-detect ID type
                if identifier.startswith('9606.ENSP'):
                    return self.id_mapping['string_to_uniprot'].get(identifier)
                elif identifier.startswith('ENSP'):
                    # Try STRING format
                    return self.id_mapping['string_to_uniprot'].get(f'9606.{identifier}')
                else:
                    # Try as symbol
                    return self.id_mapping['symbol_to_uniprot'].get(identifier)
            else:
                mapping_dict = {
                    'string': self.id_mapping['string_to_uniprot'],
                    'symbol': self.id_mapping['symbol_to_uniprot'],
                    'ensembl': self.id_mapping['ensembl_to_uniprot'],
                    'entrez': self.id_mapping['entrez_to_uniprot']
                }
                return mapping_dict.get(id_type, {}).get(identifier)
                
        def get_symbol(uniprot_id: str) -> Optional[str]:
            """Get gene symbol for UniProt ID"""
            return self.id_mapping['uniprot_to_symbol'].get(uniprot_id)
            
        def get_protein_info(uniprot_id: str) -> Dict:
            """Get protein information for UniProt ID"""
            return self.id_mapping['uniprot_info'].get(uniprot_id, {})
            
        # Attach methods to class
        self.map_to_uniprot = map_to_uniprot
        self.get_symbol = get_symbol  
        self.get_protein_info = get_protein_info

def main():
    """
    Main execution function
    """
    logger.info("Starting Mitochondrial Network Integration Pipeline")
    
    # Initialize integrator
    integrator = MitoNetIntegrator()
    
    # Execute Phase 1: File Inspection
    integrator.phase1_inspect_files()
    
    # Execute Phase 2: ID Mapping
    integrator.phase2_build_id_mapping()
    
    # Execute Phase 3: Mitochondrial Reference
    integrator.phase3_create_mitochondrial_reference()
    
    # Execute Phase 4: Network Integration
    integrator.phase4_integrate_networks()
    
    # Execute Phase 5: Node Annotation
    integrator.phase5_annotate_nodes()
    
    # Execute Phase 6: Quality Control
    integrator.phase6_quality_control()
    
    # Execute Phase 7: Generate Exports
    integrator.phase7_generate_exports()
    
    logger.info("All phases complete! Mitochondrial network pipeline finished successfully.")

if __name__ == "__main__":
    main()
