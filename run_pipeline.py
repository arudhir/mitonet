#!/usr/bin/env python3
"""
Run the complete Mitochondrial Network Integration Pipeline
"""

import sys
import time
import logging
from pathlib import Path
from main import MitoNetIntegrator

# Setup enhanced logging
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('outputs/pipeline.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def run_complete_pipeline():
    """Run all 7 phases of the pipeline with progress tracking"""
    start_time = time.time()
    
    logger.info("="*80)
    logger.info("MITOCHONDRIAL NETWORK INTEGRATION PIPELINE")
    logger.info("="*80)
    logger.info(f"Starting at: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Ensure outputs directory exists
    Path("outputs").mkdir(exist_ok=True)
    
    try:
        # Initialize integrator
        logger.info("Initializing integrator...")
        integrator = MitoNetIntegrator()
        
        # Phase 1: File Inspection (fast)
        phase_start = time.time()
        logger.info("\n🔍 Executing Phase 1: Data File Inspection...")
        integrator.phase1_inspect_files()
        logger.info(f"✅ Phase 1 complete ({time.time() - phase_start:.1f}s)")
        
        # Phase 2: ID Mapping (moderate)
        phase_start = time.time()
        logger.info("\n🔗 Executing Phase 2: Identifier Standardization...")
        integrator.phase2_build_id_mapping()
        logger.info(f"✅ Phase 2 complete ({time.time() - phase_start:.1f}s)")
        
        # Phase 3: Protein Reference (moderate)
        phase_start = time.time()
        logger.info("\n🧬 Executing Phase 3: Protein Reference Creation...")
        integrator.phase3_create_mitochondrial_reference()
        logger.info(f"✅ Phase 3 complete ({time.time() - phase_start:.1f}s)")
        
        # Phase 4: Network Integration (slow - largest phase)
        phase_start = time.time()
        logger.info("\n🕸️  Executing Phase 4: Network Integration (may take several minutes)...")
        integrator.phase4_integrate_networks()
        logger.info(f"✅ Phase 4 complete ({time.time() - phase_start:.1f}s)")
        
        # Phase 5: Node Annotation (moderate)
        phase_start = time.time()
        logger.info("\n🏷️  Executing Phase 5: Node Annotation...")
        integrator.phase5_annotate_nodes()
        logger.info(f"✅ Phase 5 complete ({time.time() - phase_start:.1f}s)")
        
        # Phase 6: Quality Control (fast)
        phase_start = time.time()
        logger.info("\n✅ Executing Phase 6: Quality Control...")
        integrator.phase6_quality_control()
        logger.info(f"✅ Phase 6 complete ({time.time() - phase_start:.1f}s)")
        
        # Phase 7: Export Generation (fast)
        phase_start = time.time()
        logger.info("\n📤 Executing Phase 7: Export Generation...")
        integrator.phase7_generate_exports()
        logger.info(f"✅ Phase 7 complete ({time.time() - phase_start:.1f}s)")
        
        # Summary
        total_time = time.time() - start_time
        logger.info("\n" + "="*80)
        logger.info("🎉 PIPELINE COMPLETED SUCCESSFULLY!")
        logger.info("="*80)
        logger.info(f"Total execution time: {total_time/60:.1f} minutes")
        logger.info(f"Network size: {integrator.network.number_of_nodes():,} nodes, {integrator.network.number_of_edges():,} edges")
        
        # List generated files
        logger.info("\n📁 Generated files in outputs/ directory:")
        outputs_dir = Path("outputs")
        for file_path in sorted(outputs_dir.glob("*")):
            if file_path.is_file():
                size_mb = file_path.stat().st_size / (1024*1024)
                logger.info(f"  📄 {file_path.name} ({size_mb:.1f} MB)")
        
        logger.info("\n🚀 Ready for analysis! Check outputs/mitonet_summary_report.md for details.")
        
        return True
        
    except Exception as e:
        logger.error(f"❌ Pipeline failed: {e}")
        logger.exception("Full error traceback:")
        return False

if __name__ == "__main__":
    success = run_complete_pipeline()
    sys.exit(0 if success else 1)