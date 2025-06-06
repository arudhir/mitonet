#!/usr/bin/env python3
"""
Memory-optimized version of the Mitochondrial Network Integration Pipeline
This version includes chunked processing and memory monitoring to avoid heap overflow
"""

import sys
import time
import logging
import gc
import os
from pathlib import Path
from main import MitoNetIntegrator

import psutil

# Setup enhanced logging with memory info
class MemoryFormatter(logging.Formatter):
    def format(self, record):
        try:
            process = psutil.Process(os.getpid())
            memory_mb = process.memory_info().rss / (1024 * 1024)
            record.memory = f"[{memory_mb:.0f}MB]"
        except:
            record.memory = "[MEM?]"
        return super().format(record)

# Configure logging
memory_formatter = MemoryFormatter('%(asctime)s %(memory)s - %(levelname)s - %(message)s')

# File handler
file_handler = logging.FileHandler('outputs/pipeline_memory_optimized.log')
file_handler.setFormatter(memory_formatter)

# Console handler
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(memory_formatter)

# Setup logger
logging.basicConfig(
    level=logging.INFO,
    handlers=[file_handler, console_handler]
)
logger = logging.getLogger(__name__)

def force_garbage_collection():
    """Force aggressive garbage collection"""
    collected = gc.collect()
    logger.debug(f"🧹 Garbage collection: {collected} objects collected")

def run_memory_optimized_pipeline():
    """Run the pipeline with memory optimizations"""
    start_time = time.time()
    
    logger.info("="*80)
    logger.info("MEMORY-OPTIMIZED MITOCHONDRIAL NETWORK INTEGRATION PIPELINE")
    logger.info("="*80)
    logger.info(f"Starting at: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Monitor initial system state
    memory = psutil.virtual_memory()
    logger.info(f"💾 System Memory: {memory.total/(1024**3):.1f}GB total, "
                f"{memory.available/(1024**3):.1f}GB available ({memory.percent:.1f}% used)")
    
    # Ensure outputs directory exists
    Path("outputs").mkdir(exist_ok=True)
    
    try:
        # Initialize integrator
        logger.info("🚀 Initializing integrator...")
        integrator = MitoNetIntegrator()
        force_garbage_collection()
        
        # Phase 1: File Inspection (fast, low memory)
        phase_start = time.time()
        logger.info("\n🔍 Executing Phase 1: Data File Inspection...")
        integrator.phase1_inspect_files()
        logger.info(f"✅ Phase 1 complete ({time.time() - phase_start:.1f}s)")
        force_garbage_collection()
        
        # Phase 2: ID Mapping (moderate, chunked)
        phase_start = time.time()
        logger.info("\n🔗 Executing Phase 2: Identifier Standardization (chunked)...")
        integrator.phase2_build_id_mapping()
        logger.info(f"✅ Phase 2 complete ({time.time() - phase_start:.1f}s)")
        force_garbage_collection()
        
        # Phase 3: Protein Reference (moderate)
        phase_start = time.time()
        logger.info("\n🧬 Executing Phase 3: Protein Reference Creation...")
        integrator.phase3_create_mitochondrial_reference()
        logger.info(f"✅ Phase 3 complete ({time.time() - phase_start:.1f}s)")
        force_garbage_collection()
        
        # Phase 4: Network Integration (memory-intensive, now chunked)
        phase_start = time.time()
        logger.info("\n🕸️  Executing Phase 4: Network Integration (chunked processing)...")
        integrator.phase4_integrate_networks()
        logger.info(f"✅ Phase 4 complete ({time.time() - phase_start:.1f}s)")
        force_garbage_collection()
        
        # Phase 5: Node Annotation (moderate)
        phase_start = time.time()
        logger.info("\n🏷️  Executing Phase 5: Node Annotation...")
        integrator.phase5_annotate_nodes()
        logger.info(f"✅ Phase 5 complete ({time.time() - phase_start:.1f}s)")
        force_garbage_collection()
        
        # Phase 6: Quality Control (fast)
        phase_start = time.time()
        logger.info("\n✅ Executing Phase 6: Quality Control...")
        integrator.phase6_quality_control()
        logger.info(f"✅ Phase 6 complete ({time.time() - phase_start:.1f}s)")
        force_garbage_collection()
        
        # Phase 7: Export Generation (fast)
        phase_start = time.time()
        logger.info("\n📤 Executing Phase 7: Export Generation...")
        integrator.phase7_generate_exports()
        logger.info(f"✅ Phase 7 complete ({time.time() - phase_start:.1f}s)")
        force_garbage_collection()
        
        # Summary
        total_time = time.time() - start_time
        logger.info("\n" + "="*80)
        logger.info("🎉 MEMORY-OPTIMIZED PIPELINE COMPLETED SUCCESSFULLY!")
        logger.info("="*80)
        logger.info(f"Total execution time: {total_time/60:.1f} minutes")
        logger.info(f"Network size: {integrator.network.number_of_nodes():,} nodes, {integrator.network.number_of_edges():,} edges")
        
        # Final memory check
        final_memory = psutil.virtual_memory()
        logger.info(f"💾 Final Memory Usage: {final_memory.percent:.1f}% ({final_memory.available/(1024**3):.1f}GB available)")
        
        # List generated files
        logger.info("\n📁 Generated files in outputs/ directory:")
        outputs_dir = Path("outputs")
        for file_path in sorted(outputs_dir.glob("*")):
            if file_path.is_file() and file_path.name != "pipeline_memory_optimized.log":
                size_mb = file_path.stat().st_size / (1024*1024)
                logger.info(f"  📄 {file_path.name} ({size_mb:.1f} MB)")
        
        logger.info("\n🚀 Ready for analysis! Check outputs/mitonet_summary_report.md for details.")
        
        return True
        
    except Exception as e:
        logger.error(f"❌ Pipeline failed: {e}")
        logger.exception("Full error traceback:")
        return False

if __name__ == "__main__":
    logger.info("🔧 Memory-optimized pipeline starting...")
    logger.info("📊 This version uses chunked processing to reduce memory usage")
    success = run_memory_optimized_pipeline()
    sys.exit(0 if success else 1)