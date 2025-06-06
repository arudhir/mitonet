#!/usr/bin/env python3
"""
Test runner script for MitoNet with different test categories
"""

import subprocess
import sys
from pathlib import Path

def run_tests(test_type="basic"):
    """Run tests based on type"""
    
    base_cmd = ["uv", "run", "pytest"]
    
    if test_type == "basic":
        # Run basic tests that are known to work
        cmd = base_cmd + [
            "tests/unit/test_database.py::TestMitoNetDatabase::test_initialize_database",
            "tests/unit/test_database.py::TestMitoNetDatabase::test_save_and_load_checkpoint",
            "tests/unit/test_cli.py::TestCLICommands::test_help_command",
            "-v"
        ]
    elif test_type == "unit":
        # Run all unit tests
        cmd = base_cmd + ["tests/unit/", "-v"]
    elif test_type == "integration":
        # Run integration tests
        cmd = base_cmd + ["tests/integration/", "-v", "-m", "integration"]
    elif test_type == "slow":
        # Run slow tests
        cmd = base_cmd + ["tests/", "-v", "-m", "slow"]
    elif test_type == "all":
        # Run all tests
        cmd = base_cmd + ["tests/", "-v"]
    elif test_type == "coverage":
        # Run with coverage report
        cmd = base_cmd + [
            "tests/unit/test_database.py::TestMitoNetDatabase::test_initialize_database",
            "tests/unit/test_database.py::TestMitoNetDatabase::test_save_and_load_checkpoint", 
            "tests/unit/test_cli.py::TestCLICommands::test_help_command",
            "--cov=mitonet",
            "--cov-report=term-missing",
            "--cov-report=html"
        ]
    else:
        print(f"Unknown test type: {test_type}")
        print("Available types: basic, unit, integration, slow, all, coverage")
        return 1
    
    print(f"Running {test_type} tests...")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=False)
        return result.returncode
    except KeyboardInterrupt:
        print("\nTests interrupted by user")
        return 1
    except Exception as e:
        print(f"Error running tests: {e}")
        return 1

if __name__ == "__main__":
    test_type = sys.argv[1] if len(sys.argv) > 1 else "basic"
    exit_code = run_tests(test_type)
    sys.exit(exit_code)