# MitoNet Pipeline Makefile
# Provides common development and deployment commands

# Variables
PYTHON := uv run python
PIP := uv pip
PYTEST := uv run pytest
PACKAGE := mitonet
TEST_DIR := tests
COV_DIR := htmlcov

# Default target
.PHONY: help
help: ## Show this help message
	@echo "MitoNet Pipeline - Available Commands:"
	@echo ""
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

# Installation and setup
.PHONY: install install-dev sync
install: ## Install the package and dependencies
	uv sync

install-dev: ## Install with development dependencies
	uv sync --extra test

sync: ## Sync dependencies (alias for install)
	uv sync

# Database operations
.PHONY: db-init db-status db-add-genes
db-init: ## Initialize the database
	$(PYTHON) -m mitonet.cli init

db-status: ## Show database status
	$(PYTHON) -m mitonet.cli status

db-add-genes: ## Add genes to database (use GENES="gene1,gene2" or UNIPROTS="P12345,Q67890")
	@if [ -n "$(GENES)" ]; then \
		$(PYTHON) -m mitonet.cli add-genes --genes "$(GENES)"; \
	elif [ -n "$(UNIPROTS)" ]; then \
		$(PYTHON) -m mitonet.cli add-genes --uniprots "$(UNIPROTS)"; \
	else \
		echo "Usage: make db-add-genes GENES='ATP1A1,MYOD1' or UNIPROTS='P12345,Q67890'"; \
	fi

# Data ingestion
.PHONY: update-all update-string update-biogrid update-mitocarta update-hpa
update-all: ## Update all data sources
	$(PYTHON) -m mitonet.cli update

update-string: ## Update STRING data
	$(PYTHON) -m mitonet.cli update --source STRING_aliases
	$(PYTHON) -m mitonet.cli update --source STRING_interactions

update-biogrid: ## Update BioGRID data
	$(PYTHON) -m mitonet.cli update --source BIOGRID

update-mitocarta: ## Update MitoCarta data
	$(PYTHON) -m mitonet.cli update --source MitoCarta

update-hpa: ## Update HPA muscle data
	$(PYTHON) -m mitonet.cli update --source HPA_muscle

# Network generation
.PHONY: generate-network export-network
generate-network: ## Generate network from current database
	$(PYTHON) -m mitonet.cli export-network --format all

export-network: ## Export network in specific format (use FORMAT=json|graphml|csv)
	@if [ -n "$(FORMAT)" ]; then \
		$(PYTHON) -m mitonet.cli export-network --format $(FORMAT); \
	else \
		echo "Usage: make export-network FORMAT=json"; \
		echo "Available formats: json, graphml, csv, all"; \
	fi

# Legacy pipeline
.PHONY: run-legacy run-memory-optimized run-simplified
run-legacy: ## Run the original main.py pipeline
	$(PYTHON) main.py

run-memory-optimized: ## Run memory-optimized pipeline
	$(PYTHON) run_memory_optimized.py

run-simplified: ## Run simplified pipeline
	$(PYTHON) run_simplified.py

# Testing
.PHONY: test test-unit test-integration test-cli test-database test-slow test-coverage
test: ## Run all tests
	$(PYTEST)

test-unit: ## Run unit tests only
	$(PYTEST) $(TEST_DIR)/unit/ -v

test-integration: ## Run integration tests only
	$(PYTEST) $(TEST_DIR)/integration/ -v -m integration

test-cli: ## Run CLI tests only
	$(PYTEST) -m cli -v

test-database: ## Run database tests only
	$(PYTEST) -m database -v

test-slow: ## Run slow tests (including performance tests)
	$(PYTEST) -m slow -v

test-coverage: ## Run tests with detailed coverage report
	$(PYTEST) --cov-report=html --cov-report=term-missing
	@echo "Coverage report available at: $(COV_DIR)/index.html"

# Code quality and maintenance
.PHONY: lint format clean clean-cache clean-outputs clean-all
lint: ## Run linting (if available)
	@echo "Note: Add linting tools (ruff, black, etc.) to pyproject.toml for this command"

format: ## Format code (if available)
	@echo "Note: Add formatting tools (black, etc.) to pyproject.toml for this command"

clean: ## Clean up temporary files and caches
	rm -rf __pycache__/ .pytest_cache/ *.pyc .coverage
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete

clean-cache: clean ## Clean caches (alias for clean)

clean-outputs: ## Clean output files
	rm -rf $(COV_DIR)/ coverage.xml
	rm -rf outputs/*.log

clean-all: clean clean-outputs ## Clean everything
	rm -rf .pytest_cache/ htmlcov/ coverage.xml

# Documentation and reporting
.PHONY: checkpoints show-stats
checkpoints: ## Show processing checkpoints
	$(PYTHON) -m mitonet.cli checkpoints

show-stats: ## Show database statistics
	$(PYTHON) -m mitonet.cli status

# Development workflow
.PHONY: dev-setup dev-reset dev-test-data
dev-setup: install-dev db-init ## Setup development environment
	@echo "Development environment ready!"
	@echo "Try: make db-add-genes GENES='ATP1A1,MYOD1,CYC1'"

dev-reset: clean db-init ## Reset development environment
	@echo "Development environment reset!"

dev-test-data: ## Add test data for development
	$(PYTHON) -m mitonet.cli add-genes --genes "ATP1A1,MYOD1,CYC1,NDUFA1,ACTA1"
	@echo "Test data added to database"

# Quick workflows
.PHONY: quick-test full-pipeline demo
quick-test: ## Run quick tests for development
	$(PYTEST) $(TEST_DIR)/unit/test_database.py::TestMitoNetDatabase::test_initialize_database -v
	$(PYTEST) $(TEST_DIR)/unit/test_cli.py::TestCLICommands::test_help_command -v

full-pipeline: dev-setup dev-test-data update-all generate-network ## Run complete pipeline end-to-end
	@echo "Full pipeline completed!"

demo: dev-setup ## Setup demo environment with sample data
	$(PYTHON) -m mitonet.cli add-genes --genes "ATP1A1,MYOD1,CYC1,NDUFA1,ACTA1,TPM1,CYCS"
	@echo "Demo environment ready with sample proteins!"
	@echo "Try: make show-stats"