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
	@echo "\033[1;34mInstallation & Setup:\033[0m"
	@grep -E '^(install|install-dev|sync|dev-setup|dev-reset):.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'
	@echo ""
	@echo "\033[1;34mDatabase Operations:\033[0m"
	@grep -E '^(db-.*|show-stats|checkpoints):.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'
	@echo ""
	@echo "\033[1;34mData Updates:\033[0m"
	@grep -E '^update-.*:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'
	@echo ""
	@echo "\033[1;34mNetwork Export & Filtering:\033[0m"
	@grep -E '^(export-.*):.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'
	@echo ""
	@echo "\033[1;34mTesting:\033[0m"
	@grep -E '^(test.*|quick-test):.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'
	@echo ""
	@echo "\033[1;34mDevelopment Workflows:\033[0m"
	@grep -E '^(dev-test-data|demo|full-pipeline):.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'
	@echo ""
	@echo "\033[1;34mLegacy Pipeline:\033[0m"
	@grep -E '^run-.*:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'
	@echo ""
	@echo "\033[1;34mMaintenance:\033[0m"
	@grep -E '^(lint|format|clean.*):.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

# Installation and setup
.PHONY: install install-dev sync
install: ## Install the package and dependencies
	uv sync

install-dev: ## Install with development dependencies
	uv sync --extra test

sync: ## Sync dependencies (alias for install)
	uv sync

# Database operations
.PHONY: db-init db-status
db-init: ## Initialize the database
	$(PYTHON) -m mitonet.cli init

db-status: ## Show database status
	$(PYTHON) -m mitonet.cli status

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

# Network export and filtering
.PHONY: export-all export-mitochondrial export-muscle export-genes export-high-confidence export-predefined
export-all: ## Export complete network (all proteins and interactions)
	$(PYTHON) -m mitonet.cli export-network --filter-type all

export-mitochondrial: ## Export mitochondrial protein network
	$(PYTHON) -m mitonet.cli export-network --filter-type mitochondrial

export-muscle: ## Export muscle-expressed protein network
	$(PYTHON) -m mitonet.cli export-network --filter-type muscle

export-genes: ## Export network for specific genes (use GENES="ATP1A1,MYOD1" NEIGHBORS=1)
	@if [ -n "$(GENES)" ]; then \
		$(PYTHON) -m mitonet.cli export-network --filter-type genes --genes "$(GENES)" --neighbors $(NEIGHBORS); \
	elif [ -n "$(UNIPROTS)" ]; then \
		$(PYTHON) -m mitonet.cli export-network --filter-type genes --uniprots "$(UNIPROTS)" --neighbors $(NEIGHBORS); \
	else \
		echo "Usage: make export-genes GENES='ATP1A1,MYOD1' NEIGHBORS=1"; \
		echo "   or: make export-genes UNIPROTS='P12345,Q67890' NEIGHBORS=1"; \
	fi

export-high-confidence: ## Export high-confidence interactions only (use CONFIDENCE=0.7)
	$(PYTHON) -m mitonet.cli export-network --filter-type high_confidence --min-confidence $(or $(CONFIDENCE),0.7)

export-predefined: ## Export predefined networks (mitochondrial, muscle, high-confidence)
	$(PYTHON) -m mitonet.cli export-predefined

# Custom network export
export-custom: ## Export custom filtered network (use FILTER_TYPE, MIN_CONF, etc.)
	@if [ -n "$(FILTER_TYPE)" ]; then \
		$(PYTHON) -m mitonet.cli export-network --filter-type $(FILTER_TYPE) \
			$(if $(MIN_CONF),--min-confidence $(MIN_CONF)) \
			$(if $(MAX_CONF),--max-confidence $(MAX_CONF)) \
			$(if $(GENES),--genes "$(GENES)") \
			$(if $(NEIGHBORS),--neighbors $(NEIGHBORS)) \
			$(if $(FORMATS),--formats "$(FORMATS)") \
			$(if $(OUTPUT),--output-prefix $(OUTPUT)); \
	else \
		echo "Usage: make export-custom FILTER_TYPE=mitochondrial MIN_CONF=0.5"; \
		echo "Available FILTER_TYPE: mitochondrial, muscle, genes, high_confidence, all"; \
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
	@echo "Try: make update-all && make export-predefined"

dev-reset: clean db-init ## Reset development environment
	@echo "Development environment reset!"

dev-test-data: ## Add sample data for development (requires data files)
	@echo "Loading sample data sources..."
	@if [ -f "networks/string/9606.protein.aliases.v12.0.txt.gz" ]; then \
		$(PYTHON) -m mitonet.cli update --source STRING_aliases; \
		echo "✅ STRING aliases loaded"; \
	else \
		echo "⚠️  STRING data not found - download to networks/string/"; \
	fi
	@if [ -f "networks/mitocarta/Human.MitoCarta3.0.xls" ]; then \
		$(PYTHON) -m mitonet.cli update --source MitoCarta; \
		echo "✅ MitoCarta data loaded"; \
	else \
		echo "⚠️  MitoCarta data not found - download to networks/mitocarta/"; \
	fi
	@echo "Sample data loading attempted - check database status with 'make db-status'"

# Quick workflows
.PHONY: quick-test full-pipeline demo
quick-test: ## Run quick tests for development
	$(PYTEST) $(TEST_DIR)/unit/test_database.py::TestMitoNetDatabase::test_initialize_database -v
	$(PYTEST) $(TEST_DIR)/unit/test_cli.py::TestCLICommands::test_help_command -v

full-pipeline: dev-setup dev-test-data update-all export-predefined ## Run complete pipeline end-to-end
	@echo "Full pipeline completed! Check outputs/ directory for networks."

demo: dev-setup dev-test-data ## Setup demo environment with sample data
	@echo "Demo environment ready!"
	@echo "Try: make show-stats"
	@echo "Then: make export-mitochondrial or make export-muscle"