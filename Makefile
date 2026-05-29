.PHONY: help validate run-legacy-scripts clean install-env

help:
	@echo "Available targets:"
	@echo "  make install-env       - Create conda environment from environment.yml"
	@echo "  make validate          - Validate repository structure"
	@echo "  make run-legacy-scripts - Execute all R scripts in Analyses/codes/"
	@echo "  make clean             - Remove R session artifacts"
	@echo "  make help              - Display this help message"

install-env:
	@echo "Creating conda environment..."
	conda env create -f environment.yml
	@echo "Environment created. Activate with: conda activate mtb-iron-deprivation"

validate:
	@echo "Validating repository structure..."
	@test -d Analyses/codes || (echo "ERROR: Analyses/codes not found"; exit 1)
	@test -d Analyses/inputs || (echo "ERROR: Analyses/inputs not found"; exit 1)
	@test -d Analyses/outputs || (echo "ERROR: Analyses/outputs not found"; exit 1)
	@echo "✓ Repository structure validation passed."

run-legacy-scripts: validate
	@echo "Executing R scripts from Analyses/codes/..."
	@for script in Analyses/codes/*.R; do \
		if [ -f "$$script" ]; then \
			echo "Running: $$script"; \
			Rscript "$$script" || exit 1; \
		fi; \
	done
	@echo "✓ All scripts executed successfully."

clean:
	@echo "Cleaning up R session artifacts..."
	@rm -f .Rhistory .RData Rplots.pdf
	@find . -name "*.Rout" -delete
	@echo "✓ Cleanup complete."
