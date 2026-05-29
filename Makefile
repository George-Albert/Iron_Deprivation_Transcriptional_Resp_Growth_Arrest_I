.PHONY: help validate run-legacy-scripts clean

help:
	@echo "Available targets:"
	@echo "  make validate           - Check that key project paths exist"
	@echo "  make run-legacy-scripts - Print and execute legacy figure scripts"
	@echo "  make clean              - Remove temporary R session artifacts"

validate:
	@test -d Analyses/codes
	@test -d Analyses/inputs
	@test -d Analyses/outputs
	@echo "Validation passed: legacy analysis directories are present."

run-legacy-scripts:
	@echo "Running legacy figure scripts from Analyses/codes"
	@for f in Analyses/codes/*.R; do \
		echo "Rscript $$f"; \
		Rscript $$f || exit 1; \
	done

clean:
	@find . -type f \( -name '.Rhistory' -o -name '.RData' \) -delete
	@echo "Removed R session artifacts."
