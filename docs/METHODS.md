# Methods and workflow details

## Analytical scope

The project analyzes transcriptomic and lipidomic responses of *Mycobacterium tuberculosis* under:

- Exponential growth vs. growth arrest conditions.
- Iron-replete (+Fe) vs. iron-deprived (-Fe) environments.

## Implementation notes

- Figure-oriented R scripts are currently located in `Analyses/codes/`.
- The new standardized location for analysis code is `analysis/scripts/`.
- Reusable utilities should be placed in `analysis/functions/`.

## Execution strategy

- Use `Makefile` targets for discoverable and repeatable command execution.
- Keep generated outputs in `results/` and avoid committing regenerated files unnecessarily.
- Keep methods documentation in English and update this file when workflows change.

## Quality and traceability

- Use explicit file paths and deterministic script parameters.
- Save run logs to `results/logs/`.
- Record environment versions via `environment.yml` updates.
