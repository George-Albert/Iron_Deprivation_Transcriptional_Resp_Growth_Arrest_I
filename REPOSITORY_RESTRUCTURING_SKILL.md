# Repository Restructuring Skill for Reproducible Research

**Version**: 1.0  
**Last Updated**: 2026-05-29  
**Status**: Production Ready

## Overview

This skill transforms disorganized repositories into structured, reproducible, and professional projects through a hierarchical 5-phase approach. It's context-aware, generating files and documentation tailored to project type, primary language, and existing structure.

---

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Phase 1: Discovery & Analysis](#phase-1-discovery--analysis)
3. [Phase 2: Planning & Design](#phase-2-planning--design)
4. [Phase 3: Content Generation](#phase-3-content-generation)
5. [Phase 4: Validation & Cleanup](#phase-4-validation--cleanup)
6. [Phase 5: Delivery & Documentation](#phase-5-delivery--documentation)
7. [Conditional Logic & Templates](#conditional-logic--templates)
8. [Edge Cases & Troubleshooting](#edge-cases--troubleshooting)

---

## Architecture Overview

```
┌─────────────────────────────────────────────────────┐
│ LEVEL 0: GENERAL OBJECTIVE                          │
│ "Transform → Reproducible, Scalable, Professional"  │
└──────────────┬──────────────────────────────────────┘
               │
    ┌──────────┴──────────┐
    │                     │
┌───▼────────────┐  ┌────▼──────────┐
│ LEVEL 1:       │  │ LEVEL 1:      │
│ 5 Phases       │  │ Key Decisions │
├────────────────┤  ├───────────────┤
│ 1. Discovery   │  │ project_type  │
│ 2. Planning    │  │ language      │
│ 3. Generation  │  │ preserve_legacy
│ 4. Validation  │  │ has_paper     │
│ 5. Delivery    │  │ has_data      │
└────┬───────────┘  └───────────────┘
     │
┌────▼──────────────────────────────┐
│ LEVEL 2: CONDITIONAL LOGIC        │
│ "If type=research then..."        │
└────┬──────────────────────────────┘
     │
┌────▼──────────────────────────────┐
│ LEVEL 3: CONTEXT-SPECIFIC FILES   │
│ Templates with injected context   │
├───────────────────────────────────┤
│ - .gitignore                      │
│ - CITATION.cff                    │
│ - environment.yml                 │
│ - Makefile                        │
│ - README.md                       │
│ - docs/*.md (4+ files)            │
│ - .github/ISSUE_TEMPLATE/         │
└────┬──────────────────────────────┘
     │
┌────▼──────────────────────────────┐
│ LEVEL 4: VALIDATION & RULES       │
│ "If R then validate .Rhistory..." │
└───────────────────────────────────┘
```

---

## Phase 1: Discovery & Analysis

### Inputs

```yaml
REQUIRED_INPUTS:
  repository_name:
    description: "Full repository name (with owner if needed)"
    example: "Iron_Deprivation_Transcriptional_Resp_Growth_Arrest_I"
  
  repository_owner:
    description: "GitHub username or organization"
    example: "George-Albert"
  
  project_type:
    description: "Type of project"
    enum:
      - research         # Data science, scientific papers, experiments
      - library          # Reusable code, SDK, toolkit
      - webapp           # Web application (frontend/backend)
      - cli              # Command-line tool
      - documentation    # Documentation-focused
      - monorepo         # Multiple projects in one repo
  
  primary_language:
    description: "Main programming language"
    enum:
      - R
      - Python
      - JavaScript
      - TypeScript
      - Java
      - Go
      - Rust
      - C++
      - Mixed
  
  brief_description:
    description: "1-2 sentence summary of the project"
    example: "RNA-seq analysis of M. tuberculosis transcriptional responses to iron deprivation"

OPTIONAL_INPUTS:
  has_paper:
    description: "Is this research backed by a published paper?"
    type: boolean
    default: false
  
  has_data:
    description: "Does the repo contain significant data files?"
    type: boolean
    default: false
  
  preserve_legacy:
    description: "Keep existing directory structure alongside new one?"
    type: boolean
    default: true
  
  authors:
    description: "List of project authors for CITATION.cff"
    type: array[string]
    example: ["Alebouyeh S", "Cárdenas-Pestana JA"]
  
  doi:
    description: "DOI of published paper (if applicable)"
    type: string
    example: "10.3389/fmicb.2022.956602"
  
  paper_title:
    description: "Full title of published paper"
    type: string
  
  paper_year:
    description: "Year of publication"
    type: integer
    example: 2022
  
  paper_journal:
    description: "Journal name"
    type: string
    example: "Frontiers in Microbiology"
```

### Automatic Analysis

Analyze the repository to understand its current state:

```yaml
DISCOVERY_TASKS:
  1. Detect Current Structure
     - List existing directories: /Analyses, /data, /src, /tests, etc.
     - Identify key files: README, LICENSE, package.json, environment.yml, etc.
     - Assess organization level: 0% (chaotic) to 100% (well-organized)
  
  2. Language Detection
     - Scan file extensions: .R, .py, .js, .ts, .java, etc.
     - Check package files: package.json, setup.py, DESCRIPTION, etc.
     - Confirm with user input if ambiguous
  
  3. Risk Assessment
     - Check for OS artifacts: .DS_Store, Thumbs.db, .AppleDouble
     - Check for generated files: __pycache__, node_modules, dist/, venv/
     - Assess data size: large files that shouldn't be in git
     - Check git history: any sensitive files already committed
  
  4. Preservation Check
     - Legacy directories to keep: Analyses/, legacy/, etc.
     - Existing documentation: README.md, docs/
     - Configuration files: .gitignore, package.json, setup.py
  
  5. Metadata Extraction (if has_paper=true)
     - Parse DOI from README or look up metadata
     - Extract authors, journal, publication date
     - Validate DOI format and accessibility

OUTPUT:
  Repository Profile:
    current_structure: {}
    detected_languages: []
    organization_score: 0-100
    issues_found: []
    recommended_cleanup: []
    preservation_items: []
    metadata_ready: boolean
```

---

## Phase 2: Planning & Design

### Key Decisions & Conditional Logic

```yaml
DECISION_TREE:
  
  IF project_type == "research" THEN:
    directories:
      - data/raw/              # Original, unmodified data
      - data/processed/        # Cleaned, processed data
      - analysis/scripts/      # Main analysis code
      - analysis/functions/    # Reusable functions
      - analysis/notebooks/    # RMarkdown, Jupyter notebooks
      - results/figures/       # Generated plots
      - results/tables/        # Output tables
      - results/logs/          # Execution logs
      - manuscript/main/       # Main paper figures
      - manuscript/supplementary/
      - manuscript/proofs/
      - docs/                  # Full documentation
    
    documentation_files:
      - docs/README.md         # Documentation index
      - docs/REPRODUCIBILITY.md    # How to run analyses
      - docs/DATA_SOURCES.md       # Data descriptions
      - docs/METHODS.md            # Technical methodology
    
    config_files:
      - environment.yml or Pipfile or requirements.txt
      - Makefile               # With validate, run-analyses targets
      - .gitignore             # Data + analysis outputs excluded
    
    github_templates:
      - ISSUE_TEMPLATE/reproducibility_report.md
    
    key_features:
      - Focus on reproducibility
      - Data management best practices
      - Methodology documentation
      - Validation automation

  ELSE IF project_type == "library" THEN:
    directories:
      - src/ or lib/           # Source code
      - tests/                 # Unit/integration tests
      - examples/              # Usage examples
      - docs/                  # API documentation
      - benchmarks/ (optional) # Performance tests
    
    documentation_files:
      - docs/README.md         # API overview
      - docs/INSTALLATION.md   # How to install
      - docs/GETTING_STARTED.md
      - docs/API.md            # Full API reference
      - docs/CONTRIBUTING.md   # How to contribute
    
    config_files:
      - package.json (JS/TS) or setup.py (Python) or Cargo.toml (Rust)
      - pytest.ini / jest.config.js / cargo.toml
      - .github/workflows/ (CI/CD)
    
    github_templates:
      - ISSUE_TEMPLATE/bug_report.md
      - ISSUE_TEMPLATE/feature_request.md
      - pull_request_template.md
    
    key_features:
      - Code quality emphasis
      - Testing infrastructure
      - API clarity
      - Contribution guidelines

  ELSE IF project_type == "webapp" THEN:
    directories:
      - frontend/              # React, Vue, Angular, etc.
      - backend/               # Node, Python, Java, etc.
      - infrastructure/        # Docker, Kubernetes configs
      - docs/                  # Architecture docs
      - scripts/               # Deployment scripts
    
    documentation_files:
      - docs/ARCHITECTURE.md   # System design
      - docs/DEPLOYMENT.md     # How to deploy
      - docs/API.md            # Backend API docs
      - docs/DEVELOPMENT.md    # Local dev setup
    
    config_files:
      - docker-compose.yml
      - .env.example
      - .github/workflows/ (CI/CD with deploy)
    
    github_templates:
      - pull_request_template.md
    
    key_features:
      - Deployment automation
      - Environment management
      - Architecture documentation
      - CI/CD pipelines

  IF language == "R" THEN:
    environment_file: environment.yml
    environment_content:
      - r-base=4.3.0
      - r-tidyverse, r-ggplot2, r-dplyr, etc.
      - r-rmarkdown, r-knitr
      - r-devtools (if library)
      - bioconductor packages (if research)
    
    makefile_targets:
      - install-env: "conda env create -f environment.yml"
      - validate: "test -d <required dirs>"
      - run-analyses: "for script in Analyses/codes/*.R do Rscript"
      - clean: "rm .Rhistory .RData *.Rout"
    
    gitignore_additions:
      - .Rhistory, .Rapp.history, .RData, .Ruserdata
      - .Rproj.user, *.Rproj
      - Rplots.pdf

  ELSE IF language == "Python" THEN:
    environment_file: requirements.txt or pyproject.toml
    environment_content:
      - numpy, pandas, scipy, scikit-learn
      - matplotlib, seaborn, plotly (if visualization-heavy)
      - pytest, coverage (if library)
      - jupyter (if notebooks used)
    
    makefile_targets:
      - install-env: "pip install -r requirements.txt"
      - test: "pytest tests/"
      - lint: "flake8 src/ tests/"
      - run-analyses: "python analysis/scripts/*.py"
    
    gitignore_additions:
      - __pycache__/, *.pyc, *.pyo
      - venv/, env/, .venv/
      - .pytest_cache/
      - dist/, build/, *.egg-info/

  ELSE IF language == "JavaScript" THEN:
    environment_file: package.json
    environment_content:
      - depends on framework (react, vue, express, etc.)
      - test runner (jest, vitest, mocha)
      - linter (eslint), formatter (prettier)
    
    makefile_targets:
      - install: "npm install"
      - dev: "npm run dev"
      - build: "npm run build"
      - test: "npm test"
      - lint: "npm run lint"
    
    gitignore_additions:
      - node_modules/, package-lock.json (optional)
      - dist/, build/, .next/, .nuxt/
      - .env.local

  IF preserve_legacy == true THEN:
    - Keep existing directory structure (e.g., Analyses/)
    - Create docs/LEGACY_STRUCTURE.md explaining old structure
    - Add cross-references in README pointing to new structure
    - Plan gradual migration in CONTRIBUTING.md
  
  ELSE:
    - Reorganize existing directories
    - More disruptive but cleaner result
    - Requires stakeholder communication

TEMPLATES_TO_APPLY:
  1. .gitignore                          # Language + OS + project_type specific
  2. CITATION.cff or codemeta.json       # If has_paper or needs attribution
  3. environment.yml / package.json      # Language-specific dependencies
  4. Makefile or scripts/setup.sh        # Automation & development workflow
  5. README.md                           # Improved with hierarchy & links
  6. docs/README.md                      # Documentation index
  7. docs/REPRODUCIBILITY.md or CONTRIBUTING.md
  8. docs/API.md or METHODS.md or ARCHITECTURE.md
  9. .github/ISSUE_TEMPLATE/             # Context-aware issue templates
  10. .github/workflows/ (optional)      # CI/CD templates

OUTPUT:
  Restructuring Plan:
    directories_to_create: []
    files_to_generate: []
    dependencies: {}
    preserve_items: []
    cleanup_items: []
    estimated_files: integer
  
  File Manifest:
    - filename: path
      type: config | documentation | template
      context_variables: {}
      dependencies: []
```

---

## Phase 3: Content Generation

### Generation Order (Atomic Commits)

```yaml
GENERATION_ORDER:

  BATCH_1: Core Configuration
    1. .gitignore
       - Inputs: language, os, project_type
       - Content: Language-specific patterns + OS artifacts + project outputs
       - Commit: "Add comprehensive gitignore for {language} {project_type}"
    
    2. CITATION.cff (if has_paper OR has_authors)
       - Inputs: authors, doi, paper_title, paper_year, paper_journal, repository_url
       - Content: CFF-format citation metadata
       - Commit: "Add citation metadata in CFF format"
    
    3. LICENSE (if not exists)
       - Inputs: recommended MIT for research, check existing
       - Content: MIT License text (or preserve existing)
       - Commit: "Add MIT License"

  BATCH_2: Environment & Automation
    4. environment.yml / package.json / pyproject.toml
       - Inputs: language, dependencies_for_type
       - Content: Package manager configuration
       - Commit: "Add {language} dependency specification"
    
    5. Makefile or scripts/
       - Inputs: language, project_type, targets_needed
       - Content: Common development/research automation targets
       - Commit: "Add Makefile with {type}-specific targets"

  BATCH_3: Primary Documentation
    6. README.md (update/improve existing)
       - Inputs: project_title, description, key_features, documentation_index
       - Content: Hierarchical, badge-rich, link to docs/
       - Commit: "Update README with reproducibility-first approach"
    
    7. docs/README.md (new)
       - Inputs: documentation_files_list, quick_nav
       - Content: Documentation index and navigation
       - Commit: "Add documentation index and navigation guide"

  BATCH_4: Specialized Documentation (3-4 files)
    For research projects:
      8. docs/REPRODUCIBILITY.md
         - Inputs: language, dependencies, steps
         - Content: Step-by-step reproduction guide
         - Commit: "Add comprehensive reproducibility guide"
      
      9. docs/DATA_SOURCES.md
         - Inputs: has_data, data_structure
         - Content: Data descriptions, sources, processing
         - Commit: "Add data sources and description documentation"
      
      10. docs/METHODS.md
          - Inputs: methodology, tools, parameters
          - Content: Technical methodology details
          - Commit: "Add comprehensive methodology documentation"
    
    For library projects:
      8. docs/INSTALLATION.md
      9. docs/API.md
      10. docs/CONTRIBUTING.md
    
    For webapp projects:
      8. docs/ARCHITECTURE.md
      9. docs/DEPLOYMENT.md
      10. docs/API.md

  BATCH_5: Directory Structure
    11. Create directory placeholders
        - data/.gitkeep (if has_data)
        - analysis/scripts/.gitkeep
        - analysis/functions/.gitkeep
        - analysis/notebooks/.gitkeep
        - results/figures/.gitkeep
        - results/tables/.gitkeep
        - results/logs/.gitkeep
        - src/.gitkeep (if library)
        - tests/.gitkeep (if library)
        - docs/.gitkeep
        - Commit: "Add directory structure with .gitkeep placeholders"

  BATCH_6: GitHub Templates
    12. .github/ISSUE_TEMPLATE/reproducibility_report.md (research)
        - Content: Form for reporting reproduction issues
        - Commit: "Add GitHub issue template for reproducibility reports"
    
    13. .github/ISSUE_TEMPLATE/bug_report.md (library)
        - Content: Standard bug report template
        - Commit: "Add GitHub issue template for bug reports"
    
    14. .github/pull_request_template.md (if library/webapp)
        - Content: PR submission guidelines
        - Commit: "Add GitHub pull request template"

  BATCH_7: CI/CD Templates (optional, advanced)
    15. .github/workflows/ci.yml (if library/webapp)
        - Content: GitHub Actions for testing & linting
        - Commit: "Add GitHub Actions CI/CD workflow"

COMMIT_MESSAGE_PATTERN:
  "{action}: {description} for {context}"
  
  Examples:
    "Add comprehensive gitignore for reproducible research project"
    "Add Conda environment specification with R and analysis packages"
    "Add Makefile with analysis automation targets"
    "Update README with reproducibility-first approach"
    "Add documentation index and navigation guide"
    "Add comprehensive reproducibility guide"
    "Add comprehensive data sources and description documentation"
    "Add comprehensive methodology documentation"
    "Add directory structure with .gitkeep placeholders"
    "Add GitHub issue template for reproducibility reports"

TOTAL_FILES_CREATED:
  Minimal setup: 7-9 files
  Standard setup: 12-15 files
  Complete setup: 18-22 files
```

### Content Templates (Context-Injected)

#### **.gitignore Template**
```yaml
TEMPLATE_VARIABLES:
  - language: string (R, Python, JavaScript, etc.)
  - os: string (macOS, Linux, Windows)
  - project_type: string (research, library, webapp)
  - has_data: boolean
  - has_notebooks: boolean

CONTENT_STRUCTURE:
  # Header with metadata
  # OS artifacts (always)
  # Language-specific patterns
  # IDE/Editor patterns
  # Project-specific patterns (based on project_type)
  # Data patterns (if has_data)
  # Build/dist patterns (if library/webapp)

EXAMPLE_SECTIONS:
  # R-specific
  .Rhistory
  .Rapp.history
  .RData
  .Ruserdata
  .Rproj.user
  *.Rproj
  
  # Python-specific
  __pycache__/
  *.pyc, *.pyo
  venv/, env/, .venv/
  .pytest_cache/
  
  # Research-specific (data, results)
  data/raw/*
  !data/raw/.gitkeep
  results/
  *.log
```

#### **CITATION.cff Template**
```yaml
TEMPLATE_VARIABLES:
  - has_paper: boolean
  - authors: array[string]
  - doi: string (optional)
  - paper_title: string (optional)
  - paper_year: integer (optional)
  - paper_journal: string (optional)
  - repository_url: string

CONTENT_STRUCTURE (if has_paper):
  cff-version: 1.2.0
  message: "If you use this research or code, please cite the published paper..."
  authors: [list]
  title: {paper_title}
  journal: {paper_journal}
  year: {paper_year}
  volume: (if available)
  article: (if available)
  doi: {doi}
  url: {repository_url}

CONTENT_STRUCTURE (if library/no paper):
  cff-version: 1.2.0
  message: "Please cite this software..."
  authors: [list]
  title: {project_title}
  repository-code: {repository_url}
  repository: {repository_url}
  type: software
  version: (from package.json or setup.py)
  license: MIT (or detected)
```

#### **environment.yml Template (R)**
```yaml
TEMPLATE_VARIABLES:
  - r_version: "4.3.0"
  - has_tidyverse: boolean (usually true)
  - has_bioconductor: boolean (if research)
  - has_rmarkdown: boolean (if notebooks)
  - has_shiny: boolean (if interactive)
  - additional_packages: array[string]

CONTENT:
  name: {project_name_slug}-env
  channels:
    - conda-forge
    - bioconda
    - defaults
  dependencies:
    - r-base={r_version}
    - r-essentials
    - r-tidyverse (if has_tidyverse)
    - r-ggplot2, r-dplyr, r-tidyr, r-readr, r-readxl
    - r-rmarkdown, r-knitr (if has_rmarkdown)
    - bioconductor-deseq2 (if has_bioconductor)
    - r-shiny (if has_shiny)
    - [additional_packages]
  
  pip:
    - numpy, pandas, scipy, scikit-learn
    - matplotlib, seaborn
```

#### **Makefile Template (Research + R)**
```makefile
.PHONY: help validate install-env run-analyses clean

help:
	@echo "Available targets:"
	@echo "  make install-env       - Create conda environment"
	@echo "  make validate          - Validate repository structure"
	@echo "  make run-analyses      - Execute all R scripts"
	@echo "  make clean             - Remove R session artifacts"

install-env:
	conda env create -f environment.yml
	@echo "Environment created. Activate: conda activate {env_name}"

validate:
	@test -d Analyses/codes || (echo "ERROR: Analyses/codes not found"; exit 1)
	@test -d data || (echo "ERROR: data directory not found"; exit 1)
	@test -f README.md || (echo "ERROR: README.md not found"; exit 1)
	@echo "✓ Repository structure validation passed."

run-analyses:
	@for script in Analyses/codes/*.R; do \
		if [ -f "$$script" ]; then \
			echo "Running: $$script"; \
			Rscript "$$script" || exit 1; \
		fi; \
	done
	@echo "✓ All analyses completed."

clean:
	@rm -f .Rhistory .RData Rplots.pdf
	@find . -name "*.Rout" -delete
	@echo "✓ Cleanup complete."
```

#### **README.md Template**
```markdown
STRUCTURE:
  1. Title + Badges (DOI, License, Language)
  2. Quick Overview (1-2 sentences)
  3. Quick Start (3 steps: clone, setup, run)
  4. Repository Structure (tree view)
  5. Documentation Links (→ docs/)
  6. Key Findings or Features
  7. Citation (BibTeX + CITATION.cff)
  8. Make Targets / Commands
  9. License
  10. Contact / Issues

TEMPLATE_VARIABLES:
  - project_title: string
  - description: string
  - has_paper: boolean
  - doi: string (optional)
  - language: string
  - quick_start_steps: array[string]
  - key_features: array[string]
```

#### **docs/REPRODUCIBILITY.md Template (Research)**
```markdown
SECTIONS:
  1. Table of Contents
  2. System Requirements (min & recommended)
  3. Environment Setup (step-by-step)
  4. Running Analyses (3 options)
  5. Expected Outputs
  6. Troubleshooting (common issues + solutions)
  7. Verification Checklist
  8. Performance Notes
  9. Session Information

TEMPLATE_VARIABLES:
  - language: string (R, Python, etc.)
  - environment_manager: string (conda, pip, renv)
  - estimated_runtime: string ("30-60 minutes")
  - memory_requirements: string ("8 GB minimum")
  - common_issues: dict
```

---

## Phase 4: Validation & Cleanup

### Validation Checklist

```yaml
VALIDATION_TASKS:

  File Existence:
    - ✓ All planned files exist in repo
    - ✓ No duplicate files
    - ✓ File sizes reasonable (no huge commits)

  Configuration Syntax:
    - ✓ .gitignore valid syntax
    - ✓ CITATION.cff valid YAML
    - ✓ environment.yml valid YAML
    - ✓ Makefile valid syntax
    - ✓ package.json valid JSON (if JavaScript)
    - ✓ setup.py valid Python (if Python)

  Content Quality:
    - ✓ All Markdown files render correctly
    - ✓ No typos in critical paths
    - ✓ Internal links (to other docs) are valid
    - ✓ Code blocks are properly formatted
    - ✓ Commit messages are descriptive

  Repository Cleanliness:
    - ✓ No .DS_Store files tracked
    - ✓ No __pycache__ directories tracked
    - ✓ No node_modules tracked
    - ✓ No venv/ directories tracked
    - ✓ No generated files tracked (dist/, build/, etc.)

  Structure Validation:
    - ✓ Directory hierarchy is logical
    - ✓ README links match actual docs/
    - ✓ All .gitkeep placeholders created
    - ✓ No empty directories without .gitkeep

  Context-Specific Validations:
    
    IF language == "R":
      - ✓ environment.yml contains r-base
      - ✓ Makefile has Rscript targets
      - ✓ .gitignore excludes .Rhistory, .RData
      - ✓ docs/REPRODUCIBILITY.md mentions conda activate
    
    IF language == "Python":
      - ✓ requirements.txt or pyproject.toml exists
      - ✓ Makefile has pip install target
      - ✓ .gitignore excludes __pycache__, venv/
      - ✓ setup.py present (if library)
    
    IF project_type == "research":
      - ✓ data/ directory exists
      - ✓ docs/DATA_SOURCES.md exists
      - ✓ docs/REPRODUCIBILITY.md exists
      - ✓ docs/METHODS.md exists
      - ✓ Makefile has validate target
    
    IF project_type == "library":
      - ✓ src/ or lib/ directory exists
      - ✓ tests/ directory exists
      - ✓ docs/API.md or docs/CONTRIBUTING.md exists
      - ✓ Package manifest (setup.py, package.json, etc.) present
    
    IF has_paper:
      - ✓ CITATION.cff contains valid DOI
      - ✓ README mentions paper with link
      - ✓ All authors listed in CITATION.cff

OUTPUT:
  Validation Report:
    passed_checks: []
    failed_checks: []
    warnings: []
    status: "PASS" | "PASS_WITH_WARNINGS" | "FAIL"
    
  Recommended Actions:
    cleanup_items: []
    fix_items: []
    optional_improvements: []
```

### Cleanup Operations

```yaml
CLEANUP_OPERATIONS:

  1. Remove Tracked OS Artifacts
     git rm -r --cached .DS_Store 2>/dev/null || true
     git rm -r --cached .AppleDouble 2>/dev/null || true
     git rm -r --cached Thumbs.db 2>/dev/null || true

  2. Remove Tracked Language Artifacts
     git rm -r --cached __pycache__ 2>/dev/null || true
     git rm -r --cached node_modules 2>/dev/null || true
     git rm -r --cached *.pyc 2>/dev/null || true

  3. Verify .gitignore Updated
     Confirm .gitignore includes all artifacts
     Run: git status (should show no unwanted files)

  4. Cleanup Commit
     git add .gitignore
     git commit -m "Clean up OS artifacts and update gitignore"

OUTPUT:
  Cleanup Report:
    files_removed: []
    artifacts_cleaned: []
    status: "CLEAN" | "WARNINGS"
    next_steps: []
```

---

## Phase 5: Delivery & Documentation

### Delivery Summary

```yaml
DELIVERY_ARTIFACTS:

  1. Restructured Repository
     - Location: main branch
     - All 15-22 files committed
     - Clean git history (no OS artifacts)
     - Ready for immediate use

  2. Change Summary Document
     Content:
       - What Changed (directories, files, purpose)
       - Why (reproducibility, best practices, scalability)
       - Impact (what users need to know)
       - Migration Path (if applicable)

  3. User Guide
     Sections:
       - Quick Start (3 steps)
       - Directory Guide (what each folder is for)
       - Documentation Map (which file for what question)
       - Recommended Next Steps (1-5 actionable items)
       - FAQ (5-10 common questions)

  4. Technical Reference
     Includes:
       - All make targets and their purpose
       - Environment variable setup
       - Key file descriptions
       - Troubleshooting matrix

OUTPUT_FORMAT:

RESTRUCTURING_COMPLETE ✓

📊 SUMMARY
─────────────────────────────────────
Project Type:        {project_type}
Primary Language:    {language}
Files Created:       {count}
Commits Made:        {count}
Structure Score:     {0-100} → {improved 0-100}

📁 NEW STRUCTURE
─────────────────────────────────────
{tree view of new structure}

🔧 KEY COMPONENTS ADDED
─────────────────────────────────────
✓ Configuration (.gitignore, CITATION.cff, environment.yml)
✓ Automation (Makefile with X targets)
✓ Documentation (4 files in docs/)
✓ Directory Structure (9 directories with .gitkeep)
✓ GitHub Templates (Issue & PR templates)

🚀 QUICK START
─────────────────────────────────────
1. make install-env
2. conda activate {env_name}
3. make validate
4. Read: docs/README.md

📚 DOCUMENTATION MAP
─────────────────────────────────────
→ How to reproduce?        docs/REPRODUCIBILITY.md
→ Data information?        docs/DATA_SOURCES.md
→ Technical methods?       docs/METHODS.md
→ General guidance?        docs/README.md

⚠️ IMPORTANT NOTES
─────────────────────────────────────
{if preserve_legacy}
- Legacy structure (Analyses/) preserved
- New structure ready: analysis/, data/, results/
- Gradual migration recommended

{if cleanup_needed}
- OS artifacts cleaned from git history
- .gitignore updated for future prevention

📝 RECOMMENDED NEXT STEPS
─────────────────────────────────────
1. Review docs/README.md for documentation index
2. Run: make validate && make install-env
3. Test: conda activate {env_name} && make validate
4. If legacy: plan migration to new structure
5. Update team about new structure

❓ FAQ
─────────────────────────────────────
Q: Where do I put raw data?
A: data/raw/ — don't modify these files
Q: How do I run analyses?
A: conda activate {env_name} && make run-analyses
Q: Can I still use the old Analyses/ directory?
A: Yes, but new structure (analysis/) is recommended
Q: How do I reproduce everything?
A: Follow docs/REPRODUCIBILITY.md step-by-step
```

---

## Conditional Logic & Templates

### Language-Specific Rules

```yaml
LANGUAGE_RULES:

  R:
    environment_manager: conda + renv (optional)
    environment_file: environment.yml
    r_base_default: "4.3.0"
    r_packages: [tidyverse, ggplot2, dplyr, tidyr, readr, readxl, rmarkdown, knitr]
    r_bioconductor: [deseq2, limma, edger] (if research)
    makefile_targets: [install-env, validate, run-analyses, clean]
    validation_checks:
      - .Rhistory, .RData in .gitignore
      - *.Rproj in .gitignore
      - environment.yml valid YAML
      - Makefile has Rscript target
    
  Python:
    environment_manager: pip | conda | poetry
    environment_file: requirements.txt | pyproject.toml | Pipfile
    python_version_default: "3.10"
    python_packages: [numpy, pandas, scipy, scikit-learn]
    test_framework: pytest (default)
    makefile_targets: [install-env, test, lint, run-analyses, clean]
    validation_checks:
      - __pycache__ in .gitignore
      - venv/ in .gitignore
      - *.pyc in .gitignore
      - setup.py (if library)
      - pytest.ini (if library)
  
  JavaScript/TypeScript:
    environment_manager: npm | yarn | pnpm
    environment_file: package.json
    node_version_default: "18.x"
    test_framework: jest | vitest | mocha
    linter: eslint
    formatter: prettier
    makefile_targets: [install, dev, build, test, lint]
    validation_checks:
      - node_modules in .gitignore
      - dist/ or build/ in .gitignore
      - .env.local in .gitignore
      - package-lock.json (or equivalent)

PROJECT_TYPE_RULES:

  research:
    mandatory_dirs: [data/raw, data/processed, analysis, results, docs]
    mandatory_docs: [REPRODUCIBILITY.md, DATA_SOURCES.md, METHODS.md]
    mandatory_files: [.gitignore, Makefile, environment.yml]
    emphasis: reproducibility, data management, methodology
    targets: validate, run-analyses
    
  library:
    mandatory_dirs: [src, tests, docs, examples]
    mandatory_docs: [API.md, CONTRIBUTING.md, INSTALLATION.md]
    mandatory_files: [package.json|setup.py, pytest.ini|jest.config, .github/workflows]
    emphasis: code quality, testing, API clarity
    targets: test, lint, build, publish
    
  webapp:
    mandatory_dirs: [frontend, backend, infrastructure, docs]
    mandatory_docs: [ARCHITECTURE.md, DEPLOYMENT.md, API.md]
    mandatory_files: [docker-compose.yml, .env.example, .github/workflows]
    emphasis: deployment automation, environment management
    targets: dev, build, deploy
    
  cli:
    mandatory_dirs: [src, tests, docs]
    mandatory_docs: [INSTALLATION.md, USAGE.md, API.md]
    mandatory_files: [package.json|setup.py|Cargo.toml]
    emphasis: usability, documentation, quick start
    targets: build, install, test

PRESERVE_LEGACY_RULES:

  IF preserve_legacy == true:
    - Keep existing: Analyses/, legacy/, old_code/
    - Create: docs/LEGACY_STRUCTURE.md
    - Add README section: "Legacy Structure & Migration"
    - Document: mapping between old and new
    - Plan: gradual migration strategy
    - Timeline: when to deprecate old structure

SPECIAL_CASES:

  IF has_paper == true:
    - Mandatory: CITATION.cff with DOI
    - Mandatory: Paper reference in README
    - Recommended: paper_url, paper_journal
    - Include: All authors in CITATION.cff
    - Add: @article BibTeX in README
  
  IF has_data == true:
    - Mandatory: docs/DATA_SOURCES.md
    - Mandatory: data/raw/, data/processed/
    - Recommended: Data dictionary (in DATA_SOURCES.md)
    - Recommended: Preprocessing scripts (in analysis/scripts/)
    - Validate: No large files in git (use .gitignore)
  
  IF monorepo == true:
    - Create: Root README with package overview
    - Create: packages/{package_name}/README.md per package
    - Create: Root Makefile with targets for all packages
    - Create: docs/ARCHITECTURE.md explaining package relations
```

---

## Edge Cases & Troubleshooting

```yaml
EDGE_CASES:

  1. Repository Already Partially Structured
     Detection: Some dirs/files exist (data/, docs/, tests/)
     Handling:
       - Preserve existing structure
       - Extend with missing pieces
       - Merge documentation
       - Update Makefile with all targets
     Risk: Conflicts in existing README
     Solution: Create README.new, let user merge
  
  2. Very Large Repository (>1GB)
     Detection: During Phase 1 analysis
     Handling:
       - Warn about .gitignore criticality
       - Recommend git-lfs for large files
       - Focus on data/ structure documentation
       - Create docs/DATA_MANAGEMENT.md
     Risk: Slow commits if large files tracked
     Solution: Pre-cleanup before restructuring
  
  3. Multiple Languages (Monorepo)
     Detection: .py, .js, .R files in root
     Handling:
       - Ask: Which is primary?
       - Create environment files for each
       - Create docs/ARCHITECTURE.md
       - Makefile targets for each language
     Risk: Conflicting dependencies
     Solution: Language-specific subdirectories
  
  4. Existing .gitignore
     Detection: .gitignore already present
     Handling:
       - Merge (not overwrite)
       - Show differences
       - Keep user's patterns
       - Add recommended patterns
     Risk: Losing important exclusions
     Solution: Create .gitignore.backup before merge
  
  5. Existing Documentation
     Detection: docs/ folder or README.md exists
     Handling:
       - Preserve existing docs
       - Create docs/OLD_README_BACKUP.md if needed
       - Restructure to match template
       - Cross-link old and new
     Risk: Loss of existing documentation
     Solution: Archive in docs/archive/
  
  6. Repository with Secrets (API keys in history)
     Detection: .env, secrets in old commits
     Handling:
       - Alert user before proceeding
       - Recommend: git-filter-branch or BFG
       - Create .env.example with sanitized keys
       - Add docs/SECRETS_MANAGEMENT.md
     Risk: Exposing sensitive data
     Solution: Require cleanup before restructuring
  
  7. No README Present
     Detection: Missing README.md
     Handling:
       - Create comprehensive README from scratch
       - Infer description from repo name/files
       - Ask user for brief description
       - Create template with TODOs
     Risk: Empty or generic README
     Solution: Interactive prompts for key fields

TROUBLESHOOTING_MATRIX:

  Problem: Files not created
  Causes:
    - Branch protection rules blocking pushes
    - Insufficient permissions
    - Quota exceeded (large file limit)
  Solution:
    - Check branch settings
    - Verify user has write access
    - Verify file size < 100MB each

  Problem: Makefile targets not working
  Causes:
    - Missing environment.yml
    - Conda not installed
    - Wrong environment name
  Solution:
    - Run: make install-env first
    - Check: conda --version
    - Verify: conda activate {env_name}

  Problem: .gitignore not preventing files
  Causes:
    - Files already tracked before .gitignore added
    - Typo in .gitignore pattern
    - Pattern syntax error
  Solution:
    - git rm --cached {file}
    - Verify pattern with online regex tools
    - Run: git status (should show removed)

  Problem: Documentation links broken
  Causes:
    - Typo in file paths
    - Case sensitivity (docs/ vs Docs/)
    - File not created
  Solution:
    - Run validation in Phase 4
    - Check case consistency
    - Verify all referenced files exist

  Problem: Commit history corrupted
  Causes:
    - Network interruption during push
    - Race condition with other users
    - Rebase/force-push conflict
  Solution:
    - Run: git fsck --full
    - Create backup branch: git branch backup
    - Reset if needed: git reset --hard {safe_commit}
```

---

## Implementation Checklist

```yaml
PRE_EXECUTION:
  ☐ User has write access to repository
  ☐ Repository backed up (or on non-critical branch)
  ☐ No active PRs with conflicting changes
  ☐ All input parameters validated
  ☐ Git history checked for secrets (optional)

EXECUTION:
  ☐ Phase 1: Discovery & Analysis complete
  ☐ Phase 2: Planning & Design approved
  ☐ Phase 3: Content Generation complete
  ☐ Phase 4: Validation & Cleanup passed
  ☐ Phase 5: Delivery Summary generated

POST_EXECUTION:
  ☐ User reviewed change summary
  ☐ User confirmed repository structure
  ☐ User read docs/README.md
  ☐ User tested make validate
  ☐ User tested make install-env
  ☐ Team notified of restructuring
  ☐ Documentation updated (if needed)
```

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-05-29 | Initial release with full 5-phase workflow |

---

## References & Related Skills

- **pr-understanding**: Entity model for comprehensive PR analysis
- **stack-trace-debugging**: Root cause analysis for errors
- **pr-summary**: Summarizing changes and impact
- **Repository Automation**: CI/CD, pre-commit hooks, linting

---

## License

This skill documentation is provided as-is for repository restructuring workflows. Apply context-aware judgment to each repository type and always validate output before merging to production.
