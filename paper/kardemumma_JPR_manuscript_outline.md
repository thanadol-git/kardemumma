# KARDEMUMMA — Manuscript Outline for *J. Proteome Res.* Special Issue "Software Tools and Resources 2027"

**Deadline:** 2 November 2026 · **Submit via:** ACS Publishing Center → Special Issue Selection: *Software Tools and Resources 2027*

---

## 0. Key decisions before writing

| Item | Recommendation | Why |
|---|---|---|
| Manuscript type | **Article** (limit 8,000 words) | KARDEMUMMA is a new, not previously published tool. The Special Issue reserves *Technical Notes* (limit 4,500) for substantial updates of already-published tools. |
| Target length | **~4,800 counted words** (well under 8,000) | The call explicitly asks authors to be concise and focus on novel functionality. |
| What counts | Introduction + Results & Discussion + Conclusions | JPR excludes the Experimental Section, Acknowledgments, SI statement and References from the word count. Abstract is capped separately at 200 words. |
| Display items | 5 figures + 2 tables in main text; rest in SI | Table of input → operation → output is **mandatory** for tools/libraries in this call. |
| Positioning | "Library of functional building blocks + CLI for targeted (PRM/SRM) plasma proteomics QC and absolute quantification" | Editors prefer GUI/web tools but accept well-documented libraries/APIs; say this explicitly and emphasise the generated Jupyter reports as the user-facing layer. |

### Word budget at a glance

| Section | Words | Counted? |
|---|---:|:---:|
| Abstract | ≤200 | separate cap |
| Introduction | ~900 | ✔ |
| Experimental Section | ~1,600–2,000 | ✘ |
| Results & Discussion | ~3,650 | ✔ |
| Conclusions | ~250 | ✔ |
| Data & Software Availability, SI list, Author info, Acknowledgments, References (~40–50) | — | ✘ |
| **Counted total** | **~4,800** | |

---

## 1. Title (pick one, ≤ ~20 words)

- *KARDEMUMMA: A Python Toolkit for Reproducible Quality Control and Absolute Quantification of Targeted Plasma Proteomics Assays*
- *KARDEMUMMA: From Skyline Exports to SDRF-Annotated, Batch-Corrected Absolute Protein Concentrations in Targeted Mass Spectrometry*
- *KARDEMUMMA: Standardized QC, Plate-Effect Correction, and qRePS-Based Absolute Quantification for PRM/SRM Plasma Proteomics*

## 2. Abstract (≤200 words)

| Sentence(s) | Content | ~Words |
|---|---|---:|
| 1–2 | Problem: targeted MS (PRM/SRM) is moving toward clinical plasma studies across many plates, but QC, normalization and absolute quantification after Skyline/OpenSWATH are done with ad hoc scripts; metadata rarely standardized. | 40 |
| 3 | We present KARDEMUMMA (version x.y.z), an open-source Python package + CLI. | 20 |
| 4–5 | What it does: imports Skyline/OpenSWATH exports, links them to targeted-SDRF metadata, three-step pipeline (preview → cutoff → report): dot-product & heavy/light filtering, intra-/inter-plate CV, ANOVA-based plate correction with pool normalization, and absolute quantification from qRePS spike-in standards retrieved automatically by lot number. | 60 |
| 6 | Demonstration: applied to *N* plasma samples across *P* plates (datasets …); key numbers (e.g., peptides retained, median inter-plate CV before/after, agreement with vendor software). Biological outcomes with examples | 45 |
| 7 | Availability: PyPI, GitHub (MIT), Zenodo DOI; reports as Jupyter notebooks. | 25 |

## 3. Keywords (≤10)

targeted proteomics; parallel reaction monitoring (PRM); selected reaction monitoring (SRM); plasma proteomics; quality control; absolute quantification; stable isotope-labeled standards; SDRF-Proteomics; batch correction; Python

## 4. TOC / Abstract graphic

Left: Skyline/OpenSWATH export + targeted-SDRF + qRePS lot number → centre: three-step CLI (preview / cutoff / report) → right: QC dashboard thumbnails + absolute concentration (e.g., fmol/µL plasma) output. Include the cardamom logo.

---

## 5. Introduction (~900 words, 4 paragraphs)

| ¶ | Topic | Tentative content | ~Words | Status |
|---|---|---|---:|:---:|
| 1 | Why targeted plasma proteomics | PRM/SRM with stable-isotope-labeled (SIL) standards as the bridge from discovery to clinical assays; precision, multiplexing, absolute quantification; growth of large cohort studies processed in 96-well plates. | 200 | ✅ |
| 2 | The gap after peak integration | Skyline (and OpenSWATH) do integration well, but downstream steps (transition/peptide QC, plate/batch effects, pool-based normalization, conversion to absolute concentrations) are lab-specific scripts → poor reproducibility, hard to audit for clinical translation. Metadata usually lives in spreadsheets not linked to raw files. | 250 | ✅ |
| 3 | Existing tools and their limits | Briefly compare (full comparison in Table 2): Skyline/Panorama & AutoQC (instrument/system suitability QC), MSstats / MSstatsQC (statistical modeling, longitudinal QC), TargetedMSQC, community QC standards (mzQC, pmultiqc), SDRF-Proteomics (mostly DDA/DIA-oriented). None connects SDRF metadata, plate-aware QC, and vendor-lot-driven absolute quantification in one reproducible workflow. **[verify each tool's current scope before citing]** | 250 | ✅ |
| 4 | This work | Introduce KARDEMUMMA: name/acronym, version, language, license; core contributions as 3–4 bullets-in-prose (targeted-SDRF integration, three-step CLI with cutoff guidance, ANOVA-based plate correction, automated qRePS absolute quantification, notebook report generation); developed and used at KTH/SciLifeLab; summary of demonstration datasets. | 200 |   |

---

## 6. Experimental Section (~1,600–2,000 words; *not* counted)

Write in full sentences (JPR forbids bullet format here).

| Subsection | Tentative content | ~Words |
|---|---|---:|
| 6.1 Software architecture and implementation | Python ≥3.10; package layout under `src/kardemumma` (importer, sdrf, prm, prm_plots, proteomedge, openswath modules); dependencies (pandas, statsmodels/patsy, pyteomics, matplotlib/seaborn, sdrf-pipelines parser…); CLI entry points `kardemumma-preview`, `-cutoff`, `-report`, `-ipynb`; design as composable functions re-exported at `kdm.*`. Refer to Figure 1. | 300 |
| 6.2 Input data and metadata model | Required Skyline CSV columns (13 columns incl. `Normalized Area`, `RatioLightToHeavy`, `Library/Ratio Dot Product`); OpenSWATH TSV support; targeted-SDRF: `comment[data file]` join key, `characteristics[...]`, `factor value[...]`, `comment[ProteomEdge]` lot column; file-name normalization and Skyline↔SDRF cross-check; SDRF validation via `parse_sdrf`. | 300 |
| 6.3 Quality-control metrics and filtering | Library dot-product threshold; heavy/light detection count filters; missing-value flags; retention-time deviation (observed vs predicted/iRT); ratio dot product; default thresholds and how `preview` suggests cutoffs. | 300 |
| 6.4 Precision, normalization and plate-effect correction | Intra-plate CV (pool replicates) and inter-plate CV (plate means); lowest-CV peptide selection by percentile; two-way ANOVA `ratio ~ Plate + Peptide` (log scale); plate conversion factors from median-centered ratios; adjusted L/H ratios. Give equations. | 350 |
| 6.5 Absolute quantification with qRePS standards | Reading lot number from SDRF; retrieval of the qRePS target table and FASTA from the public ProteomEdge lot page; peptide-to-protein mapping; conversion of L/H ratio × spiked amount → concentration; units. Note network dependency and caching/offline option. | 250 |
| 6.6 Datasets used for demonstration | For each dataset (e.g., DA4K, MORPHEUS, PD — **confirm names/allowed disclosure**): cohort, sample numbers, plates, pools, PRM/SRM instrument & method, Skyline version, spike-in lot. **Ethics approval number + informed consent statement (mandatory for human plasma).** ProteomeXchange/PanoramaPublic identifiers + reviewer credentials. | 350 |
| 6.7 Benchmarking | Benchmark notebook: comparison with manual/legacy workflow or vendor results (concordance of concentrations, CV improvement, runtime vs number of samples/peptides, memory). | 150 |
| 6.8 Software development practices and use of AI tools | pytest unit tests, GitHub Actions CI, versioned releases, PyPI publishing, Zenodo archiving, API docs (HTML/PDF) per release, issue tracking. **Explicitly state whether/how LLMs were used for code generation and how that code was validated** (the call requires this). | 150 |

---

## 7. Results and Discussion (~3,650 words)

| Subsection | Tentative content | Display item | ~Words |
|---|---|---|---:|
| 7.1 Overview of KARDEMUMMA | Problem → solution walkthrough; three user layers (Python API, CLI, generated notebooks); how a user goes from exports to report in three commands. | **Figure 1**: architecture/workflow diagram. **Table 1** (required): each module/command → input → operations → output. | 500 |
| 7.2 Targeted-SDRF as the backbone of reproducibility | Why SDRF for targeted experiments (link to proposed targeted-SDRF template); automated validation, whitespace cleaning, Skyline↔SDRF file cross-check catching mislabelled/missing runs; plate derivation; qRePS lot captured as metadata so absolute quantification is traceable. Show an example of an error caught. | Figure S1 (SDRF example + validation report) | 450 |
| 7.3 Step 1 — `preview`: data-driven QC and cutoff guidance | Dot-product distributions, heavy/light count scatter, RT deviation, missingness; how suggested cutoffs are derived; QC-sample detection. | **Figure 2**: multi-panel QC (dot-product KDE, heavy vs light density scatter, RT deviation, missingness heatmap). | 600 |
| 7.4 Step 2 — `cutoff`: filtering, pool normalization and plate-effect correction | Peptides/proteins retained at each filter (funnel); intra- vs inter-plate CV before and after correction; ANOVA variance explained by plate; cumulative peptides vs CV threshold. Quantify improvement (e.g., median inter-plate CV from X% to Y%). | **Figure 3**: pool boxplots by plate before/after, CV KDEs, cumulative peptide-by-CV curve. | 700 |
| 7.5 Step 3 — `report`: absolute quantification | Automatic qRePS retrieval by lot; concentrations per protein/peptide; agreement between peptides of the same protein; comparison with reference/vendor values or literature plasma concentrations (e.g., HPPP/published ranges). | **Figure 4**: concentration dynamic range plot + peptide-concordance scatter. | 500 |
| 7.6 Application to cohort data | Short biological use case (e.g., group comparison via `kardemumma-ipynb --group-a/--group-b`) showing the tool delivers interpretable results; keep biology brief — the tool is the focus. Mention OpenSWATH/other-export input to show generality. | **Figure 5**: group comparison / per-peptide concentration by group. | 500 |
| 7.7 Comparison with existing tools and performance | Feature matrix vs other tools; runtime and scalability from benchmark notebook; reproducibility (same input → identical output across versions). | **Table 2**: feature comparison. Figure S2: runtime scaling. | 450 |
| 7.8 Biological results: demonstration of discovery potential | Example case: application of KARDEMUMMA on a multi-plate plasma cohort revealed differential expression of clinically relevant proteins. Summarize main biological findings (e.g., identification of established and novel plasma biomarkers, concordance with previous literature, power to distinguish group phenotypes via volcano plot and dimensionality reduction). Describe improved quantification precision and batch correction enabling detection of subtle group differences. Note reproducibility and traceability of quantitative results via metadata linkage. | Figure 6 (e.g., volcano plot, t-SNE/UMAP sample separation, boxplots of key proteins) | 450 |
| 7.9 Limitations and roadmap | Currently tied to Skyline column schema and qRePS vendor standards; web scraping dependency; SRM module in development; LOD/LOQ not yet estimated; no GUI yet. Roadmap: mzQC export, pmultiqc integration, SRM, Nextflow pipeline/nf-core. Community contributions welcome. | — | 450 |

## 8. Conclusions (~250 words)

Interpretation, not a recap: what KARDEMUMMA changes for labs running targeted plasma assays (auditable, metadata-linked, plate-aware QC → comparable absolute concentrations across studies/sites), and how it supports translation from research to clinical settings.

---

## 9. End-matter (not counted)

- **Data and Software Availability** — software name, version described, GitHub URL, PyPI, Zenodo DOI (10.5281/zenodo.22882804 → mint a DOI for the *exact* version in the paper), license (MIT), documentation link, test data location; ProteomeXchange/PanoramaPublic IDs + reviewer login.
- **Supporting Information** (non-sentence descriptions with file types), e.g.: "Targeted SDRF template and validation report (PDF)"; "Full parameter reference for CLI commands (PDF)"; "Example Jupyter report generated by kardemumma-ipynb (HTML)"; "Benchmark results and runtime scaling (PDF)"; "Peptide/protein-level results tables (XLSX)".
- **Author Information** — corresponding author, ORCIDs linked in the system (not typed in text).
- **Author Contributions** — CRediT statement (six key developers listed in README).
- **Notes** — competing interests (declare any ProteomEdge relationship).
- **Acknowledgments** — funding (KTH, SciLifeLab, grants), **AI-tool disclosure**.
- **References** — ~40–50, ACS style with full titles.

---

## 10. Pre-submission checklist (gaps found in the repo)

**Required or strongly expected by the call**
- [ ] Tag a release matching the version in the paper and archive it on Zenodo (version-specific DOI). (Thanadol)
- [x] Add a `CITATION.cff` (or `codemeta.json`) to the repo.
- [ ] Deposit raw files + Skyline documents in ProteomeXchange (PRIDE or PanoramaPublic) — required by JPR; author-hosted links are not accepted. (Sathya)
- [ ] Provide small public **test data** + an end-to-end "quick start" that reviewers can run in minutes. (Thanadol)
- [ ] Write the LLM code-generation/validation statement. (Thanadol)
- [ ] Ethics approval number and consent statement for plasma samples. (Fredrik)

**Repository clean-up reviewers will notice**
- [ ] README installation step references `environment.yml`, but the repo contains `config.yml` — fix.
- [ ] "Available pipelines" checkboxes are all unchecked — update to reflect what works.
- [ ] Rename legacy `skyline_qc` / `openms_qc` references in `function_summary.md` to `kardemumma`.
- [ ] Remove or move project-specific notebooks (DA4K, MORPHEUS, PD) and internal sample names from the root README if not intended for public release.
- [ ] Resolve open "Issues" in README (iRT Biognosys sequences, oxidation handling) or document them as known limitations.
- [ ] Complete Phase 1 item 5–6 (consistent naming, type hints, docstrings, stable public API) so Table 1 matches the actual API.
- [ ] Consider a lightweight GUI (e.g., Streamlit or a hosted Voilà/JupyterLite report) — editors state a preference for graphical/web interfaces.

**Cover letter points**
- State it is for the *Software Tools and Resources* Special Issue and why it fits (novel library + CLI for targeted proteomics QC/absolute quant; free, open-source, documented, tested).
- Confirm reviewers can install from PyPI and run the test dataset without charge.
- Note any preprint (bioRxiv/ChemRxiv).

---

### Suggested writing timeline (≈5 weeks)

| Week | Focus |
|---|---|
| 28 Sep – 4 Oct | Repo clean-up, test dataset, CITATION.cff, start PRIDE/PanoramaPublic deposit |
| 5 – 11 Oct | Generate final figures 1–5, Tables 1–2; run benchmark |
| 12 – 18 Oct | Draft Experimental + Results & Discussion |
| 19 – 25 Oct | Introduction, Conclusions, Abstract, SI; co-author review |
| 26 Oct – 1 Nov | Tag release + Zenodo DOI, final checks, cover letter, submit |
