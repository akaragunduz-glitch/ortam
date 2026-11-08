# Reference documents

Every equation in ORTAM is implemented from a document in this folder.
Code comments cite the file's short name, page and equation number, e.g.
`# RAGS Part B (1991), p. 23, Eq. 4`.

All documents are US EPA publications (public domain).
Downloaded on 2026-11-07.

## Risk Assessment Guidance for Superfund (RAGS), Volume I

| File | Short name | Document | Pages | Source |
|------|------------|----------|-------|--------|
| `rags_a_1989.pdf` | RAGS Part A (1989) | Human Health Evaluation Manual, Part A. Interim Final, EPA/540/1-89/002, December 1989 | 291 | https://www.epa.gov/system/files/documents/2024-10/rags_a_508.pdf |
| `rags_b_1991.pdf` | RAGS Part B (1991) | Part B, Development of Risk-based Preliminary Remediation Goals. Interim, Publication 9285.7-01B, December 1991 | 64 | https://semspub.epa.gov/src/document/HQ/156739 |
| `rags_e_2004.pdf` | RAGS Part E (2004) | Part E, Supplemental Guidance for Dermal Risk Assessment. Final, EPA/540/R/99/005, OSWER 9285.7-02EP, July 2004 | 156 | https://19january2025snapshot.epa.gov/sites/default/files/2015-09/documents/part_e_final_revision_10-03-07.pdf |
| `rags_e_2004_impl_memo.pdf` | RAGS Part E memo (2004) | OSWER implementation memorandum for Part E, 16 August 2004 | 4 | https://19january2025snapshot.epa.gov/sites/default/files/2015-09/documents/part_e_impl_2004_final_supp.pdf |
| `rags_f_2009.pdf` | RAGS Part F (2009) | Part F, Supplemental Guidance for Inhalation Risk Assessment. Final, EPA-540-R-070-002, OSWER 9285.7-82, January 2009 | 68 | https://semspub.epa.gov/src/document/HQ/140530 |
| `rags_f_2009_transmittal_memo.pdf` | RAGS Part F memo (2009) | OSWER transmittal memorandum for Part F, 28 January 2009 | 3 | https://semspub.epa.gov/src/document/HQ/140529 |

## Regional Screening Levels (RSLs)

RSL tables version: **November 2024** (worksheet names end in `1124`).

| File | Short name | Content | Source |
|------|------------|---------|--------|
| `rsl_users_guide_2026-03.pdf` | RSL User's Guide | Web page saved as PDF; page last updated 13 March 2026 (268 pages) | https://www.epa.gov/risk/regional-screening-levels-rsls-users-guide |
| `rsl_default_exposure_parameters_2024-11.xlsx` | RSL Default Exposure Parameters | Exposure parameters, PEF, VF, soil to groundwater (12 worksheets) | https://semspub.epa.gov/src/document/HQ/404323 |
| `rsl_summary_table_2024-11.xlsx` | RSL Summary Table (Nov 2024) | Screening levels and toxicity values (TR = 1E-06, THQ = 1) | https://www.epa.gov/risk/regional-screening-levels-rsls-generic-tables |
| `rsl_chemical_parameters_2024-11.xlsx` | RSL Chemical Parameters (Nov 2024) | Chemical-specific parameters | https://www.epa.gov/risk/regional-screening-levels-rsls-generic-tables |
| `rsl_subchronic_toxicity_2024-11.xlsx` | RSL Subchronic Toxicity (Nov 2024) | Subchronic RfD and RfC values | https://www.epa.gov/risk/regional-screening-levels-rsls-generic-tables |

## Other guidance

| File | Short name | Document | Pages | Source |
|------|------------|----------|-------|--------|
| `rais_chemical_risk_users_guide_2026-11.pdf` | RAIS Chemical Risk Guide | RAIS Chemical Risk Calculator User's Guide (ORNL); web page saved as PDF on 2026-11-06; farmer food chain in §4.8 | 159 | https://rais.ornl.gov/tools/rais_chemical_risk_guide.html |
| `epa_exposure_assessment_guidelines_2019.pdf` | EPA Exposure Guidelines (2019) | Guidelines for Human Exposure Assessment, EPA/100/B-19/001, October 2019 | 223 | https://www.epa.gov/risk |

## Soil screening

| File | Short name | Document | Pages | Source |
|------|------------|----------|-------|--------|
| `ssg_supplemental_2002.pdf` | SSG 2002 | Supplemental Guidance for Developing Soil Screening Levels for Superfund Sites, OSWER 9355.4-24, December 2002 (main text and references; the appendices are separate documents) | 106 | https://semspub.epa.gov/src/document/hq/175878 |

## Vapor intrusion

| File | Short name | Document | Pages | Source |
|------|------------|----------|-------|--------|
| `je_documentation_v6_2017.pdf` | J&E v6 (2017) | Documentation for EPA's Implementation of the Johnson and Ettinger Model to Evaluate Site Specific Vapor Intrusion into Buildings, Version 6.0 | 82 | https://semspub.epa.gov/work/HQ/100000489.pdf |
| `visl_users_guide_2014.pdf` | VISL UG (2014) | Vapor Intrusion Screening Level (VISL) Calculator User's Guide, May 2014 | 11 | https://semspub.epa.gov/work/HQ/196703.pdf |
| `oswer_vi_technical_guide_2015.pdf` | OSWER VI Guide (2015) | OSWER Technical Guide for Assessing and Mitigating the Vapor Intrusion Pathway from Subsurface Vapor Sources to Indoor Air, OSWER 9200.2-154, June 2015 | 267 | https://www.epa.gov/sites/default/files/2015-09/documents/oswer-vapor-intrusion-technical-guide-final.pdf |
| `je_spreadsheet_v6_2017.xlsm` | J&E spreadsheet v6 | EPA J&E model spreadsheet, Version 6.0, September 2017 (downloaded as 100000499.xlsm); data tables imported to `data/reference/je/v6-2017/` | - | https://semspub.epa.gov/work/HQ/100000499.xlsm |
| `visl_calculator_resident_2026-11-08.xlsx` | VISL output (2026-11-08) | Output of the EPA VISL calculator, resident, THQ 0.1, benzene/PCE/TCE/vinyl chloride; reference values for `vapor.visl` | - | https://epa-visl.ornl.gov/ |

The downloaded SSG and J&E documentation files had damaged cross-reference
tables (text could be read from page 1 only); they were repaired with `qpdf`
(content unchanged).

## Food chain (HHRAP 2005)

| File | Short name | Document | Pages | Source |
|------|------------|----------|-------|--------|
| `hhrap_2005_ch5.pdf` | HHRAP ch. 5 | Human Health Risk Assessment Protocol for Hazardous Waste Combustion Facilities, Final, EPA530-R-05-006, September 2005: Chapter 5, Estimation of Media Concentrations | 90 | https://archive.epa.gov/epawaste/hazard/tsd/td/web/pdf/05hhrap5.pdf |
| `hhrap_2005_ch6.pdf` | HHRAP ch. 6 | Chapter 6, Quantifying Exposure | 24 | https://archive.epa.gov/epawaste/hazard/tsd/td/web/pdf/05hhrap6.pdf |
| `hhrap_2005_ch7.pdf` | HHRAP ch. 7 | Chapter 7, Risk and Hazard Characterization | 17 | https://archive.epa.gov/epawaste/hazard/tsd/td/web/pdf/05hhrap7.pdf |
| `hhrap_2005_appendix_a.pdf` | HHRAP App. A | Appendix A, Chemical-Specific Parameter Values | 81 | https://archive.epa.gov/epawaste/hazard/tsd/td/web/pdf/05hhrapapa.pdf |
| `hhrap_2005_appendix_b.pdf` | HHRAP App. B | Appendix B, Equations and Default Values for Estimating Media Concentrations | 285 | https://archive.epa.gov/epawaste/hazard/tsd/td/web/pdf/05hhrapapb.pdf |
| `hhrap_2005_appendix_c.pdf` | HHRAP App. C | Appendix C, Equations and Default Values for Cancer Risk and Non-Cancer Effects | 53 | https://archive.epa.gov/epawaste/hazard/tsd/td/web/pdf/05hhrapapc.pdf |
| `hhrap_2005_chemical_database.mdb` | HHRAP database | HHRAP Companion Database (chemical-specific parameters, MS Access) | - | https://archive.epa.gov/epawaste/hazard/tsd/td/web/mdb/05hhrapchemdat.mdb |
| `hhrap_2005_cattle_btf_report.pdf` | Cattle BTF report | Methodology for Predicting Cattle Biotransfer Factors (HHRAP support document) | - | https://archive.epa.gov/epawaste/hazard/tsd/td/web/pdf/btfreportfull05.pdf |

## Vadose zone transport (VLEACH)

| File | Short name | Document | Pages | Source |
|------|------------|----------|-------|--------|
| `vleach_2.2_users_guide.pdf` | VLEACH UG | Ravi, V. and Johnson, J.A., VLEACH: A One-Dimensional Finite Difference Vadose Zone Leaching Model, Version 2.2, EPA CSMOS | 71 | https://www.epa.gov/sites/production/files/2014-08/documents/vleach.pdf |
| `vleach_2.2a_setup.exe` | VLEACH setup | VLEACH 2.2a Windows setup (May 2007 release); contains the executable and the EPA sample problem used as reference (see `modules/vadose/tests/generate_vleach_reference.py`) | - | https://19january2025snapshot.epa.gov/sites/default/files/2014-07/vleachsetup.exe |

## Notes

- Part E links point to the EPA web snapshot of 19 January 2025; the
  document content is the final July 2004 guidance.
- When RSL tables are updated, add the new files with the new version
  suffix; do not overwrite older versions, so that past results remain
  reproducible.
