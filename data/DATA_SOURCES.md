# Data Sources

All GWAS summary statistics used in this study are publicly available, but the full files are not redistributed in this package. Download them from the original sources and place them in a `gwas_data/` directory before running manuscript-scale analyses. The package includes only the small demo subset in `../demo/`.

## Exposure GWAS (Lipid Traits)

| Source | Population | N | Traits | Build | URL | Reference |
|--------|-----------|---|--------|-------|-----|-----------|
| KoGES | Korean | 72,298 | LDL, TC, TG, Total lipids | GRCh37 | https://pheweb.snu.ac.kr | Kim et al. 2023, Cell Genomics |
| BBJ | Japanese | ~260,000 | LDL-C, TC, TG, HDL-C | GRCh37 | https://biobankportal.hgc.jp | Sakaue et al. 2021, Nat Genet |
| TPMI | Taiwanese | 142,676–165,202 | LDL-C, TC, HDL-C | GRCh38 | https://tpmi.ibms.sinica.edu.tw | Yang et al. 2025, Nature |

## Outcome GWAS (Glaucoma)

| Source | Population | N | Cases | Phenotype | Build | URL | Reference |
|--------|-----------|---|-------|-----------|-------|-----|-----------|
| BBJ Glaucoma | Japanese | 177,351 | 5,765 | ICD glaucoma | GRCh37 | https://biobankportal.hgc.jp | Sakaue et al. 2021 |
| TPMI Glaucoma | Taiwanese | 318,597 | 12,092 | Phecode 365 | GRCh38 | https://tpmi.ibms.sinica.edu.tw | Yang et al. 2025 |
| FinnGen POAG | Finnish | 342,499 | — | H7_GLAUCOMA_POAG | GRCh37 | https://r8.finngen.fi | Kurki et al. 2023, Nature |
| FinnGen NTG | Finnish | 342,499 | — | H7_GLAUCOMA_NTG | GRCh37 | https://r8.finngen.fi | Kurki et al. 2023 |
| FinnGen XFG | Finnish | 342,499 | — | H7_GLAUCOMA_XFG | GRCh37 | https://r8.finngen.fi | Kurki et al. 2023 |
| Luben PACG | Multi-ancestry | 797,502 | 9,217 | PACG | GRCh38 | GWAS Catalog GCST90832110 | Luben et al. 2025, Nat Commun |
| GBMI POAG | EAS meta | 269,960 | 9,715 | POAG | GRCh38 | https://globalbiobankmeta.org | Faro et al. 2024, Cell Rep Med |

## LD Reference Panel

| Resource | Description | URL |
|----------|-------------|-----|
| 1000 Genomes Phase 3 EAS | LD reference for clumping (504 EAS samples) | ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ |

## Standardized File Format

All GWAS files should be standardized to TSV with these columns:

```
SNP  effect_allele  other_allele  beta  se  pval  eaf  n  trait  source  CHR  POS
```

- `beta`: per-allele effect size (log-OR for binary traits)
- `se`: standard error of beta
- `eaf`: effect allele frequency
- `n`: sample size (can be NA if constant across variants)
