# Data

`smarter_anonymised_data.csv`: anonymised individual-level data from the SMARTER cluster-randomised trial of a village doctor-led mobile health intervention for cardiovascular risk reduction in rural China.

- Source: Dryad Digital Repository, Li, X. (2025). doi:[10.5061/dryad.tmpg4f58w](https://doi.org/10.5061/dryad.tmpg4f58w)
- Trial publication: Zhang, X. et al. (2025). A village doctor-led mobile health intervention for cardiovascular risk reduction in rural China: cluster randomised controlled trial. *BMJ* 389, e082765.
- Trial registration: ClinicalTrials.gov NCT05645640.

Variables used in this repository:

| Variable | Meaning |
|---|---|
| `village_id` | village (unit of randomisation) |
| `group` | 1 = intervention, 2 = usual care (control) |
| `SBP0`, `SBP4` | systolic blood pressure (mmHg) at baseline and at 12 months |
| `risk0`, `risk4` | predicted 10-year ASCVD risk at baseline and at 12 months |

The analysis outcome is the SBP reduction `SBP0 - SBP4`, restricted to participants with a positive reduction (3917 participants in 127 villages).

If the data are removed from this folder, download the file from the Dryad link above and save it here under the same name.
