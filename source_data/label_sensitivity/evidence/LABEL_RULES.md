# Label rules and their interpretation

The original MILK metadata contains paired image-level diagnoses that agree within each of the 5,240 lesion pairs. After the unchanged histopathology, overlap and duplicate rules, the additional 261 components consist of 99 explicitly named `Solar or actinic keratosis` and 162 `Keratoacanthoma` records. The frozen provider-derived binary field is 0 for the former and 1 for the latter. This experiment tests reinstating these records on those original binary sides; it does not re-adjudicate histopathology.

For the low-flag subtype heads, AK belongs to the binary-negative/background class. SCCIS remains its own class. KA follows the source SCC-family class. In source threshold calibration, the original SCC-family 95% floor therefore also applies to the reinstated KA cases. AK is never included in the SCCIS sensitivity denominator.

AK as an alternative non-cancer binary endpoint is supported by the description of AK as a precancerous lesion by the [Canadian Dermatology Association](https://dermatology.ca/public-patients/diseases-conditions/skin-conditions/actinic-keratoses/) and the [American Academy of Dermatology](https://www.aad.org/public/diseases/skin-cancer/actinic-keratosis-treatment). A non-cancer endpoint label is not a claim that the lesion requires no clinical attention.

KA classification should not be inferred solely from a broad `SCCKA` training category. Here the original fine diagnosis and original binary field are retained explicitly. [DermNet](https://dermnetnz.org/topics/keratoacanthoma) describes the diagnostic overlap with SCC and use of the SCC, KA-type terminology; the [WHO fifth-edition update](https://pmc.ncbi.nlm.nih.gov/articles/PMC12793825/) discusses criteria distinguishing KA from cutaneous SCC. The experiment's provider-aligned positive-side convention is consequently a declared sensitivity definition, not a new assertion that all KA are invasive SCC.

In the external frozen cohorts, the prior raw-field audit identified 92 explicit AK records and no explicitly labelled KA records. The original 107 other ambiguous/unspecified excluded endpoints are not assigned invented diagnoses. They remain excluded in `AK_noncancer` and are removed from the resource population in `endpoint_complete_case`.

SCCIS is already represented independently in the original complete criterion with a 90% floor. Renaming that same stratum IEC without changing membership or floor changes neither the integer deficits nor any feasibility margin. The equality check is reported as a mapping audit, not counted as a second empirical robustness experiment.

Sources consulted on 2026-11-09. The raw-field and image hashes, exact component labels, and all original fold assignments are retained in the experiment metadata and validation files.
