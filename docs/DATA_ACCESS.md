# Data, inputs and access boundaries

Provider source identifiers, as cited in the manuscript:

- MILK10k: https://doi.org/10.1016/j.jid.2025.06.1594
- MRA-MIDAS dataset: https://doi.org/10.71718/15nz-jv40 ; descriptor: https://doi.org/10.1056/AIdbp2400732
- HIBA: https://doi.org/10.1038/s41597-023-02630-0
- UQ: https://doi.org/10.1038/s41597-025-05880-2
- MCR-SL dataset: https://doi.org/10.5281/zenodo.17306338 ; descriptor: https://doi.org/10.3390/data10100166

Original images and externally licensed pretrained weights are not distributed here. Obtain them through the providers under their applicable access and reuse terms for image-level retraining.

The numerical modules contain aggregate cells and also coded record/group identifiers, diagnostic labels, frozen predictions/selections, and raw diagnostic-field slices. They support exact joins and grouped replay. This is therefore **not an aggregate-only release**, even though older module documentation described an earlier aggregate-only package. The repository is public. Provider-derived records remain subject to their source access and reuse terms; repository visibility does not assign a new license to those records.

No new clinical adjudication was performed. Missing confirmation-type metadata is reported as unspecified, not as proof of absent clinical assessment. Endpoint eligibility was fixed using existing records before replay; eligible diagnoses/subtypes are then withheld from the verifier until a group query. Excluded records are not assigned hypothetical cancer labels in compatible completions.

Historical validation JSON files can contain workstation paths and original timestamps. These are provenance, not required runtime paths. The root runner explicitly marks the original-workspace preservation check unavailable in working copies; scientific inputs and definitions remain unchanged. The quota script is adapted to packaged relative inputs and does not access its original recovery location.
