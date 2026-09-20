# OpenFold3 caveats

- Current work is an unqualified implementation design with native reference evidence. No managed Hub run or release has succeeded.
- Target NIM 1.6.0 and the actual hosted deployment are not proven equivalent. Hosted OpenAPI metadata incorrectly names AlphaFold2 and reports version 1.0.0; template fields and nullable scores differ from current NIM/observations.
- Exact weights, code, container, CCD and MSA database/index digests and artifact-specific rights are unresolved. Upstream preview weights are not assumed identical to the NVIDIA deployment.
- Native HTTP200, finite coordinates and confidence are not measured structural accuracy, biological function, affinity or clinical suitability.
- Historical 1CRN/1ZIH/1BNA fixtures may overlap training. Query-only RNA alignment contains no homology hits. Same-protein templates are not held-out evidence.
- Chain, residue and ligand atom correspondence must be preserved and independently reviewed. Ligand composition alone cannot establish bonds, stereochemistry or atom-index compatibility.
- All compatibility additions remain proposals until reviewed, released and supported by the actual local/managed runtimes. No Boltz score profile applies to OpenFold3 confidence.
- Handwritten scientific parsing/mapping and acceptance criteria require manual review. Missing benchmark thresholds block qualification instead of being assigned arbitrary values.
