# OpenFold3 captured reference evidence

BIO-007, 2026-09-20. Native technical evidence only. The Biosimulant adapter, manifests and runtime tests are not yet implemented; no managed Hub run or release exists for this entry.

| Case | NVIDIA request ID | Result |
|---|---|---|
| openfold3-native-reference-001 | 827a8f89-1e77-4b2f-bcd4-506590ae81be | HTTP200; two protein samples, each 46 residues/326 atoms; full query identity and finite coordinates checked. |
| openfold3-rna-native-001 | a2363751-e3f4-4008-8359-9982f7166ade | HTTP200; one complete 12-base RNA sample; query-only MSA explicitly zero-hit. |
| openfold3-dna-native-001 | 42fc1d6f-3a15-40be-8f44-d4d8deb8d3e6 | HTTP504 after 32.169 seconds; exact existing UUID subsequently returned404. No result or live handle obtained; no duplicate submitted. |
| openfold3-ccd-native-001 | 8a9bb522-f92d-4779-84d8-07b5681fc510 | HTTP200; one ATP sample with all 31 heavy-atom names/elements matching the pinned CCD record. Requested chain L became A. |

The independent audit imports no adapter. It verifies original request/response hashes and lengths, echoed request/input IDs, sample counts/order, finite coordinates, unique atom keys, exact protein/RNA residue sequence and backbone atoms, native score domains and availability. It retains every original CIF sample. ATP mapping is derived from the single requested entity, unique CCD identity and exact heavy-atom name/element match; it is not a native chain-ID echo. Ligand bond orders are absent from the output and stereochemistry has not been independently validated. Full graph-dependent downstream compatibility remains unqualified.

All three successful response types have `iptm_score: null`, contradicting the published numeric-only hosted schema. Preserve unavailable interface confidence rather than substituting zero. Native source label `trtbnm` is not a deployment digest. No confidence value is called measured accuracy. Empty runtime_metrics is retained as empty, not an invented provider timing.

Full receipts and inspection reports are in the four named case folders under the tracking output directory. `openfold3-evidence-summary.json` summarizes the independent audit. Source records are complete public PDB1CRN/1ZIH/1BNA and CCD ATP; these are technical fixtures, not held-out scientific benchmarks.

Remaining: individually implement all proposed runtime ports and manifests; complete supported SMILES, DNA, paired, template, CSV/PDB and mixed-entity checks; inspect native/adapter parity and fresh execution; resolve chain multiplicity/renaming and chemical graph identity; pin exact deployment/artifact rights; establish reviewed benchmark thresholds; verify managed credentials/runtime and obtain immutable source/run approvals. No gate is passed solely because the direct NVIDIA endpoint returned HTTP200.
