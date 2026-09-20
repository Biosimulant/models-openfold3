# OpenFold3 captured reference evidence

BIO-007, 2026-09-20. Individually authored adapter, manifests and66 local tests are complete. Native and local BioWorld technical evidence exists; managed runtime and scientific qualification remain incomplete.

| Case | NVIDIA request ID | Result |
|---|---|---|
| openfold3-native-reference-001 | 827a8f89-1e77-4b2f-bcd4-506590ae81be | HTTP200; two protein samples, each 46 residues/326 atoms; full query identity and finite coordinates checked. |
| openfold3-rna-native-001 | a2363751-e3f4-4008-8359-9982f7166ade | HTTP200; one complete 12-base RNA sample; query-only MSA explicitly zero-hit. |
| openfold3-dna-native-001 | 42fc1d6f-3a15-40be-8f44-d4d8deb8d3e6 | HTTP504 after 32.169 seconds; exact existing UUID subsequently returned404. No result or live handle obtained; no duplicate submitted. |
| openfold3-ccd-native-001 | 8a9bb522-f92d-4779-84d8-07b5681fc510 | HTTP200; one ATP sample with all 31 heavy-atom names/elements matching the pinned CCD record. Requested chain L became A. |
| openfold3-smiles-native-001 | 7aed7c05-6d94-4b2a-b30a-7e9809698d25 | HTTP200 but failed chemical composition: input C10N5O13P3, output C10N5O9P2. Adapter rejects it. |
| openfold3-adapter-live-001 | 5ec77dce-1064-4837-99f7-4538877eabd7 | BioWorld completed, but evidence driver used unsupported signal.unit attribute and failed before saving results. Capture failure retained. |
| openfold3-adapter-live-002 | f7b75219-bdcc-4bae-9666-93059077bde8 | Protein BioWorld execution; five terminal outputs durably saved; native response byte-equal to reference. |
| openfold3-rna-adapter-001 | 778944a4-564b-411d-b7ec-fca33ba8381a | RNA BioWorld execution; five terminal outputs saved; native response byte-equal to reference. |
| openfold3-ccd-adapter-001 | 2318363d-1e87-49f1-b877-f146c06e186e | CCD ATP BioWorld execution; five outputs saved; native response byte-equal to reference. |
| openfold3-adapter-live-003 | 93123ff0-460b-4c30-aeea-bc9ac0d6b917 | Protein BioWorld with stronger standard-heavy-atom checks; all outputs saved and native response byte-equal to reference. |

The independent audit imports no adapter. It verifies original request/response hashes and lengths, echoed request/input IDs, sample counts/order, finite coordinates, unique atom keys, exact protein/RNA residue sequence and backbone atoms, native score domains and availability. It retains every original CIF sample. ATP mapping is derived from the single requested entity, unique CCD identity and exact heavy-atom name/element match; it is not a native chain-ID echo. Ligand bond orders are absent from the output and stereochemistry has not been independently validated. Full graph-dependent downstream compatibility remains unqualified.

All three successful response types have `iptm_score: null`, contradicting the published numeric-only hosted schema. Preserve unavailable interface confidence rather than substituting zero. Native source label `trtbnm` is not a deployment digest. No confidence value is called measured accuracy. Empty runtime_metrics is retained as empty, not an invented provider timing.

Full receipts and inspection reports are in the ten named case folders under the tracking output directory. `openfold3-evidence-summary.json` summarizes the independent audit, including five terminal signals and raw/result hashes where obtained. `openfold3-adapter-contract-summary.json` separately verifies byte-exact request parity for all five native cases and current parser acceptance/rejection; that parity check imports the adapter and is not independent algorithm verification. Source records are complete public PDB1CRN/1ZIH/1BNA and CCD ATP; these are technical fixtures, not held-out scientific benchmarks.

Current source SHA-256: `580075d026bb3ea4055de51b45c6763803070f071ae27f9bbe4dfed7f876cabb`. Protein/RNA/CCD durable live calls first used `1fde08e4a26e0ffa6534b0185a66bd7590a4069a9ae6fc826891a653ac657108`; stronger atom-check protein live003 used `6a998738c3a5435183232f8fb317b15748df33b5be0b4ea64fa03bf7e4c1502c`. The final single-CCD-versus-polymer classification refinement is tested locally and against the original captured outputs; do not attribute earlier live runs to a later source hash.

Remaining: resolve provider SMILES composition and DNA failures; complete paired/template/CSV/PDB/mixed-entity native integration; verify chemical graph/stereochemistry and ambiguous copy mapping for qualified consumers; pin exact deployment/artifact rights; establish reviewed benchmark criteria; verify managed credentials/runtime and immutable source/run approvals. No gate is passed solely because the NVIDIA endpoint returned HTTP200.
