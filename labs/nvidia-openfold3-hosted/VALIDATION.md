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
| openfold3-pdb-native-001 | c198576b-8aa5-4c80-aabf-d3710b114520 | Two PDB samples pass independent fixed-column identity and coordinate checks, plus current adapter parsing. |
| openfold3-csv-native-001 | a01c64ca-b503-4e27-8713-4c97d768b5d4 | Two CIF samples accepted from CSV MSA. All 101 original A3M sequences, insertion/gap characters and order preserved; original headers mapped to local CSV keys. |
| openfold3-template-native-001 | 1eb0ccb5-5f4c-41a8-8d83-ec334c979b69 | Two CIF samples accepted with the complete pinned 1CRN template. Provider template use is not attested; this self-template fixture cannot measure predictive accuracy. |
| openfold3-pdb-adapter-001 | 505491a7-69ad-4ed0-81d9-5bf8356db425 | Current source BioWorld PDB execution; five terminal outputs saved and native response byte-equal to PDB reference. |

The independent audit imports no adapter. It verifies original request/response hashes and lengths, echoed request/input IDs, sample counts/order, finite coordinates, unique atom keys, exact protein/RNA residue sequence and backbone atoms, native score domains and availability. It retains every original CIF sample. ATP mapping is derived from the single requested entity, unique CCD identity and exact heavy-atom name/element match; it is not a native chain-ID echo. Ligand bond orders are absent from the output and stereochemistry has not been independently validated. Full graph-dependent downstream compatibility remains unqualified.

All three successful response types have `iptm_score: null`, contradicting the published numeric-only hosted schema. Preserve unavailable interface confidence rather than substituting zero. Native source label `trtbnm` is not a deployment digest. No confidence value is called measured accuracy. Empty runtime_metrics is retained as empty, not an invented provider timing.

Full receipts and inspection reports are in the fourteen named case folders under the tracking output directory. Eleven attempts pass their recorded technical checks, two fail native interface/provider checks, and one fails durable capture after inference. `openfold3-evidence-summary.json` summarizes the independent audit, including five terminal signals and raw/result hashes where obtained. `openfold3-adapter-contract-summary.json` separately verifies byte-exact request parity for all eight native cases and current parser acceptance/rejection; that parity check imports the adapter and is not independent algorithm verification. Source records are complete public PDB1CRN/1ZIH/1BNA and CCD ATP; these are technical fixtures, not held-out scientific benchmarks.

Current source SHA-256: `580075d026bb3ea4055de51b45c6763803070f071ae27f9bbe4dfed7f876cabb`. Protein/RNA/CCD durable live calls first used `1fde08e4a26e0ffa6534b0185a66bd7590a4069a9ae6fc826891a653ac657108`; stronger atom-check protein live003 used `6a998738c3a5435183232f8fb317b15748df33b5be0b4ea64fa03bf7e4c1502c`. The final single-CCD-versus-polymer classification refinement is tested locally and against the original captured outputs; do not attribute earlier live runs to a later source hash.

The PDB BioWorld call uses current source `580075d026bb3ea4055de51b45c6763803070f071ae27f9bbe4dfed7f876cabb`. It saves results SHA-256 `18a3e31cde1ea24c3a4108784078827c860932d4d96a628b39a0d277fd17868a` (507846 bytes). CSV keys are local row identifiers, not taxonomic associations. The exact 1CRN template is SHA-256 `23787562c427d7c1abe5420e86d5f1d0a6c7007dec1e8ce85645a6d69c32e8ba`. Successful template request handling does not establish causal template use.

Remaining: resolve provider SMILES composition and DNA failures; complete paired/mixed-entity native integration and verify template consumption; verify chemical graph/stereochemistry and ambiguous copy mapping for qualified consumers; pin exact deployment/artifact rights; establish reviewed benchmark criteria; verify managed credentials/runtime and immutable source/run approvals. No gate is passed solely because the NVIDIA endpoint returned HTTP200.

Queue decision: explicitly parked after these fourteen attempts. Source plan `5e2abaac-f409-430c-92c1-ccddcea2a04e` expired at 2026-09-20 03:18:30 UTC without exact confirmation or application. It covers earlier documents, not this evidence update. No managed compute or release occurred. Resume once a concrete unresolved gate can be addressed; do not renew an expired plan merely to wait again.
