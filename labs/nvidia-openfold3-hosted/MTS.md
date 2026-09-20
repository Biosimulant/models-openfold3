# OpenFold3 technical design

BIO-007, 2026-09-20. Individually implemented local design; requirements review and managed runtime validation remain pending.

One `NvidiaOpenFold3` BioModule uses `ExecutionPolicy.ONCE_BEFORE_RUN`. It validates one complex, submits one canonical request to `https://health.api.nvidia.com/v1/biology/openfold/openfold3/predict`, observes the existing NVIDIA UUID, validates all returned samples and emits five terminal reports. The runtime timestamp is a computational boundary; the diffusion sample index is not biological time. Dependencies will be pinned to the actually tested environment, separately from local source checkout and the unverified managed runtime.

## Declared ports

All ports below are mixed report or identity envelopes with signal unit `1`. Every numerical scientific field carries its own units; no connection may interpret the whole envelope as one concentration or coordinate quantity. The manifest, Python SignalSpec and Lab mappings are implemented and locally tested. Scientific compatibility lifecycle remains Proposed pending review and target-runtime support.

| Direction | Port | Typed envelope | Scientific contract |
|---|---|---|---|
| Input | molecular_entities | record `{entities:json}` | Required ordered entity records with stable entity_id, type, explicit chain_ids, sequence or exactly one ligand SMILES/CCD specification, and source provenance. Optional CCD chemical_reference is original ccd_cif text, sha256 and source. No default entity. |
| Input | alignments | record `{entries:json}` | Optional empty at signal layer, conditionally required for protein/RNA. Each entry names entity_id, unpaired/paired kind, database, format, unchanged text, sha256 and provenance. |
| Input | templates | record `{entries:json}` | Optional empty list. Each entry names protein entity, original CIF text/hash, source/name, optional native chain_id and provenance. |
| Input | run_options | record `{parameters:json}` | Defaults diffusion_samples=1 and output_format=cif. No exposed native seed is invented. |
| Input | query_metadata | record `{query_id:str,request_id:str,provenance:json}` | Explicit caller IDs and provenance; deterministic content-derived fallback IDs if blank. |
| Output | predicted_complex | record `{entities:json,structures:json}` | All native ranked samples; original coordinate text/format/hash, explicit entity-to-chain/residue/atom maps and mapping limitations. |
| Output | confidence_scores | record `{records:json,semantics:json}` | Native per-sample ranking, complex pLDDT/PDE, pTM/ipTM and explicit availability; preserve named ModelCIF QA metrics separately. |
| Output | inference_metrics | record `{native:json,semantics:json}` | Original native runtime metrics; `{}` or null stays unavailable rather than a fabricated duration. |
| Output | run_provenance | record `{request:json,provider:json,inputs:json,backend_identity:json,limitations:json}` | Exact request/response hashes, provider ID, original inputs/fixture checks, observation bounds and version gaps. |
| Output | native_response | scalar str | Original successful JSON response. Error bodies remain failed-attempt evidence rather than scientific signals. |

## Entity and artifact validation

The hosted schema permits one input, 1–32 molecule records and sequence length 2–4096. The card says ≤1000 and current NIM examples say 1–4096: retain this discrepancy and do not imply that a small successful fixture validates the full domain. The wrapper additionally bounds total chain copies to32 and total polymer residues across copies to4096; these are client operational bounds, not NVIDIA accuracy/capacity guarantees. Chain IDs must be unique and native-valid. PDB output can represent only single-character chain IDs, while CIF is the preferred lossless identity representation. Standard protein, DNA and RNA alphabets are validated by type. Nucleic acid sequences use 5′→3′ order; double-stranded DNA requires both strands explicitly.

Protein and RNA alignments must begin with the exact ungapped query. Preserve A3M insertions, row headers and column mapping; CSV keys require a separate parser and pairing semantics. Maximum three native database entries per alignment kind. Query-only rows are accepted only when explicitly supplied and labeled zero-hit. The BIO-005 collection needs an entity/chain selector and mapping adapter; paired collection rows are provider-asserted unless species evidence establishes them. No universal biological pairing claim follows from equal lengths or row counts.

Protein templates require parseable original CIF, source hash and selected/available chain identity. Current NIM documents identity×coverage selection and configurable minimum score; the client must not invent its own alignment or pretend that the selected template chain was returned. Hosted schema omits structural_templates, so capability and actual-use evidence are required. Compressed/PDB templates and templates for DNA/RNA/ligands remain unsupported native contracts.

For coordinate outputs, parse every sample and preserve original text. Check finite xyz values in angstrom, unique atom keys, model/chain multiplicity, polymer residue identity/order and required backbone atoms. Do not infer ligand atom order from SMILES order. CCD ligand atom-name/element correspondence needs a pinned chemical component reference; SMILES requires a reviewed graph/stereochemistry adapter. Missing bond-order/stereochemistry evidence stays unresolved and blocks qualified ligand-compatible wiring. Never infer chemical bonds solely from distances without an explicit scientifically reviewed method.

## Confidence and observations

Native fields are distinct: confidence_score is sample ranking (formula/scale unresolved), complex_plddt_score is mean predicted local quality on 0–100 scale, complex_pde_score is predicted distance-error score (exact aggregation/units require upstream verification), ptm_score and iptm_score are dimensionless predicted TM-style quantities. Actual monomer iptm_score is null despite the published numeric schema. Preserve null and applicability. PAE/PDE are different; never relabel one as the other. Native coordinate B-factor fields are not experimental uncertainty; explicit ModelCIF metric declarations provide a separate basis for interpretation.

The request uses canonical JSON bytes with SHA-256. Credentials are read only from secure worker `NVIDIA_API_KEY`; no `.env.agent` content is embedded. HTTPS destinations are fixed NVIDIA endpoints, redirects disabled. One POST only; HTTP202 leads to bounded GET polling of the exact UUID. Maximum client observation and response byte limits are explicit. Failed observation preserves UUID/hash and last safe status; no automatic resubmission. Resume requires both UUID and matching original request digest. Native 413/422/429/500/503 are captured as distinct failed outcomes.

## Verification and remaining work

Native protein reference uses complete 1CRN chain A and genuine 101-row BIO-005 A3M, two samples, CIF output. Native RNA fixture uses complete 12-base 1ZIH sequence with explicit zero-hit query-only alignment. DNA fixture uses complete two-strand 1BNA sequence and explicit A/B chain copies. These old structures can overlap training; they are technical fixtures, not held-out benchmarks. Ligand, paired, template, CSV, PDB and mixed-entity cases require separate evidence. Exact native/adapter request parity and fresh inference comparison must be measured before claiming local runtime success.

There are66 passing local tests covering canonical manifest/SignalSpec/Lab parity; conditional inputs; type/shape/unit/identity failures; malformed native output; null confidence; sample ordering; source hashes; transport bounds; no automatic duplicate submission; chemical composition loss; standard side-chain completeness; single CCD nucleotides versus polymers; and all five terminal BioWorld outputs. Dependencies are biosimulant0.0.34, httpx0.28.1, biopython1.84, numpy1.26.4 and rdkit2025.3.3 on Python3.10.20, exactly pinned in the manifest. The last two parser refinements are checked against captured native responses and unit tests; live-run receipts retain the exact earlier source hashes they executed. Scientific acceptance remains blocked by absent reviewed benchmark panel and tolerances. Managed credentials, immutable source/run approval and exact release identity remain independent gates.

Original sequence/molecule identity is checked before output acceptance. Standard protein side-chain heavy atoms and nucleic-acid sugar/base atoms are required; internal nucleotides also require phosphate atoms. Terminal chemistry such as protein OXT and first-nucleotide phosphate is not resolved by this contract and needs separate reference evidence. Copy mappings are candidate sets checked for a complete one-to-one assignment; ambiguous identical copies/isomers are retained rather than arbitrarily assigned. A ligand matching only element counts remains chemically unqualified. Actual ATP SMILES output lost atoms and fails even the composition prerequisite.

Official sources are listed in MRS.md. Current hosted response identifies `source: trtbnm`; that string is not a code, checkpoint or container digest.
