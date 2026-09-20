"""Individually authored NVIDIA OpenFold3 adapter; qualification is separate."""
import collections,csv,hashlib,io,json,math,os,re,time,uuid
from datetime import datetime,timezone
import httpx
from Bio.PDB import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.PDBData import protein_letters_3to1
from rdkit import Chem
from biosim import BioModule,ExecutionPolicy,SignalSpec

ENDPOINT='https://health.api.nvidia.com/v1/biology/openfold/openfold3/predict'
POLL='https://api.nvcf.nvidia.com/v2/nvcf/pexec/status/'
AA='ACDEFGHIKLMNPQRSTVWY'
SIDECHAIN={
 'ALA':'CB','ARG':'CB CG CD NE CZ NH1 NH2','ASN':'CB CG OD1 ND2','ASP':'CB CG OD1 OD2',
 'CYS':'CB SG','GLN':'CB CG CD OE1 NE2','GLU':'CB CG CD OE1 OE2','GLY':'',
 'HIS':'CB CG ND1 CD2 CE1 NE2','ILE':'CB CG1 CG2 CD1','LEU':'CB CG CD1 CD2',
 'LYS':'CB CG CD CE NZ','MET':'CB CG SD CE','PHE':'CB CG CD1 CD2 CE1 CE2 CZ',
 'PRO':'CB CG CD','SER':'CB OG','THR':'CB OG1 CG2','TRP':'CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
 'TYR':'CB CG CD1 CD2 CE1 CE2 CZ OH','VAL':'CB CG1 CG2'}
BASE_ATOMS={'A':'N9 C8 N7 C5 C6 N6 N1 C2 N3 C4','G':'N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4','C':'N1 C2 O2 N3 C4 N4 C5 C6','U':'N1 C2 O2 N3 C4 O4 C5 C6','T':'N1 C2 O2 N3 C4 O4 C5 C6 C7'}
ENTITIES={'entities':'json'};ENTRIES={'entries':'json'};OPTIONS={'parameters':'json'}
META={'query_id':'str','request_id':'str','provenance':'json'}
STRUCTURES={'entities':'json','structures':'json'};REPORT={'records':'json','semantics':'json'}
METRICS={'native':'json','semantics':'json'}
PROVENANCE={'request':'json','provider':'json','inputs':'json','backend_identity':'json','limitations':'json'}

def require(condition,message):
    if not condition:raise ValueError(message)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def encode(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
def finite(x):return type(x) in (int,float) and math.isfinite(x)
def exact(value,keys,label):require(isinstance(value,dict) and set(value)==set(keys),label+' schema mismatch')
def strict_json(raw):
    def pairs(items):
        d={}
        for key,value in items:require(key not in d,'Duplicate JSON key');d[key]=value
        return d
    value=json.loads(raw,object_pairs_hook=pairs);encode(value);return value
def identifier(value,limit=128):return isinstance(value,str) and 0<len(value)<=limit and not any(ord(c)<32 for c in value)
def table(d,prefix,columns):
    require(all(prefix+c in d for c in columns),'Missing CIF columns')
    values=[d[prefix+c] for c in columns];require(len({len(x) for x in values})==1,'CIF column length mismatch')
    return [dict(zip(columns,row)) for row in zip(*values)]

def check_alignment(text,fmt,sequence,kind):
    require(isinstance(text,str) and 0<len(text.encode())<=5000000,'Alignment byte bound')
    rows=[];keys=[]
    if fmt=='a3m':
        current=None;parts=[]
        for line in text.splitlines():
            if not line:continue
            if line.startswith('>'):
                if current is not None:require(bool(parts),'Empty alignment row');rows.append(''.join(parts));keys.append(current)
                current=line[1:];parts=[];require(bool(current.strip()) and len(current)<=4096 and not any(ord(c)<32 and c!='\t' for c in current),'Invalid alignment header')
                require(not current.startswith(('ss_pred','ss_conf','ss_dssp','sa_')),'Annotation track requires adapter')
            else:require(current is not None and not any(c.isspace() for c in line),'Malformed A3M');parts.append(line)
        require(current is not None and parts,'Empty alignment');rows.append(''.join(parts));keys.append(current)
    elif fmt=='csv':
        reader=csv.DictReader(io.StringIO(text));require(reader.fieldnames==['key','sequence'],'Expected CSV key,sequence header')
        for row in reader:
            require(set(row)=={'key','sequence'} and isinstance(row['sequence'],str) and identifier(row['key'],256),'Malformed CSV row')
            keys.append(row['key']);rows.append(row['sequence'])
    else:raise ValueError('Only native A3M and reviewed key/sequence CSV supported')
    require(1<=len(rows)<=20000 and rows[0]==sequence,'Alignment query/row count mismatch')
    alphabet=AA+'XBZJUO' if kind=='protein' else 'ACGUN'
    allowed=set(alphabet+alphabet.lower()+'-.')
    for row in rows:
        require(set(row)<=allowed,'Alignment alphabet mismatch')
        require(len(''.join(c for c in row if not c.islower() and c!='.'))==len(sequence),'Alignment match columns differ from query')
    return {'rows':len(rows),'query_only':len(rows)==1,'row_keys':keys,'sha256':sha(text.encode()),'bytes':len(text.encode()),'format':fmt}

def chemical_reference(value,code):
    if value is None:return None
    exact(value,('format','text','sha256','source'),'Chemical reference')
    text=value['text'];require(value['format']=='ccd_cif' and isinstance(text,str) and 0<len(text.encode())<=2000000,'Expected bounded CCD CIF')
    require(value['sha256']==sha(text.encode()) and identifier(value['source'],2048),'Chemical reference hash/source mismatch')
    d=MMCIF2Dict(io.StringIO(text));require(d.get('_chem_comp.id')==[code],'Chemical reference CCD mismatch')
    atoms=table(d,'_chem_comp_atom.',['atom_id','type_symbol']);heavy={x['atom_id']:x['type_symbol'].upper() for x in atoms if x['type_symbol'].upper() not in ('H','D')}
    require(heavy and len({x['atom_id'] for x in atoms})==len(atoms),'Duplicate/empty CCD atoms')
    bonds=table(d,'_chem_comp_bond.',['atom_id_1','atom_id_2','value_order']) if '_chem_comp_bond.atom_id_1' in d else []
    require(all(b['atom_id_1'] in {x['atom_id'] for x in atoms} and b['atom_id_2'] in {x['atom_id'] for x in atoms} for b in bonds),'CCD bond atom missing')
    return {'sha256':value['sha256'],'source':value['source'],'heavy_atoms':heavy,'reference_bonds':bonds,'stereochemistry_validated':False}

def prepare_request(entities,alignments,templates,options,metadata):
    exact(entities,ENTITIES,'Entities');items=entities['entities'];require(isinstance(items,list) and 1<=len(items)<=32,'Expected 1–32 entity records')
    exact(alignments,ENTRIES,'Alignments');exact(templates,ENTRIES,'Templates');exact(options,OPTIONS,'Options');exact(metadata,META,'Metadata')
    require(isinstance(metadata['provenance'],dict),'Provenance must be object')
    for key in ('query_id','request_id'):require(metadata[key]=='' or identifier(metadata[key]),'Invalid request/query ID')
    settings=options['parameters'];require(isinstance(settings,dict) and set(settings)<= {'diffusion_samples','output_format'},'Unknown run option')
    settings={'diffusion_samples':1,'output_format':'cif',**settings}
    require(type(settings['diffusion_samples']) is int and 1<=settings['diffusion_samples']<=5,'Invalid sample count')
    require(settings['output_format'] in ('cif','pdb'),'Invalid output format')
    native=[];identity=[];by_id={};all_chains=set();total_residues=0
    for item in items:
        require(isinstance(item,dict) and set(item)<= {'entity_id','type','chain_ids','sequence','smiles','ccd_codes','provenance','chemical_reference'},'Unknown entity field')
        require({'entity_id','type','chain_ids','provenance'}<=set(item),'Missing entity identity')
        eid=item['entity_id'];kind=item['type'];chains=item['chain_ids']
        require(identifier(eid) and eid not in by_id and kind in ('protein','dna','rna','ligand'),'Invalid entity ID/type')
        require(isinstance(chains,list) and 1<=len(chains)<=32 and all(isinstance(x,str) and re.fullmatch('[A-Za-z0-9]{1,4}',x) for x in chains),'Invalid explicit chain IDs')
        require(len(set(chains))==len(chains) and not all_chains.intersection(chains),'Duplicate chain ID');all_chains.update(chains)
        require(len(all_chains)<=32 and isinstance(item['provenance'],dict),'Too many copies or invalid provenance')
        require(settings['output_format']!='pdb' or all(len(x)==1 for x in chains),'PDB requires one-character chain IDs')
        molecule={'type':kind,'id':chains[0] if len(chains)==1 else chains}
        info={'entity_id':eid,'type':kind,'chain_ids':chains,'provenance':item['provenance']}
        if kind=='ligand':
            require('sequence' not in item and (('smiles' in item)^('ccd_codes' in item)),'Ligand requires exactly one SMILES/CCD')
            if 'smiles' in item:
                require('chemical_reference' not in item,'SMILES chemical reference requires a reviewed graph adapter')
                smiles=item['smiles'];require(isinstance(smiles,str) and 0<len(smiles)<=20000,'Invalid SMILES size');mol=Chem.MolFromSmiles(smiles)
                require(mol is not None and 1<=mol.GetNumHeavyAtoms()<=2000,'Invalid ligand SMILES')
                molecule['smiles']=smiles;info.update({'smiles':smiles,'canonical_isomeric_smiles':Chem.MolToSmiles(mol,isomericSmiles=True),'heavy_element_counts':dict(collections.Counter(a.GetSymbol().upper() for a in mol.GetAtoms() if a.GetAtomicNum()>1)),'atom_graph_mapping':'Unverified until native atom graph correspondence is established.'})
            else:
                code=item['ccd_codes'];require(isinstance(code,str) and re.fullmatch('[A-Z0-9]{1,5}',code),'Invalid CCD code')
                molecule['ccd_codes']=code;info.update({'ccd_codes':code,'chemical_reference':chemical_reference(item.get('chemical_reference'),code)})
        else:
            require('sequence' in item and not set(item).intersection(('smiles','ccd_codes','chemical_reference')),'Polymer entity fields mismatch')
            seq=item['sequence'];alphabet=AA if kind=='protein' else 'ACGT' if kind=='dna' else 'ACGU'
            require(isinstance(seq,str) and 2<=len(seq)<=4096 and set(seq)<=set(alphabet),'Invalid standard polymer sequence')
            total_residues+=len(seq)*len(chains);require(total_residues<=4096,'Client total polymer residue bound exceeded')
            molecule['sequence']=seq;info.update({'sequence':seq,'sequence_sha256':sha(seq.encode()),'length':len(seq)})
        by_id[eid]=len(native);native.append(molecule);identity.append(info)
    require(isinstance(alignments['entries'],list) and len(alignments['entries'])<=192,'Invalid alignment entries')
    alignment_info=[];seen=set();paired=[]
    for entry in alignments['entries']:
        exact(entry,('entity_id','kind','database','format','text','sha256','provenance'),'Alignment entry')
        require(entry['entity_id'] in by_id and entry['kind'] in ('unpaired','paired'),'Unknown alignment entity/kind')
        index=by_id[entry['entity_id']];entity=native[index];kind=entity['type'];field='msa' if entry['kind']=='unpaired' else 'paired_msa'
        require(kind in ('protein','rna') and not (kind=='rna' and field=='paired_msa'),'Unsupported alignment type')
        db=entry['database'];fmt=entry['format'];key=(index,field,db,fmt)
        require(identifier(db) and key not in seen and isinstance(entry['provenance'],dict),'Duplicate alignment or invalid provenance');seen.add(key)
        checked=check_alignment(entry['text'],fmt,entity['sequence'],kind);require(entry['sha256']==checked['sha256'],'Alignment checksum mismatch')
        entity.setdefault(field,{}).setdefault(db,{})[fmt]={'alignment':entry['text'],'format':fmt};require(len(entity[field])<=3,'Too many MSA databases')
        record={'entity_id':entry['entity_id'],'kind':entry['kind'],'database':db,'provenance':entry['provenance'],**checked};alignment_info.append(record)
        if field=='paired_msa':paired.append(record)
    for entity in native:
        if entity['type']=='protein':require(bool(entity.get('msa') or entity.get('paired_msa')),'Protein requires explicit alignment')
        if entity['type']=='rna':require(bool(entity.get('msa')),'RNA requires explicit alignment')
    if paired:
        require(len({x['entity_id'] for x in paired})>=2,'Paired alignment requires multiple explicit entities')
        for db in {x['database'] for x in paired}:
            group=[x for x in paired if x['database']==db];require(len({x['rows'] for x in group})==1,'Paired row count mismatch')
    require(isinstance(templates['entries'],list) and len(templates['entries'])<=32,'Invalid template entries')
    template_info=[]
    for entry in templates['entries']:
        require(isinstance(entry,dict) and set(entry)<= {'entity_id','text','sha256','source','name','chain_id','provenance'} and {'entity_id','text','sha256','source','name','provenance'}<=set(entry),'Template schema mismatch')
        require(entry['entity_id'] in by_id,'Unknown template entity');entity=native[by_id[entry['entity_id']]];require(entity['type']=='protein','Templates supported for protein only')
        text=entry['text'];require(isinstance(text,str) and 0<len(text.encode())<=2000000 and sha(text.encode())==entry['sha256'],'Template size/hash mismatch')
        require(identifier(entry['name']) and identifier(entry['source'],2048) and isinstance(entry['provenance'],dict),'Invalid template provenance')
        d=MMCIF2Dict(io.StringIO(text));atoms=table(d,'_atom_site.',['label_asym_id','label_comp_id','Cartn_x','Cartn_y','Cartn_z']);require(atoms and all(math.isfinite(float(a[k])) for a in atoms for k in ('Cartn_x','Cartn_y','Cartn_z')),'Invalid template coordinates')
        record={'structure':text,'format':'cif','name':entry['name']}
        if 'chain_id' in entry:
            cid=entry['chain_id'];require(isinstance(cid,str) and re.fullmatch('[A-Za-z0-9]{1,10}',cid) and cid in {a['label_asym_id'] for a in atoms},'Template chain unavailable');record['chain_id']=cid
        entity.setdefault('structural_templates',[]).append(record);require(len(entity['structural_templates'])<=10,'Too many templates per entity')
        template_info.append({k:v for k,v in entry.items() if k!='text'})
    core={'molecules':native,**settings};fallback='sha256:'+sha(encode(core));core['input_id']=metadata['query_id'] or fallback
    request={'request_id':metadata['request_id'] or fallback,'inputs':[core]};require(len(encode(request))<=10000000,'Request byte bound exceeded')
    return request,{'entities':identity,'alignments':alignment_info,'templates':template_info,'metadata':metadata,'paired_association':'Provider-asserted row association; species identity not independently verified.' if paired else None}

def read_atoms(text,fmt):
    require(isinstance(text,str) and 0<len(text.encode())<=5000000,'Invalid structure size')
    atoms=[];qa={}
    if fmt=='cif':
        d=MMCIF2Dict(io.StringIO(text))
        columns=['label_asym_id','auth_asym_id','label_seq_id','auth_seq_id','label_comp_id','label_atom_id','type_symbol','Cartn_x','Cartn_y','Cartn_z','pdbx_PDB_model_num']
        rows=table(d,'_atom_site.',columns)
        for i,a in enumerate(rows):
            alt=d.get('_atom_site.label_alt_id',['.']*len(rows))[i];require(alt in ('.','?',''),'Alternate atom location needs reviewed mapping')
            ins=d.get('_atom_site.pdbx_PDB_ins_code',['?']*len(rows))[i]
            atoms.append({'chain':a['label_asym_id'],'auth_chain':a['auth_asym_id'],'residue_id':a['label_seq_id'] if a['label_seq_id'] not in ('.','?') else a['auth_seq_id'],'auth_residue_id':a['auth_seq_id'],'insertion_code':ins,'monomer':a['label_comp_id'],'atom_name':a['label_atom_id'],'element':a['type_symbol'].upper(),'xyz':[float(a[k]) for k in ('Cartn_x','Cartn_y','Cartn_z')],'model_id':a['pdbx_PDB_model_num']})
        qa={k:v for k,v in d.items() if k.startswith(('_ma_qa_metric.','_ma_qa_metric_local.','_ma_qa_metric_global.'))}
    elif fmt=='pdb':
        parsed=PDBParser(QUIET=True,PERMISSIVE=False).get_structure('prediction',io.StringIO(text));models=list(parsed.get_models());require(len(models)==1,'Expected one coordinate model per sample')
        for chain in models[0]:
            for residue in chain:
                require(not residue.is_disordered(),'Disordered residue needs adapter')
                for atom in residue:
                    require(not atom.is_disordered(),'Disordered atom needs adapter')
                    atoms.append({'chain':chain.id,'auth_chain':chain.id,'residue_id':str(residue.id[1]),'auth_residue_id':str(residue.id[1]),'insertion_code':residue.id[2].strip(),'monomer':residue.resname,'atom_name':atom.name,'element':atom.element.upper(),'xyz':[float(x) for x in atom.coord],'model_id':str(models[0].id)})
    else:raise ValueError('Unsupported native structure format')
    require(atoms and len(atoms)<=200000 and len({a['model_id'] for a in atoms})==1,'Empty/oversized/multiple-model structure')
    require(all(all(math.isfinite(x) for x in a['xyz']) for a in atoms),'Nonfinite coordinate')
    keys=[(a['chain'],a['residue_id'],a['insertion_code'],a['monomer'],a['atom_name']) for a in atoms];require(len(keys)==len(set(keys)),'Duplicate atom identity')
    return atoms,qa

def map_structure(atoms,entities):
    chains=collections.OrderedDict()
    for atom in atoms:chains.setdefault(atom['chain'],[]).append(atom)
    require(len(chains)==sum(len(x['chain_ids']) for x in entities),'Chain multiplicity mismatch')
    candidates={};descriptions={}
    for chain,group in chains.items():
        residues=collections.OrderedDict()
        for atom in group:residues.setdefault((atom['residue_id'],atom['insertion_code']),[]).append(atom)
        monomers=[]
        for residue in residues.values():
            kinds={a['monomer'] for a in residue};require(len(kinds)==1,'Conflicting residue monomers');monomers.append(next(iter(kinds)))
        if len(monomers)>=2 and all(m in protein_letters_3to1 for m in monomers):kind='protein';seq=''.join(protein_letters_3to1[m] for m in monomers)
        elif len(monomers)>=2 and all(m in ('DA','DC','DG','DT') for m in monomers):kind='dna';seq=''.join(m[1:] for m in monomers)
        elif len(monomers)>=2 and all(m in ('A','C','G','U') for m in monomers):kind='rna';seq=''.join(monomers)
        else:kind='ligand';seq=None;require(len(residues)==1,'Multiple unknown/nonstandard residues require explicit adapter')
        if kind in ('protein','dna','rna'):
            for position,residue in enumerate(residues.values()):
                names={a['atom_name'] for a in residue};required={'N','CA','C','O'} if kind=='protein' else {"C1'","C2'","C3'","C4'","C5'","O3'","O4'","O5'"}
                monomer=residue[0]['monomer']
                if kind=='protein':required.update(SIDECHAIN[monomer].split())
                else:
                    required.update(BASE_ATOMS[monomer[-1]].split())
                    if position>0:required.update(('P','OP1','OP2'))
                    if kind=='rna':required.add("O2'")
                require(required<=names,'Missing standard polymer heavy atom')
        observed_heavy={a['atom_name']:a['element'] for a in group if a['element'] not in ('H','D')}
        counts=dict(collections.Counter(a['element'] for a in group if a['element'] not in ('H','D')))
        matches=[];details=[]
        for entity in entities:
            if kind!=entity['type']:continue
            if kind!='ligand':ok=seq==entity['sequence'];basis='Exact polymer type and complete residue sequence'
            elif 'ccd_codes' in entity:
                ok=monomers==[entity['ccd_codes']];reference=entity.get('chemical_reference')
                if reference:ok=ok and observed_heavy==reference['heavy_atoms']
                basis='CCD identity and exact pinned heavy-atom names/elements' if reference else 'CCD residue label only; atom graph/CCD snapshot unverified'
            else:
                ok=counts==entity['heavy_element_counts'];basis='Heavy-element composition only; native graph and stereo mapping unresolved'
            if ok:
                matches.extend((entity['entity_id'],copy) for copy in entity['chain_ids']);details.append({'entity_id':entity['entity_id'],'basis':basis})
        require(matches,'Native chain does not match any requested entity sequence/chemical composition')
        candidates[chain]=matches;descriptions[chain]={'output_chain':chain,'auth_chains':sorted({a['auth_chain'] for a in group}),'type':kind,'sequence':seq,'residue_count':len(residues),'heavy_element_counts':counts,'candidate_input_copies':[{'entity_id':eid,'chain_id':cid} for eid,cid in matches],'identity_evidence':details,'atoms':group}
    # Verify a complete one-to-one assignment exists without claiming one arbitrary
    # assignment resolves indistinguishable copies or isomers.
    assigned={}
    def augment(chain,seen):
        for candidate in candidates[chain]:
            if candidate in seen:continue
            seen.add(candidate)
            if candidate not in assigned or augment(assigned[candidate],seen):assigned[candidate]=chain;return True
        return False
    require(all(augment(chain,set()) for chain in chains),'No complete input/output entity assignment')
    for chain,record in descriptions.items():
        record['copy_mapping_status']='unique_by_identity' if len(candidates[chain])==1 else 'ambiguous_equivalent_candidates'
        record['chemical_graph_qualified']=False if record['type']=='ligand' else None
        record['mapping_limitations']='Candidate sets preserve ambiguity; matching element counts alone does not establish chemical graph or stereochemistry.'
    return list(descriptions.values())

def parse_response(raw,request,identity):
    data=strict_json(raw);require(isinstance(data,dict) and set(data)=={'request_id','outputs'},'Native response schema mismatch')
    require(data['request_id']==request['request_id'] and isinstance(data['outputs'],list) and len(data['outputs'])==1,'Native request echo/output count mismatch')
    out=data['outputs'][0];require(isinstance(out,dict) and {'input_id','structures_with_scores'}<=set(out) and set(out)<={'input_id','structures_with_scores','runtime_metrics'},'Native output schema mismatch')
    original=request['inputs'][0];require(out['input_id']==original['input_id'],'Native input ID mismatch')
    samples=out['structures_with_scores'];require(isinstance(samples,list) and len(samples)==original['diffusion_samples'],'Native sample count mismatch')
    structures=[];scores=[];ranks=[]
    fields=('confidence_score','complex_plddt_score','complex_pde_score','ptm_score','iptm_score')
    for index,sample in enumerate(samples,1):
        require(isinstance(sample,dict) and {'structure','format',*fields}<=set(sample) and set(sample)<= {'structure','format','name','source',*fields},'Native sample schema mismatch')
        require(sample['format']==original['output_format'],'Native output format mismatch')
        atoms,qa=read_atoms(sample['structure'],sample['format']);mapping=map_structure(atoms,identity['entities'])
        values={key:sample[key] for key in fields};require(all(v is None or finite(v) for v in values.values()),'Nonfinite native confidence')
        require(finite(values['confidence_score']),'Ranking confidence unavailable')
        require(values['complex_plddt_score'] is None or 0<=values['complex_plddt_score']<=100,'pLDDT outside documented range')
        require(values['complex_pde_score'] is None or values['complex_pde_score']>=0,'Negative PDE score')
        for key in ('ptm_score','iptm_score'):require(values[key] is None or 0<=values[key]<=1,'TM confidence outside range')
        ranks.append(values['confidence_score']);artifact=sample['structure'].encode()
        record={'rank_1based':index,'name':sample.get('name'),'native_source':sample.get('source'),'format':sample['format'],'text':sample['structure'],'artifact':{'sha256':sha(artifact),'bytes':len(artifact),'coordinate_unit':'angstrom'},'chain_mappings':mapping}
        structures.append(record);scores.append({'rank_1based':index,'name':sample.get('name'),'values':values,'available':{key:value is not None for key,value in values.items()},'native_modelcif_qa':qa})
    require(ranks==sorted(ranks,reverse=True),'Native rank order contradicts confidence')
    metrics=out.get('runtime_metrics');require(metrics is None or isinstance(metrics,dict),'Invalid runtime metrics')
    return {'predicted_complex':{'entities':identity['entities'],'structures':structures},'confidence_scores':{'records':scores,'semantics':{'confidence_score':'Native sample ranking; formula/scale unresolved. Not measured accuracy.','complex_plddt_score':'Native average pLDDT score, 0–100.','complex_pde_score':'Native predicted distance-error score; exact aggregation/units unresolved. Not PAE.','ptm_score':'Dimensionless predicted TM score, 0–1.','iptm_score':'Dimensionless predicted interface TM score, 0–1; null means unavailable.','native_modelcif_qa':'Preserved named native metric declarations and values; no automatic B-factor interpretation.','qualification':'No affinity, biological function or scientific accuracy established by confidence.'}},'inference_metrics':{'native':metrics,'semantics':{'units':'Native fields require version-specific meanings; empty/null metrics stay unavailable.','temporal':'Finite prediction; diffusion sample is not biological time.'}},'native_response':raw.decode()}

class NvidiaOpenFold3(BioModule):
    execution_policy=ExecutionPolicy.ONCE_BEFORE_RUN
    def __init__(self,observation_timeout_s=600.0,max_response_bytes=20000000,resume_request_id=None,resume_request_sha256=None):
        require(finite(observation_timeout_s) and 1<=observation_timeout_s<=1200,'Invalid observation limit')
        require(type(max_response_bytes) is int and 1024<=max_response_bytes<=50000000,'Invalid response bound')
        require((resume_request_id is None)==(resume_request_sha256 is None),'Resume requires UUID and input hash')
        if resume_request_sha256 is not None:require(isinstance(resume_request_sha256,str) and re.fullmatch('[0-9a-f]{64}',resume_request_sha256),'Invalid resume digest')
        self.timeout=float(observation_timeout_s);self.max_bytes=max_response_bytes;self.resume_id=str(uuid.UUID(resume_request_id)) if resume_request_id else None;self.resume_sha=resume_request_sha256;self.state={'status':'not_started'}
    def inputs(self):
        return {'molecular_entities':SignalSpec.record(schema=ENTITIES,accepted_units=('1',),required=True),'alignments':SignalSpec.record(schema=ENTRIES,accepted_units=('1',),required=False,default={'entries':[]}),'templates':SignalSpec.record(schema=ENTRIES,accepted_units=('1',),required=False,default={'entries':[]}),'run_options':SignalSpec.record(schema=OPTIONS,accepted_units=('1',),required=False,default={'parameters':{}}),'query_metadata':SignalSpec.record(schema=META,accepted_units=('1',),required=False,default={'query_id':'','request_id':'','provenance':{}})}
    def outputs(self):
        return {'predicted_complex':SignalSpec.record(schema=STRUCTURES,emitted_unit='1'),'confidence_scores':SignalSpec.record(schema=REPORT,emitted_unit='1'),'inference_metrics':SignalSpec.record(schema=METRICS,emitted_unit='1'),'run_provenance':SignalSpec.record(schema=PROVENANCE,emitted_unit='1'),'native_response':SignalSpec.scalar(dtype='str',emitted_unit='1')}
    def snapshot(self):return dict(self.state)
    def execute(self,inputs,*,context):
        values={k:v.value for k,v in inputs.items()};require('molecular_entities' in values,'Missing entities')
        values.setdefault('alignments',{'entries':[]});values.setdefault('templates',{'entries':[]});values.setdefault('run_options',{'parameters':{}});values.setdefault('query_metadata',{'query_id':'','request_id':'','provenance':{}})
        request,identity=prepare_request(values['molecular_entities'],values['alignments'],values['templates'],values['run_options'],values['query_metadata'])
        body=encode(request);digest=sha(body);require(not self.resume_id or self.resume_sha==digest,'Resume request digest mismatch')
        key=os.environ.get('NVIDIA_API_KEY');require(isinstance(key,str) and bool(key.strip()),'Secure NVIDIA_API_KEY not provisioned on worker')
        started=time.monotonic();deadline=started+self.timeout;self.state={'status':'observing','request_sha256':digest,'provider_request_id':self.resume_id}
        method,url=('GET',POLL+self.resume_id) if self.resume_id else ('POST',ENDPOINT)
        try:
            with httpx.Client(follow_redirects=False,trust_env=False) as client:
                while True:
                    remaining=deadline-time.monotonic();require(remaining>0,'Observation deadline reached')
                    headers={'Authorization':'Bearer '+key,'Content-Type':'application/json','NVCF-POLL-SECONDS':str(min(30,max(1,int(remaining)-1)))}
                    with client.stream(method,url,headers=headers,content=body if method=='POST' else None,timeout=min(60,remaining)) as response:
                        rid=response.headers.get('nvcf-reqid') or response.headers.get('nvcf-request-id')
                        if rid:
                            rid=str(uuid.UUID(rid));require(self.state['provider_request_id'] in (None,rid),'Provider UUID changed');self.state['provider_request_id']=rid
                        safe={k:v for k,v in response.headers.items() if k.lower() in {'date','content-type','nvcf-reqid','nvcf-request-id','retry-after'}}
                        raw=b''
                        for chunk in response.iter_bytes():require(len(raw)+len(chunk)<=self.max_bytes and time.monotonic()<deadline,'Response observation bound exceeded');raw+=chunk
                        status=response.status_code
                    self.state.update({'last_http_status':status,'response_sha256':sha(raw),'response_bytes':len(raw),'safe_headers':safe})
                    if status==202:
                        require(bool(self.state['provider_request_id']),'Accepted request missing UUID');method,url='GET',POLL+self.state['provider_request_id'];time.sleep(min(1,max(0,deadline-time.monotonic())));continue
                    require(status==200,'NVIDIA HTTP failure');require(key.encode() not in raw,'Unexpected secret in output');break
            self.state['phase']='output_contract';outputs=parse_response(raw,request,identity)
        except Exception:
            self.state['status']='observation_or_contract_failed'
            raise RuntimeError('OpenFold3 observation/contract failed. Inspect safe snapshot; preserve existing UUID/hash and do not automatically resubmit.') from None
        self.state['status']='completed'
        outputs['run_provenance']={'request':{'body':request,'sha256':digest},'provider':{'endpoint':ENDPOINT,'request_id':self.state['provider_request_id'],'response_sha256':sha(raw),'response_bytes':len(raw),'safe_headers':safe,'observed_at':datetime.now(timezone.utc).isoformat(),'wall_seconds':time.monotonic()-started,'resumed':bool(self.resume_id)},'inputs':{'values':values,'checked_identity':identity},'backend_identity':{'documented_hosted_schema':'1.0.0 (stale AlphaFold2 title)','self_hosted_target':'NIM 1.6.0','deployed_version':None,'checkpoint_digest':None,'container_digest':None,'ccd_snapshot_digest':None,'provider_gpu':None,'server_seed':None},'limitations':['Exact release identity, artifact rights and managed/scientific qualification remain unresolved.','Chemical graph/stereochemistry and ambiguous equivalent copy maps are not silently resolved.','HTTP success and model confidence are not scientific accuracy; all requested sample artifacts retained.','Template selection/use is internal and not proven by request acceptance; paired species identity remains provider-asserted.']}
        return outputs
