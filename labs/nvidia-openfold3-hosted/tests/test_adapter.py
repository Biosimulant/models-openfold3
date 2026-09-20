"""Synthetic interface/transport tests; real captures are audited separately."""
import copy,importlib.util,io,json
from pathlib import Path
import httpx,pytest,yaml
from Bio.PDB import MMCIFIO
from biosim import BioModule,BioWorld,ExecutionContext,ExecutionPolicy,SignalSpec,make_signal
from biosim.pack import validate_lab_source
LAB=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('openfold3_model',LAB/'owned/models/main/model.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
CTX=ExecutionContext(policy=ExecutionPolicy.ONCE_BEFORE_RUN,run_start=0,run_end=1)
RID='12345678-1234-5678-1234-567812345678'
def pdb(chain='A'):
    lines=[];serial=0
    for number,resname in [(1,'ALA'),(2,'CYS')]:
        for offset,name in enumerate(['N','CA','C','O','CB']+(['SG'] if resname=='CYS' else [])):
            serial+=1;lines.append(f'ATOM  {serial:5d} {name:^4} {resname} {chain}{number:4d}    {float(serial):8.3f}{float(offset):8.3f}{0.0:8.3f}  1.00 75.00')
    return '\n'.join(lines)+'\n'
@pytest.fixture
def values():
    text='>query\nAC\n>hit\nAc-\n'
    return {'molecular_entities':{'entities':[{'entity_id':'protein-1','type':'protein','chain_ids':['A'],'sequence':'AC','provenance':{'fixture':'synthetic'}}]},'alignments':{'entries':[{'entity_id':'protein-1','kind':'unpaired','database':'test-db','format':'a3m','text':text,'sha256':m.sha(text.encode()),'provenance':{}}]},'templates':{'entries':[]},'run_options':{'parameters':{'output_format':'pdb','diffusion_samples':2}},'query_metadata':{'query_id':'test-query','request_id':'test-request','provenance':{}}}
def prepared(v):return m.prepare_request(v['molecular_entities'],v['alignments'],v['templates'],v['run_options'],v['query_metadata'])
def response():return {'request_id':'test-request','outputs':[{'input_id':'test-query','runtime_metrics':{},'structures_with_scores':[{'structure':pdb(),'format':'pdb','name':'s1','source':'synthetic','confidence_score':0.8,'complex_plddt_score':80,'complex_pde_score':1,'ptm_score':0.7,'iptm_score':None},{'structure':pdb(),'format':'pdb','name':'s2','source':'synthetic','confidence_score':0.7,'complex_plddt_score':75,'complex_pde_score':2,'ptm_score':0.6,'iptm_score':None}]}]}
def parse(v,d):request,identity=prepared(v);return m.parse_response(m.encode(d),request,identity)
def signals(v):return {k:make_signal(spec=m.NvidiaOpenFold3().inputs()[k],value=x,source='fixture',name=k,emitted_at=0) for k,x in v.items()}
def mock(monkeypatch,handler):
    real=httpx.Client;monkeypatch.setattr(m.httpx,'Client',lambda **kw:real(transport=httpx.MockTransport(handler),**kw));monkeypatch.setenv('NVIDIA_API_KEY','unit-test-key');monkeypatch.setattr(m.time,'sleep',lambda _:None)
def test_manifest_runtime_lab_parity():
    model=m.NvidiaOpenFold3();manifest=yaml.safe_load((LAB/'owned/models/main/model.yaml').read_text());lab=yaml.safe_load((LAB/'lab.yaml').read_text())
    for direction,ports in [('inputs',model.inputs()),('outputs',model.outputs())]:
        entries={p['name']:p for p in manifest['io'][direction]};assert set(entries)==set(ports)
        assert {p['maps_to'] for p in lab['io'][direction]}=={'main.'+k for k in ports}
        for name,port in ports.items():
            for key in ['signal_type','schema','dtype','required','default','emitted_unit']:assert entries[name].get(key)==getattr(port,key)
            assert tuple(entries[name].get('accepted_units',()))==(port.accepted_units or ())
    result=validate_lab_source(LAB);assert result.valid,result.errors
    assert len(lab['models'])==1 and model.execution_policy==ExecutionPolicy.ONCE_BEFORE_RUN
@pytest.mark.parametrize('sequence',['','A','AX','ac','A C','AC*','A'*4097,None])
def test_invalid_protein(values,sequence):
    values['molecular_entities']['entities'][0]['sequence']=sequence
    with pytest.raises(ValueError):prepared(values)
@pytest.mark.parametrize('mutation',['duplicate','bad-chain','unknown-field','no-msa','bad-entity','empty','pdb-long-chain','dna-wrong-base','rna-missing-msa'])
def test_entity_failures(values,mutation):
    entity=values['molecular_entities']['entities'][0]
    if mutation=='duplicate':entity['chain_ids']=['A','A']
    elif mutation=='bad-chain':entity['chain_ids']=['!']
    elif mutation=='unknown-field':entity['pocket']={}
    elif mutation=='no-msa':values['alignments']['entries']=[]
    elif mutation=='bad-entity':entity['entity_id']=''
    elif mutation=='empty':values['molecular_entities']['entities']=[]
    elif mutation=='pdb-long-chain':entity['chain_ids']=['AB']
    elif mutation=='dna-wrong-base':entity.update(type='dna',sequence='AU');values['alignments']['entries']=[]
    elif mutation=='rna-missing-msa':entity.update(type='rna',sequence='AU');values['alignments']['entries']=[]
    with pytest.raises(ValueError):prepared(values)
@pytest.mark.parametrize('key,value',[('diffusion_samples',0),('diffusion_samples',6),('diffusion_samples',True),('output_format','mmcif'),('seed',17)])
def test_bad_options(values,key,value):
    values['run_options']['parameters'][key]=value
    with pytest.raises(ValueError):prepared(values)
@pytest.mark.parametrize('mutation',['hash','query','format','columns','duplicate','entity','paired-single','annotation'])
def test_alignment_failures(values,mutation):
    entry=values['alignments']['entries'][0]
    if mutation=='hash':entry['sha256']='0'*64
    elif mutation=='query':entry['text']='>q\nAA\n'
    elif mutation=='format':entry['format']='fasta'
    elif mutation=='columns':entry['text']='>q\nAC\n>h\nAAA\n'
    elif mutation=='duplicate':values['alignments']['entries'].append(copy.deepcopy(entry))
    elif mutation=='entity':entry['entity_id']='missing'
    elif mutation=='paired-single':entry['kind']='paired'
    elif mutation=='annotation':entry['text']='>ss_pred\nAC\n'
    with pytest.raises(ValueError):prepared(values)
def test_csv_query_keys_and_original_text(values):
    e=values['alignments']['entries'][0];e.update(format='csv',text='key,sequence\n-1,AC\n10,A-\n');e['sha256']=m.sha(e['text'].encode());r,i=prepared(values)
    assert r['inputs'][0]['molecules'][0]['msa']['test-db']['csv']['alignment']==e['text']
    assert i['alignments'][0]['row_keys']==['-1','10']
def test_native_header_tabs_preserved(values):
    e=values['alignments']['entries'][0];e['text']='>query\tmetadata\nAC\n>hit\t1\nA-\n';e['sha256']=m.sha(e['text'].encode())
    r,i=prepared(values);assert i['alignments'][0]['row_keys']==['query\tmetadata','hit\t1']
def test_dna_copies_and_no_implicit_strand(values):
    values['molecular_entities']['entities'][0].update(type='dna',sequence='ACGT',chain_ids=['A','B']);values['alignments']['entries']=[]
    r,i=prepared(values);assert r['inputs'][0]['molecules'][0]['id']==['A','B'] and len(i['entities'])==1
def test_rna_query_only(values):
    values['molecular_entities']['entities'][0].update(type='rna',sequence='ACGU');e=values['alignments']['entries'][0];e['text']='>q\nACGU\n';e['sha256']=m.sha(e['text'].encode());r,i=prepared(values);assert i['alignments'][0]['query_only']
def test_smiles_input_and_chemical_exclusivity(values):
    values['molecular_entities']['entities']=[{'entity_id':'lig','type':'ligand','chain_ids':['L'],'smiles':'CCO','provenance':{}}];values['alignments']['entries']=[]
    r,i=prepared(values);assert i['entities'][0]['heavy_element_counts']=={'C':2,'O':1}
    values['molecular_entities']['entities'][0]['ccd_codes']='ATP'
    with pytest.raises(ValueError):prepared(values)
def test_template_hash_and_chain(values):
    parser=m.PDBParser(QUIET=True);structure=parser.get_structure('x',io.StringIO(pdb()));writer=MMCIFIO();writer.set_structure(structure);buf=io.StringIO();writer.save(buf);text=buf.getvalue()
    values['templates']['entries']=[{'entity_id':'protein-1','text':text,'sha256':m.sha(text.encode()),'source':'synthetic','name':'test','chain_id':'A','provenance':{}}]
    r,i=prepared(values);assert r['inputs'][0]['molecules'][0]['structural_templates'][0]['format']=='cif'
    values['templates']['entries'][0]['chain_id']='Z'
    with pytest.raises(ValueError):prepared(values)
def test_all_samples_and_null_confidence(values):
    out=parse(values,response());assert len(out['predicted_complex']['structures'])==2
    assert out['confidence_scores']['records'][0]['values']['iptm_score'] is None
    assert out['confidence_scores']['records'][0]['available']['iptm_score'] is False
    assert out['inference_metrics']['native']=={}
def test_renamed_unique_polymer_mapping(values):
    data=response()
    for s in data['outputs'][0]['structures_with_scores']:s['structure']=pdb('B')
    out=parse(values,data);mapping=out['predicted_complex']['structures'][0]['chain_mappings'][0]
    assert mapping['output_chain']=='B' and mapping['candidate_input_copies']==[{'entity_id':'protein-1','chain_id':'A'}]
def test_identical_copy_ambiguity_not_invented(values):
    values['molecular_entities']['entities'][0]['chain_ids']=['A','B'];data=response()
    for s in data['outputs'][0]['structures_with_scores']:s['structure']=pdb('C')+'TER\n'+pdb('D')
    out=parse(values,data);maps=out['predicted_complex']['structures'][0]['chain_mappings']
    assert all(x['copy_mapping_status']=='ambiguous_equivalent_candidates' and len(x['candidate_input_copies'])==2 for x in maps)
@pytest.mark.parametrize('mutation',['request','input','missing-sample','format','sequence','chain-count','backbone','rank','plddt','ptm','pde'])
def test_output_contract_failures(values,mutation):
    d=response();s=d['outputs'][0]['structures_with_scores'][0]
    if mutation=='request':d['request_id']='wrong'
    elif mutation=='input':d['outputs'][0]['input_id']='wrong'
    elif mutation=='missing-sample':d['outputs'][0]['structures_with_scores'].pop()
    elif mutation=='format':s['format']='cif'
    elif mutation=='sequence':s['structure']=s['structure'].replace('CYS','GLY')
    elif mutation=='chain-count':s['structure']+=pdb('B')
    elif mutation=='backbone':s['structure']='\n'.join(x for x in s['structure'].splitlines() if x[12:16].strip()!='CA')+'\n'
    elif mutation=='rank':s['confidence_score']=0.1
    elif mutation=='plddt':s['complex_plddt_score']=101
    elif mutation=='ptm':s['ptm_score']=-1
    elif mutation=='pde':s['complex_pde_score']=-1
    with pytest.raises((ValueError,KeyError)):parse(values,d)
def test_chemical_composition_loss_rejected():
    atoms=[{'chain':'A','auth_chain':'A','residue_id':'1','auth_residue_id':'1','insertion_code':'','monomer':'LIG','atom_name':'C1','element':'C','xyz':[0,0,0],'model_id':'1'}]
    expected=[{'entity_id':'lig','type':'ligand','chain_ids':['L'],'smiles':'CO','heavy_element_counts':{'C':1,'O':1}}]
    with pytest.raises(ValueError,match='chemical composition'):m.map_structure(atoms,expected)
def test_missing_sidechain_atom_rejected(values):
    d=response();s=d['outputs'][0]['structures_with_scores'][0];s['structure']='\n'.join(x for x in s['structure'].splitlines() if x[12:16].strip()!='SG')+'\n'
    with pytest.raises(ValueError,match='heavy atom'):parse(values,d)
def test_single_ccd_nucleotide_is_not_invented_polymer():
    atoms=[{'chain':'A','auth_chain':'A','residue_id':'1','auth_residue_id':'1','insertion_code':'','monomer':'A','atom_name':'C1','element':'C','xyz':[0,0,0],'model_id':'1'}]
    mapped=m.map_structure(atoms,[{'entity_id':'lig','type':'ligand','chain_ids':['L'],'ccd_codes':'A','chemical_reference':None}])
    assert mapped[0]['type']=='ligand' and mapped[0]['chemical_graph_qualified'] is False
@pytest.mark.parametrize('status',[401,413,422,429,500,503,504])
def test_http_failure_not_retried(values,monkeypatch,status):
    calls=[]
    def handler(req):calls.append(req.method);return httpx.Response(status,content=b'error',headers={'nvcf-reqid':RID})
    mock(monkeypatch,handler);model=m.NvidiaOpenFold3()
    with pytest.raises(RuntimeError):model.execute(signals(values),context=CTX)
    assert calls==['POST'] and model.snapshot()['provider_request_id']==RID
def test_pending_polls_existing_handle(values,monkeypatch):
    calls=[]
    def handler(req):
        calls.append((req.method,str(req.url)))
        return httpx.Response(202,content=b'pending',headers={'nvcf-reqid':RID}) if len(calls)==1 else httpx.Response(200,json=response(),headers={'nvcf-reqid':RID})
    mock(monkeypatch,handler);out=m.NvidiaOpenFold3().execute(signals(values),context=CTX)
    assert [x[0] for x in calls]==['POST','GET'] and calls[1][1]==m.POLL+RID and len(out)==5
def test_exact_resume(values,monkeypatch):
    r,_=prepared(values);calls=[]
    def handler(req):calls.append(req.method);return httpx.Response(200,json=response(),headers={'nvcf-reqid':RID})
    mock(monkeypatch,handler);model=m.NvidiaOpenFold3(resume_request_id=RID,resume_request_sha256=m.sha(m.encode(r)));model.execute(signals(values),context=CTX);assert calls==['GET']
    bad=m.NvidiaOpenFold3(resume_request_id=RID,resume_request_sha256='0'*64)
    with pytest.raises(ValueError):bad.execute(signals(values),context=CTX)
    assert calls==['GET']
def test_missing_credentials(values,monkeypatch):
    monkeypatch.delenv('NVIDIA_API_KEY',raising=False)
    with pytest.raises(ValueError,match='NVIDIA_API_KEY'):m.NvidiaOpenFold3().execute(signals(values),context=CTX)
def test_response_size_bound(values,monkeypatch):
    mock(monkeypatch,lambda req:httpx.Response(200,content=b'x'*2048,headers={'nvcf-reqid':RID}))
    with pytest.raises(RuntimeError):m.NvidiaOpenFold3(max_response_bytes=1024).execute(signals(values),context=CTX)
def test_bioworld_once_all_terminal_ports(values,monkeypatch):
    calls=[]
    def handler(req):calls.append(req.method);return httpx.Response(200,json=response(),headers={'nvcf-reqid':RID})
    mock(monkeypatch,handler)
    class Source(BioModule):
        execution_policy=ExecutionPolicy.ONCE_BEFORE_RUN
        def outputs(self):return {'molecular_entities':SignalSpec.record(schema=m.ENTITIES,emitted_unit='1'),'alignments':SignalSpec.record(schema=m.ENTRIES,emitted_unit='1'),'templates':SignalSpec.record(schema=m.ENTRIES,emitted_unit='1'),'run_options':SignalSpec.record(schema=m.OPTIONS,emitted_unit='1'),'query_metadata':SignalSpec.record(schema=m.META,emitted_unit='1')}
        def execute(self,inputs,*,context):return values
    model=m.NvidiaOpenFold3();world=BioWorld(communication_step=1);world.add_biomodule('source',Source());world.add_biomodule('main',model)
    for name in values:world.connect('source.'+name,'main.'+name)
    world.run(duration=1)
    assert calls==['POST'] and set(model.get_outputs())==set(model.outputs())
    assert all(signal.emitted_at==0 for signal in model.get_outputs().values())
    archived=json.loads(json.dumps({k:s.to_dict() for k,s in model.get_outputs().items()},allow_nan=False))
    assert set(archived)==set(model.outputs()) and all(x['emitted_at']==0 for x in archived.values())
