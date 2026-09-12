"""Check archived complete search coverage and codeword witnesses, not rerun DFS."""
import collections,hashlib,json,zipfile,sys
from pathlib import Path

def check(path):
 with zipfile.ZipFile(path) as z:
  cs=json.loads(z.read('cases.json'))
  for c in cs:
   J,L,P,E=c['J'],c['L'],c['P'],c['E'];assert J==4 and len(E)==J and all(len(r)==L for r in E)
   H=[0]*(J*P)
   for i,row in enumerate(E):
    for l,e in enumerate(row):
     for t in range(P):H[i*P+(t+e)%P]^=1<<(t*L+l)
   assert all(h.bit_count()==L for h in H)
   piv={}
   for h in H:
    while h:
     q=h.bit_length()-1
     if q in piv:h^=piv[q]
     else:piv[q]=h;break
   assert c['n']==L*P and c['k']==L*P-len(piv)
   v=c['witness'];assert len(v)==len(set(v))==c['d_upper'] and v and all(0<=q<L*P for q in v)
   w=sum(1<<q for q in v);assert all((w&h).bit_count()%2==0 for h in H)
   if L==5:
    e=json.loads(z.read(c['enumeration_file']));assert e['complete'] and e['nonzero_words_checked']==2**c['k']-1 and e['d']==c['d_lower']==c['d_upper']==30 and e['E']==E and e['count']==23
   else:
    seen=collections.defaultdict(set);cfg={}
    for a in c['certificates']:
     b=z.read(a['file']);assert hashlib.sha256(b).hexdigest()==a['sha256'];d=json.loads(b)
     assert d['complete'] and not d['found'] and d['stats']['kernel_circuits']==0
     assert (d['mode'],d['side'],d['stabilizer_rows'],d['instance'],d['n'],d['P'],d['base_columns'],d['check_rows'],d['max_weight'],d['effective_max_weight'])==('classical','X',0,c['name'],L*P,P,L,J*P,c['d_lower']-2,c['d_lower']-2)
     assert d['kernel_even_weight'] and d['qc_root_reduction'] and d['qc_root_partition']=='earliest_base_column_type_v2'
     r=d['root'];sh=d['shard_index'];assert 0<=r<L and sh not in seen[r];seen[r].add(sh)
     pair=(d['shard_count'],d['split_weight']);assert cfg.setdefault(r,pair)==pair
    assert set(seen)==set(range(L)) and all(seen[r]==set(range(cfg[r][0])) for r in range(L))
   assert c['d']==(c['d_lower'] if c['d_lower']==c['d_upper'] else None)
   lines=z.read(c['matrix_file']).decode().splitlines();assert lines[0]=='CPM_CSS_V1' and lines[1].split()==[c['name'],str(P),str(L),'4','0'] and [list(map(int,t.split())) for t in lines[2:]]==E
  manifest=json.loads(z.read('source_manifest.json'))
  for p,sha in manifest.items():assert hashlib.sha256(z.read(p)).hexdigest()==sha
  if 'frontier_parts_verification.zip' in z.namelist():check_frontier_archive(z.read('frontier_parts_verification.zip'),cs,z)
  return cs
def check_frontier_archive(data, cases, main_zip):
 import io,struct
 def parse(b):
  assert b[:8]==b'CPMCHK2\0' and struct.unpack_from('<I',b,8)[0]==3
  nw,mw=struct.unpack_from('<ii',b,40);off=149+8*nw;count=struct.unpack_from('<Q',b,off)[0];stride=5+8*(2*nw+mw)
  assert len(b)==off+8+count*stride
  return off,[b[off+8+i*stride:off+8+(i+1)*stride] for i in range(count)],struct.unpack_from('<8Q',b,80),b[144]
 c=next(c for c in cases if c['L']==7);assert c['d_lower']==32
 main={}
 for a in c['certificates']:
  d=json.loads(main_zip.read(a['file']));main[(d['root'],d['shard_index'])]=a['sha256']
 with zipfile.ZipFile(io.BytesIO(data)) as z:
  report=json.loads(z.read('manifest.json'));assert report['complete'] and not report['remaining_parents']
  ids=[p['parent_id'] for p in report['parents']]+report['already_complete_parent_ids'];assert len(ids)==len(set(ids))==7*2048 and set(ids)==set(range(7*2048))
  fragments=set()
  for p in report['parents']:
   u=p['unit'];assert u['case']==c['name'] and u['cutoff']==30 and u['shards']==2048 and p['parent_id']==u['id']==u['root']*2048+u['shard']
   b=z.read(p['file']);assert hashlib.sha256(b).hexdigest()==p['sha256'];off,nodes,stats,found=parse(b);assert not found and stats[6]==0 and nodes and all(n[4]==1 for n in nodes)
   assert struct.unpack_from('<7i',b,20)==(u['root'],30,5,2048,u['shard'],5,3)
   k=len(p['fragments']);assert k==len(p['finished']) and k>0;tot=[0]*8
   for j,(f,g) in enumerate(zip(p['fragments'],p['finished'])):
    assert f['id']==g['id'] and f['id'] not in fragments;fragments.add(f['id'])
    a=z.read(f['file']);assert hashlib.sha256(a).hexdigest()==f['sha256'];head=bytearray(b[:off]);ns=nodes[j::k]
    if j:struct.pack_into('<8Q',head,80,*([0]*8));struct.pack_into('<d',head,72,0.0)
    assert a==bytes(head)+struct.pack('<Q',len(ns))+b''.join(ns)
    done=z.read(g['checkpoint']);_,remaining,ds,found=parse(done);assert done[:56]==b[:56] and not remaining and not found and ds[6]==0
    assert all(x>=y for x,y in zip(ds,parse(a)[2]));d=json.loads(z.read(g['result']));assert d['complete'] and not d['found'] and tuple(d['stats'].values())==ds
    assert (d['mode'],d['side'],d['stabilizer_rows'],d['instance'],d['root'],d['shard_index'],d['shard_count'],d['max_weight'],d['effective_max_weight'],d['n'],d['P'],d['base_columns'],d['check_rows'])==('classical','X',0,c['name'],u['root'],u['shard'],2048,30,30,301,43,7,172)
    tot=[x+y for x,y in zip(tot,ds)]
   merged=z.read(p['merged_checkpoint']);_,remaining,ds,found=parse(merged);assert merged[:56]==b[:56] and not remaining and not found and tuple(tot)==ds
   result=z.read(p['merged_result']);d=json.loads(result);assert d['complete'] and not d['found'] and tuple(d['stats'].values())==ds
   assert hashlib.sha256(result).hexdigest()==p['merged_result_sha256']==main[(u['root'],u['shard'])]

if __name__=='__main__':
 path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]/'supporting_material/J4_distance_refinement_verification.zip'
 cs=check(path);print(json.dumps([{k:c[k] for k in ['J','L','P','n','k','d_lower','d_upper','d']} for c in cs]))
