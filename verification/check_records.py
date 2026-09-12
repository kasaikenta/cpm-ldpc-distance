"""Verify archived run records, coverage, matrices, witnesses, and table values.
This does NOT independently reprove nonexistence: rerun the solver for that.
"""
import hashlib, json, zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def archive(j):
 name='J3_L5to8_distance24_certificates.zip' if j==3 else 'J4_L5to8_distance26_verification.zip'
 path=ROOT/'supporting_material'/name
 if not path.is_file():path=ROOT.parent/('descending_lift' if j==3 else 'distance26_j4_smaller')/'results'/name
 return zipfile.ZipFile(path)
def main():
 table=json.loads((ROOT/'joint_examples.json').read_text());summary=[]
 for j in (3,4):
  with archive(j) as z:
   records=json.loads(z.read('results/best_verified.json'))
   cases={c['name']:c for c in json.loads(z.read('data/all_dfs_cases.json'))}
   assert z.read('src/distance.cpp')==(ROOT/'verification/distance.cpp').read_bytes()
   for r in records:
    P,L,E=r['P'],r['L'],r['E'];seen={i:set() for i in range(L)};configs={}
    tc=next(c for c in table if (c['J'],c['L'])==(j,L))
    for key in ['J','L','P','E','n','k']:assert tc[key]==r[key],key
    if j==3:assert tc['d']==r['d']==24
    assert cases[r['name']]['E']==E and cases[r['name']]['P']==P
    for a in r['certificates']:
     b=z.read(a['file']);assert hashlib.sha256(b).hexdigest()==a['sha256']
     d=json.loads(b);assert d['complete'] and not d['found']
     assert (d['mode'],d['side'],d['stabilizer_rows'])==('classical','X',0)
     assert (d['instance'],d['P'],d['base_columns'],d['check_rows'],d['max_weight'])==(r['name'],P,L,j*P,22 if j==3 else 24)
     assert d['kernel_even_weight'] and d['qc_root_reduction']
     assert d['qc_root_partition']=='earliest_base_column_type_v2'
     root=d['root'];assert root in seen and d['shard_index'] not in seen[root]
     seen[root].add(d['shard_index'])
     config=(d['shard_count'],d['split_weight']);assert configs.setdefault(root,config)==config
    assert all(seen[root]==set(range(configs[root][0])) for root in range(L))
    H=[0]*(j*P)
    for l in range(L):
     for t in range(P):
      for i in range(j):H[i*P+(t+E[i][l])%P]^=1<<(t*L+l)
    assert all(h.bit_count()==L for h in H)
    v=r['weight24_witness'] if j==3 else r.get('minimum_witness',r.get('upper_bound_witness'))
    expected=24 if j==3 else r['d_upper']
    assert len(v)==len(set(v))==expected and all(0<=a<P*L for a in v)
    if j==4 and L==5:
     enum=json.loads(z.read(r['enumeration_file']))
     assert enum['complete'] and enum['d']==30 and enum['nonzero_words_checked']==2**r['k']-1
     assert enum['support']==v and enum['E']==E and enum['count']==r['minimum_weight_word_count']
    w=sum(1<<a for a in v);assert all((h&w).bit_count()%2==0 for h in H)
    piv={}
    for h in H:
     while h:
      k=h.bit_length()-1
      if k in piv:h^=piv[k]
      else:piv[k]=h;break
    assert len(piv)==r['rank'] and P*L-len(piv)==r['k']
    summary.append({'J':j,'L':L,'P':P,'k':r['k'],'record_count':sum(map(len,seen.values())),'witness_weight':expected,'witness_checks_passed':True,'d_lower':24 if j==3 else r['d_lower'],'exact_distance':r['d']})
 from check_refined_j4 import check
 refined=check(ROOT/'supporting_material/J4_distance_refinement_verification.zip')
 for c in refined:
  tc=next(t for t in table if (t['J'],t['L'])==(4,c['L']))
  for k in ['E','P','n','k','d','d_lower','d_upper']:assert tc[k]==c[k]
  old=next(s for s in summary if (s['J'],s['L'])==(4,c['L']))
  old.update(d_lower=c['d_lower'],d_upper=c['d_upper'],exact_distance=c['d'],witness_weight=c['d_upper'],refined_record_count=len(c.get('certificates',[])))
 (ROOT/'verification/record_checks.json').write_text(json.dumps(summary,indent=2)+'\n')
 print(json.dumps(summary))
if __name__=='__main__':main()
