"""Independently enumerate all nonzero words of the displayed J=4,L=5 example."""
from pathlib import Path
import json,subprocess,time
R=Path(__file__).resolve().parents[1];c=next(c for c in json.loads((R/'joint_examples.json').read_text()) if c['J']==4 and c['L']==5)
J,L,P,E=c['J'],c['L'],c['P'],c['E'];H=[0]*(J*P)
for i,row in enumerate(E):
 for l,e in enumerate(row):
  for t in range(P):H[i*P+(t+e)%P]|=1<<(t*L+l)
piv={}
for h in H:
 while h:
  q=h.bit_length()-1
  if q in piv:h^=piv[q]
  else:piv[q]=h;break
basis=[]
for j in range(L*P):
 if j in piv:continue
 v=1<<j
 for q in sorted(piv):
  if (v&piv[q]).bit_count()%2:v|=1<<q
 assert all((v&h).bit_count()%2==0 for h in H);basis.append(v)
assert L*P<=128 and len(basis)<=28 and len(basis)==c['k']
folder=R/'verification/recomputed/J4_L5_enumeration';folder.mkdir(parents=True,exist_ok=True)
f=folder/'basis.txt';f.write_text(f'{L*P} {len(basis)}\n'+'\n'.join(f'{v&((1<<64)-1)} {v>>64}' for v in basis))
exe=folder/'enumerate128';subprocess.run(['c++','-O3','-std=c++17',str(R/'verification/enumerate128.cpp'),'-o',str(exe)],check=True)
start=time.time();d=json.loads(subprocess.check_output([str(exe),str(f)],text=True));v=sum(1<<q for q in d['support']);assert all((v&h).bit_count()%2==0 for h in H)
assert d['complete'] and d['d']==c['d']==30 and d['nonzero_words_checked']==2**c['k']-1 and d['count']==c['minimum_weight_word_count']
d.update(E=E,J=J,L=L,P=P,elapsed_sec=time.time()-start,method='Independent GF(2) kernel basis followed by exhaustive Gray-code enumeration')
(R/'verification/j4_l5_enumeration.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d))
