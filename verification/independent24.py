"""Independent exact distance for the P=24 example using two half-kernel lists."""
from pathlib import Path
import json,time
start=time.time();P=24;E=[[0,0,0,0],[0,1,2,4],[0,7,15,5]];rows=[]
for row in E:
 for t in range(P):rows.append(sum(1<<(j*P+(t-row[j])%P) for j in range(4)))
piv={}
for h in rows:
 while h:
  j=h.bit_length()-1
  if j in piv:h^=piv[j]
  else:piv[j]=h;break
basis=[]
for j in range(4*P):
 if j in piv:continue
 v=1<<j
 for q in sorted(piv):
  if (v&piv[q]).bit_count()%2:v|=1<<q
 assert all((v&h).bit_count()%2==0 for h in rows)
 basis.append(v)
assert len(basis)==26
left=[0];right=[0]
for v in basis[:13]:left += [x^v for x in left]
for v in basis[13:]:right += [x^v for x in right]
best=4*P;count=0
for a in left:
 for b in right:
  v=a^b
  if not v:continue
  w=v.bit_count()
  if w<best:best=w;count=1
  elif w==best:count+=1
assert best==24 and count==668
out=dict(P=P,E=E,n=96,k=len(basis),d=best,minimum_weight_word_count=count,nonzero_words_checked=len(left)*len(right)-1,elapsed_sec=time.time()-start,method='Independent Python implementation, high-bit echelon kernel, two half-list exhaustive enumeration')
Path(__file__).resolve().with_name('independent24.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
