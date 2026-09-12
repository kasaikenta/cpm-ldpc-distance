"""Independent finite checks of the archived classical search and cofactor words.
Uses Python's own GF(2) elimination and full kernel enumeration as the oracle.
No run-record assertion is used as an oracle. These finite checks are not a proof
of the general theorem or of the distances of the large table examples.
"""
from pathlib import Path
import hashlib, itertools, json, math, random, subprocess, tempfile, time
ROOT=Path(__file__).resolve().parent

def rows(E,P):
 L=len(E[0]); return [sum(1<<(t*L+l) for l in range(L) for t in range(P) if (t+er[l])%P==s) for er in E for s in range(P)]
def exact(E,P):
 H=rows(E,P); n=P*len(E[0]); piv={}
 for h in H:
  while h:
   i=h.bit_length()-1
   if i in piv:h^=piv[i]
   else:piv[i]=h;break
 basis=[]
 for j in range(n):
  if j in piv:continue
  v=1<<j
  for i in sorted(piv):
   if (piv[i]&v).bit_count()%2:v^=1<<i
  assert all((h&v).bit_count()%2==0 for h in H)
  basis.append(v)
 assert len(basis)<=18
 words=[0]
 for v in basis:words += [w^v for w in words]
 return min(w.bit_count() for w in words[1:]),len(words)-1

def run(exe,dr,E,P,cut,extra=()):
 L=len(E[0]);inp=dr/'case.txt';out=dr/'out.json'
 inp.write_text('CPM_CSS_V1\nreview\n%d %d %d 0\n'%(P,L,len(E))+'\n'.join(' '.join(map(str,row)) for row in E)+'\n')
 cmd=[str(exe),'--input',str(inp),'--output',str(out),'--mode','classical','--side','X','--max-weight',str(cut),'--threads','1','--quiet',*extra]
 subprocess.run(cmd,check=True,capture_output=True,timeout=60)
 d=json.loads(out.read_text())
 assert d['mode']=='classical' and d['stabilizer_rows']==0
 if d['found']:
  v=d['witness']['support'];assert 0<len(v)<=cut and len(v)==len(set(v))
  w=sum(1<<i for i in v); assert all((h&w).bit_count()%2==0 for h in rows(E,P))
 else:assert d['complete']
 return d

def main():
 start=time.time();rng=random.Random(20260912);cases=[]
 for J in (2,3,4):
  for L in (J+1,J+2):
   cases.append(([[0]*L for _ in range(J)],2))
   for P in (2,3,4,5):
    for _ in range(2):cases.append(([[rng.randrange(P) for l in range(L)] for j in range(J)],P))
 count=words=0;sharded=0
 with tempfile.TemporaryDirectory() as td:
  dr=Path(td);exe=dr/'distance'
  subprocess.run(['c++','-std=c++17','-O2','-pthread',str(ROOT/'distance.cpp'),'-o',str(exe)],check=True,capture_output=True)
  chosen=[]
  for E,P in cases:
   d,nwords=exact(E,P);words+=nwords
   for cut in (max(1,d-2),d):
    for extra in ((),('--no-qc-reduce',)):
     r=run(exe,dr,E,P,cut,extra);assert r['found']==(cut>=d),(E,P,d,cut,extra,r)
     count+=1
   if d>=4 and len(chosen)<3:chosen.append((E,P,d))
  for ci,(E,P,d) in enumerate(chosen):
   for cut in (d-2,d):
    found=False
    for root in range(len(E[0])):
     for shard in range(3):
      stem=dr/f'shard_{ci}_{cut}_{root}_{shard}'
      extra=('--root-start',str(root),'--root-end',str(root+1),'--checkpoint',str(stem)+'.chk','--progress',str(stem)+'.json','--shard-count','3','--shard-index',str(shard),'--split-weight',str(min(3,cut)))
      r=run(exe,dr,E,P,cut,extra);found|=r['found'];sharded+=1
    assert found==(cut>=d)
 # Sparse cofactor words for the enormous explicit exponents: no expanded H.
 cofactors=[]
 for J in (2,3,4):
  for L in (J+1,J+3):
   U=math.factorial(J+1);E=[[U**(i*L+l) for l in range(L)] for i in range(J)]
   R=U**((J-1)*L)*(U**(L-1)-1)
   for P in ((U-2)*R+1,(U-2)*R+2):
    support=[]
    for omit in range(J+1):
     cols=[l for l in range(J+1) if l!=omit]
     exps=[sum(E[i][perm[i]] for i in range(J))%P for perm in itertools.permutations(cols)]
     assert len(exps)==len(set(exps))==math.factorial(J)
     support.extend((omit,t) for t in exps)
    assert len(support)==U
    for er in E:
     parity=set()
     for l,t in support:
      a=(t+er[l])%P
      if a in parity:parity.remove(a)
      else:parity.add(a)
     assert not parity
    cofactors.append({'J':J,'L':L,'P':str(P),'weight':U})
 out={'seed':20260912,'small_instances':len(cases),'independent_nonzero_kernel_words':words,'unsharded_solver_calls':count,'sharded_solver_calls':sharded,'cofactor_checks':cofactors,'solver_sha256':hashlib.sha256((ROOT/'distance.cpp').read_bytes()).hexdigest(),'elapsed_sec':time.time()-start,'all_passed':True}
 (ROOT/'review_checks.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k!='cofactor_checks'}))
if __name__=='__main__':main()
