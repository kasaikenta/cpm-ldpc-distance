"""Recompute a reported lower bound with resumable classical DFS.
Example: python3 verification/reproduce.py --J 4 --L 6 --threads 4 --seconds 60
Repeat the same command to continue. L=5,J=4 uses its fast full enumeration.
"""
import argparse,collections,concurrent.futures,hashlib,json,os,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--J',type=int,choices=(3,4),required=True);p.add_argument('--L',type=int,choices=range(5,9),required=True);p.add_argument('--threads',type=int,default=1);p.add_argument('--seconds',type=int,default=60);p.add_argument('--shards',type=int,default=32);p.add_argument('--output-dir',type=Path,default=ROOT/'verification/recomputed');a=p.parse_args()
assert a.threads>0 and a.seconds>0 and a.shards>0
if (a.J,a.L)==(4,5):
 subprocess.run([os.sys.executable,str(ROOT/'verification/enumerate_j4_l5.py')],check=True);raise SystemExit(0)
c=next(c for c in json.loads((ROOT/'joint_examples.json').read_text()) if (c['J'],c['L'])==(a.J,a.L));cutoff=(24 if a.J==3 else c['d_lower'])-2
case_name=c.get('name',f'J{a.J}_L{a.L}_P{c["P"]}')
d=a.output_dir.resolve()/f'J{a.J}_L{a.L}_P{c["P"]}_w{cutoff}_s{a.shards}';d.mkdir(parents=True,exist_ok=True)
source=ROOT/'verification/distance.cpp';sha=hashlib.sha256(source.read_bytes()).hexdigest();assert sha=='0dd279be2e41cd81ce7d70b767854012bd106116d57fd455dc9262bbb1aa3620'
exe=d/'distance'
if not exe.exists():subprocess.run(['c++','-std=c++17','-O3','-pthread',str(source),'-o',str(exe)],check=True)
f=d/'instance.txt';f.write_text('CPM_CSS_V1\n'+case_name+' %d %d %d 0\n'%(c['P'],c['L'],c['J'])+'\n'.join(' '.join(map(str,r)) for r in c['E'])+'\n')
def read(p):
 try:return json.loads(p.read_text())
 except (OSError,ValueError):return {}
def atomic(path,value):
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n');os.replace(tmp,path)
units=[(r,s) for r in range(c['L']) for s in range(a.shards)];old=read(d/'summary.json');start=time.time();deadline=time.monotonic()+a.seconds;procs={};stop=False
for r,s in units:
 q=read(d/f'r{r}_s{s}.result.json')
 if q.get('found'):raise RuntimeError('Stored counterexample: inspect its result file.')
def is_done(u):return read(d/f'r{u[0]}_s{u[1]}.result.json').get('complete',False)
def shutdown(*args):
 global stop
 stop=True
 for proc in list(procs.values()):
  if proc.poll() is None:proc.send_signal(signal.SIGTERM)
for sig in (signal.SIGINT,signal.SIGTERM):signal.signal(sig,shutdown)
def save():
 complete=[list(u) for u in units if is_done(u)];missing=[list(u) for u in units if list(u) not in complete];stats=[read(d/f'r{r}_s{s}.progress.json') for r,s in units]
 v=dict(experiment_id=case_name,status='complete' if not missing else 'timed_out' if stop or time.monotonic()>=deadline else 'running',side='classical',max_weight=cutoff,requested_scope=[list(u) for u in units],completed_scope=complete,incomplete_scope=missing,current_cursor_or_frontier=missing,units_processed=sum(q.get('units_processed',0) for q in stats),witness_found=False,best_weight=None,started_at=old.get('started_at',start),updated_at=time.time(),elapsed_sec=time.time()-start,h_rt_sec=None,soft_deadline_sec=a.seconds,last_checkpoint_at=time.time(),resume_count=old.get('resume_count',-1)+1,termination_reason='search_exhausted' if not missing else 'time_budget' if stop or time.monotonic()>=deadline else None,no_witness_through=cutoff if not missing else None,source_sha256=sha)
 atomic(d/'summary.json',v);return v
def work(u):
 r,s=u;stem=d/f'r{r}_s{s}';left=int(deadline-time.monotonic())
 if stop or left<1:return
 cmd=[str(exe),'--input',str(f),'--output',str(stem.with_suffix('.result.json')),'--mode','classical','--side','X','--max-weight',str(cutoff),'--root-start',str(r),'--root-end',str(r+1),'--threads','1','--shard-count',str(a.shards),'--shard-index',str(s),'--split-weight',str(min(5,cutoff)),'--checkpoint',str(stem.with_suffix('.chk')),'--progress',str(stem.with_suffix('.progress.json')),'--soft-deadline-sec',str(min(60,left)),'--checkpoint-interval-sec','2','--quiet']
 proc=subprocess.Popen(cmd,stdout=subprocess.DEVNULL);procs[u]=proc;ret=proc.wait();procs.pop(u,None)
 if ret:raise RuntimeError(f'DFS failed for {u}; saved files are retained.')
 if read(stem.with_suffix('.result.json')).get('found'):raise RuntimeError(f'Counterexample for {u}; inspect the result file.')
pending=collections.deque(u for u in units if not is_done(u));active={};last=0
with concurrent.futures.ThreadPoolExecutor(max_workers=a.threads) as pool:
 while pending or active:
  while pending and not stop and time.monotonic()+1<deadline and len(active)<a.threads:
   u=pending.popleft();active[pool.submit(work,u)]=u
  if time.time()-last>5:save();last=time.time()
  if not active:break
  finished,_=concurrent.futures.wait(active,timeout=1,return_when=concurrent.futures.FIRST_COMPLETED)
  for fut in finished:
   u=active.pop(fut);fut.result()
   if not is_done(u) and time.monotonic()+1<deadline and not stop:pending.append(u)
stop=True;v=save();print(json.dumps({k:v[k] for k in ['status','max_weight','units_processed','no_witness_through']},ensure_ascii=False));print('Complete partitions:',len(v['completed_scope']),'/',len(units));print('Saved state:',d)
