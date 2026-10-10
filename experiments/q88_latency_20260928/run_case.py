"""Isolated full-show benchmark. Never changes the production socket."""
import argparse, json, os, pathlib, subprocess, time
p=argparse.ArgumentParser()
p.add_argument('name'); p.add_argument('--root', default='/root/autodl-tmp/pq_ac')
p.add_argument('--reps',type=int,default=4); p.add_argument('settings',nargs='*')
a=p.parse_args(); root=pathlib.Path(a.root)
out=root/'OURS/experiments/q88_latency_20260928/results'; out.mkdir(parents=True,exist_ok=True)
sock=f'/tmp/revac-audit-{os.getpid()}.sock'
env=os.environ.copy()
for key in list(env):
    if key.startswith(('REVAC_', 'OMP_', 'GOMP_')): del env[key]
env.update(PYTHONPATH=f'{root}/OURS/code/entities/zk_lazer:{root}/lazer/python',
           LD_LIBRARY_PATH=f'{root}/lazer', OMP_NUM_THREADS='12',
           REVAC_LABRADOR_Q96_SOCKET=sock,REVAC_LABRADOR_ONLINE='1',
           REVAC_LABRADOR_Q96_PREWARM='1',REVAC_LABRADOR_USE_SDE='0',
           REVAC_Q96_SHOW_REPETITIONS=str(a.reps),REVAC_LABRADOR_CUDA_AUDIT='1',
           REVAC_LABRADOR_CUDA_LIBRARY=f'{root}/lazer/librevac_labrador_cuda.so')
for setting in a.settings:
    key,value=setting.split('=',1); env[key]=value
(out/f'{a.name}.env.json').write_text(json.dumps({k:v for k,v in env.items() if k.startswith(('REVAC_', 'OMP_', 'GOMP_'))},indent=2))
with (out/f'{a.name}.worker.log').open('w') as log:
    worker=subprocess.Popen(['/root/miniconda3/envs/sage/bin/python',str(root/'OURS/code/entities/zk_lazer/revac_acc96_memver_worker.py'),'--socket',sock],env=env,stdout=log,stderr=subprocess.STDOUT)
    try:
        deadline=time.monotonic()+20
        while not pathlib.Path(sock).exists():
            if worker.poll() is not None or time.monotonic()>deadline: raise RuntimeError('worker start failed')
            time.sleep(.1)
        with (out/f'{a.name}.show.log').open('w') as demo:
            result=subprocess.run([str(root/'OURS/code/entities/build-pqac8/revocable_ac_demo')],env=env,stdout=demo,stderr=subprocess.STDOUT,timeout=500)
        print(a.name,'exit',result.returncode,flush=True)
    finally:
        worker.terminate(); worker.wait(timeout=20)
for line in (out/f'{a.name}.show.log').read_text().splitlines():
    if 'q96 show 第' in line: print(line)
for line in (out/f'{a.name}.worker.log').read_text().splitlines():
    if '[CUDA audit]' in line: print(line)
