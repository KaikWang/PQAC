#!/usr/bin/env python3
"""Measure unchanged production CRT worker using a captured local witness.
Private fixtures stay on the experiment server; outputs contain public data.
"""
import argparse, hashlib, json, struct, time
from pathlib import Path
import revac_acc96_memver_worker as worker
from revac_acc96_wire import PROVER_BYTES, POLY_BYTES, INSTANCE_MAGIC, InstancePayload, LinkPayload, ProverPayload

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('fixture',type=Path)
    ap.add_argument('--repetitions',type=int,default=3); ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--verify-artifact',type=Path); a=ap.parse_args()
    public=(a.fixture/'public.bin').read_bytes(); dynamic=(a.fixture/'prover_link.bin').read_bytes()
    entry={'pp_id':hashlib.shake_256(public).digest(worker.HASH_BYTES),'public_raw':public,
           'public':worker.PublicParameters.decode(public),'prove_statement':worker.Q96StatementCache(),
           'verify_statement':worker.Q96StatementCache()}
    prover=ProverPayload(dynamic[:PROVER_BYTES]); link=LinkPayload.decode(dynamic[PROVER_BYTES:])
    instance=(INSTANCE_MAGIC+struct.pack('<II',prover.sequence,link.window_id)+prover.raw[16:16+POLY_BYTES]+
              link.context_hash+link.parameter_digest+link.accumulator_digest+
              b''.join(struct.pack('<256I',*row) for row in link.commitment))
    InstancePayload.decode(instance)
    started=time.perf_counter(); worker.prepare_verifier(instance,entry); verifier_prep=time.perf_counter()-started
    if a.verify_artifact:
        proof=a.verify_artifact.read_bytes()
        ok,_=worker.verify(instance+struct.pack('<I',len(proof))+proof,entry)
        if not ok: raise ValueError('cross-backend verify failed')
        print('CROSS_BACKEND_VERIFY_PASS',flush=True);return
    started=time.perf_counter(); worker.prepare_prover(dynamic,entry); prover_prep=time.perf_counter()-started
    ps=entry['prove_statement'].ps
    if int(ps.simple_verify(ps.statement,ps.witness))!=1: raise ValueError('clear verify failed')
    entry['last_prover']=prover
    compact=worker.REUSE_MAGIC+struct.pack('<II',prover.sequence,prover.handle)+prover.raw[16:16+POLY_BYTES]+dynamic[PROVER_BYTES:]
    rows=[]
    for i in range(a.repetitions):
        started=time.perf_counter(); proof,_=worker.prove(compact,entry,reuse=True); generation=time.perf_counter()-started
        started=time.perf_counter(); ok,_=worker.verify(instance+struct.pack('<I',len(proof))+proof,entry); verification=time.perf_counter()-started
        if not ok: raise ValueError('Pack/reload verify failed')
        rows.append({'generation_seconds':generation,'verification_seconds':verification,'proof_bytes':len(proof)})
        print('MEASURED '+json.dumps(rows[-1]),flush=True)
    a.output.write_text(json.dumps({'profile':'production CRT unchanged','fixture':str(a.fixture),
        'prewarm_prover_seconds':prover_prep,'prewarm_verifier_seconds':verifier_prep,
        'clear_verify':True,'pack_reload_verify':True,'measurements':rows},indent=2)+'\n')
    a.output.with_suffix('.proof.bin').write_bytes(proof)
if __name__=='__main__': main()
