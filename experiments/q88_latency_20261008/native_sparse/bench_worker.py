#!/usr/bin/env python3
"""Measure an isolated production CRT compiler profile with a captured witness.
Private fixtures stay on the experiment server; outputs contain public data.
"""
import argparse, hashlib, json, os, struct, time
from pathlib import Path
import revac_acc96_memver_worker as worker
from revac_acc96_wire import PROVER_BYTES, POLY_BYTES, INSTANCE_MAGIC, Q19, InstancePayload, LinkPayload, ProverPayload

def check_rejections(instance, proof, entry):
    """Check the compiled public relation and canonical envelope separately."""
    changed_rhs=bytearray(instance)
    changed_rhs[len(INSTANCE_MAGIC)+8]^=1  # first accumulator coefficient
    changed_link=bytearray(instance)
    link_offset=len(INSTANCE_MAGIC)+8+POLY_BYTES+3*32
    coefficient=struct.unpack_from('<I',changed_link,link_offset)[0]
    struct.pack_into('<I',changed_link,link_offset,(coefficient+1)%Q19)
    for name, public_input, envelope in (
        ('changed accumulator RHS',bytes(changed_rhs),proof),
        ('changed canonical q19 link RHS',bytes(changed_link),proof),
        ('truncated proof envelope',instance,proof[:-1]),
        ('trailing proof envelope data',instance,proof+b'\x00'),
    ):
        try:
            ok,_=worker.verify(public_input+struct.pack('<I',len(envelope))+envelope,entry)
        except ValueError as exc:
            print('NEGATIVE_PASS '+name+': '+str(exc),flush=True)
        else:
            if ok: raise ValueError('tampering unexpectedly accepted: '+name)
            print('NEGATIVE_PASS '+name+': Pack Verify returned false',flush=True)
    ok,_=worker.verify(instance+struct.pack('<I',len(proof))+proof,entry)
    if not ok: raise ValueError('valid proof failed after negative checks')
    print('VALID_REBIND_PASS',flush=True)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('fixture',type=Path)
    ap.add_argument('--repetitions',type=int,default=3); ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--parameters-only',action='store_true'); ap.add_argument('--verify-artifact',type=Path)
    ap.add_argument('--expect-rejection',action='store_true')
    ap.add_argument('--negative-checks',action='store_true'); a=ap.parse_args()
    if a.expect_rejection and not a.verify_artifact: ap.error('--expect-rejection requires --verify-artifact')
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
        try:
            ok,_=worker.verify(instance+struct.pack('<I',len(proof))+proof,entry)
        except ValueError as exc:
            if not a.expect_rejection: raise
            print('EXPECTED_REJECTION '+str(exc),flush=True);return
        if a.expect_rejection:
            if ok: raise ValueError('invalid profile unexpectedly accepted')
            print('EXPECTED_REJECTION Pack Verify returned false',flush=True);return
        if not ok: raise ValueError('cross-backend verify failed')
        print('CROSS_BACKEND_VERIFY_PASS',flush=True)
        if a.negative_checks: check_rejections(instance,proof,entry)
        return
    started=time.perf_counter(); worker.prepare_prover(dynamic,entry); prover_prep=time.perf_counter()-started
    ps=entry['prove_statement'].ps
    if int(ps.simple_verify(ps.statement,ps.witness))!=1: raise ValueError('clear verify failed')
    if a.parameters_only:
        from lazer import lib
        lib.labrador38_py_free_params(ps.params)
        if ps.gen_params(ps.params,ps.statement,1,1): raise ValueError('parameter audit failed')
        print('PARAMETERS_CLEAR_PASS',flush=True)
        return
    entry['last_prover']=prover
    compact=worker.REUSE_MAGIC+struct.pack('<II',prover.sequence,prover.handle)+prover.raw[16:16+POLY_BYTES]+dynamic[PROVER_BYTES:]
    rows=[]
    for i in range(a.repetitions):
        started=time.perf_counter(); proof,_=worker.prove(compact,entry,reuse=True); generation=time.perf_counter()-started
        started=time.perf_counter(); ok,_=worker.verify(instance+struct.pack('<I',len(proof))+proof,entry); verification=time.perf_counter()-started
        if not ok: raise ValueError('Pack/reload verify failed')
        rows.append({'generation_seconds':generation,'verification_seconds':verification,'proof_bytes':len(proof)})
        print('MEASURED '+json.dumps(rows[-1]),flush=True)
    if a.negative_checks: check_rejections(instance,proof,entry)
    group_size=int(os.environ.get('REVAC_Q88_NORM_GROUP_SIZE','1'))
    a.output.write_text(json.dumps({'profile':('production CRT baseline' if group_size==1 else 'production CRT isolated accumulator norm grouping v2'),
        'norm_group_size':group_size,'approx_order':os.environ.get('REVAC_Q88_APPROX_ORDER','baseline'),'fixture':str(a.fixture),
        'prewarm_prover_seconds':prover_prep,'prewarm_verifier_seconds':verifier_prep,
        'clear_verify':True,'pack_reload_verify':True,'negative_checks':a.negative_checks,'measurements':rows},indent=2)+'\n')
    a.output.with_suffix('.proof.bin').write_bytes(proof)
if __name__=='__main__': main()
