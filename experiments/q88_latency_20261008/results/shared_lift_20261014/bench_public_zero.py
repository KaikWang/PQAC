import hashlib,json,os,struct,time,sys
from pathlib import Path
import revac_acc96_memver_worker as worker
from revac_acc96_wire import PROVER_BYTES,POLY_BYTES,INSTANCE_MAGIC,InstancePayload,LinkPayload,ProverPayload
fixture=Path(sys.argv[1]);proof=Path(sys.argv[2]).read_bytes()
public=(fixture/'public.bin').read_bytes();dynamic=(fixture/'prover_link.bin').read_bytes()
entry={'pp_id':hashlib.shake_256(public).digest(worker.HASH_BYTES),'public_raw':public,'public':worker.PublicParameters.decode(public),'prove_statement':worker.Q96StatementCache(),'verify_statement':worker.Q96StatementCache()}
prover=ProverPayload(dynamic[:PROVER_BYTES]);link=LinkPayload.decode(dynamic[PROVER_BYTES:])
instance=(INSTANCE_MAGIC+struct.pack('<II',prover.sequence,link.window_id)+prover.raw[16:16+POLY_BYTES]+link.context_hash+link.parameter_digest+link.accumulator_digest+b''.join(struct.pack('<256I',*row) for row in link.commitment))
worker.prepare_verifier(instance,entry)
payload=instance+struct.pack('<I',len(proof))+proof
# Audit the optimized public zero metadata, then interleave warm repeated
# verification of the SAME proof. This is a diagnostic, not a Show replay.
os.environ['REVAC_Q88_PUBLIC_ZERO_TERMS']='2'
assert worker.verify(payload,entry)[0]
rows=[]
for repeat in range(5):
 for mode in ((0,1) if repeat%2==0 else (1,0)):
  os.environ['REVAC_Q88_PUBLIC_ZERO_TERMS']=str(mode)
  started=time.perf_counter();assert worker.verify(payload,entry)[0]
  rows.append({'repeat':repeat,'mode':mode,'seconds':time.perf_counter()-started})
Path(sys.argv[3]).write_text(json.dumps(rows,indent=2)+'\n')
print('PUBLIC_ZERO_INTERLEAVED_SAME_PROOF_PASS')
