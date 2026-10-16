"""Create a standalone sequence benchmark from the current lifecycle demo.
Production source is read only; all samples prove and verify full packets.
"""
from pathlib import Path
base = Path(__file__).resolve().parents[3] / 'code/entities/test/revocable_ac_demo.c'
s = base.read_text()
helper = r'''
static int sequence_show_samples(const revac_public_state_t *state,
    revac_uav_t *uav, const revac_gcs_t *gcs, revac_show_context_t *context,
    revac_gcs_verify_guard_t *guard, const uint8_t msg[PARAM_M * PARAM_N / 8]) {
  for (unsigned run=1; run<=3; run++) {
    revac_show_proof_t candidate;
    uint8_t nonce[REVAC_SHOW_NONCE_BYTES];
    uint8_t *wire=NULL;
    size_t bytes=0;
    double started, generated, verified;
    int ok=0;
    revac_show_proof_init(&candidate);
    randombytes(nonce,sizeof(nonce));
    started=monotonic_seconds();
    if (!revac_uav_show_prove(&candidate,state,uav,msg,nonce) ||
        !revac_show_proof_pack_compact(&wire,&bytes,&candidate,1)) goto sample_done;
    generated=monotonic_seconds()-started;
    started=monotonic_seconds();
    if (!revac_gcs_verify_show_wire_guarded_parallel(state,gcs,context,guard,wire,bytes))
      goto sample_done;
    verified=monotonic_seconds()-started;
    printf("[q88 sequence Show measurement] {\"sequence\":%u,\"run\":%u,"
           "\"generation_seconds\":%.9f,\"verification_seconds\":%.9f,"
           "\"accumulator_bytes\":%zu,\"packet_bytes\":%zu}\n",
           state->acc96_seq,run,generated,verified,candidate.acc_zk_proof_len,bytes);
    ok=1;
sample_done:
    free(wire);
    revac_show_proof_clear(&candidate);
    if (!ok) return 0;
  }
  return 1;
}
'''
s=s.replace('int main(void) {', helper+'\nint main(void) {',1)
marker='  exact_lift = malloc(sizeof(*exact_lift));'
s=s.replace(marker,'''  if (revac_labrador_online_enabled() && !sequence_show_samples(
      &public_state,&uav2,&gcs,&show_context,&verify_guard,msg)) goto cleanup;
'''+marker,1)
marker='  ok = 1;\n  printf("q88 候选 D 实体接入测试通过'
extra=r'''  if (revac_labrador_online_enabled()) {
    for (unsigned step=0;step<2;step++) {
      revac_uav_t extra;
      acc96_relation_update_t next_update;
      uint8_t *prefetch=NULL,*manifest=NULL;
      size_t prefetch_len=0,manifest_len=0;
      int step_ok=0;
      revac_uav_init(&extra);
      uint32_t revoke_handle=uav_local.handle;
      if(step==1){
        uint64_t id=UINT64_C(700000);
        do {
          revac_uav_keygen(&extra,id++,ta.ta_pk.seed);
        } while(extra.handle==uav1.handle || extra.handle==uav2.handle ||
                extra.handle==uav_local.handle);
        if(!revac_ta_register_uav(&ta,&extra) || !revac_ta_issue_uav(&ta,&extra,msg))
          goto extra_done;
        revoke_handle=extra.handle;
      }
      if(!revac_ta_export_q96_tail4_prefetch(&ta,revoke_handle,&prefetch,&prefetch_len) ||
         !revac_uav_cache_q96_tail4_prefetch(&uav2,&public_state,prefetch,prefetch_len) ||
         !revac_ta_revoke_uav_q96(&ta,revoke_handle,&next_update)) goto extra_done;
      revac_ta_public_state(&public_state,&ta);
      if(!revac_ta_export_q96_update_manifest(&ta,next_update.sequence,&manifest,&manifest_len) ||
         !revac_uav_apply_update_q96_manifest(&uav2,&public_state,manifest,manifest_len) ||
         !revac_uav_memver_q96(&uav2,&public_state)) goto extra_done;
      revac_show_context_from_public(&show_context,&public_state,uav2.pk.seed);
      if(!sequence_show_samples(&public_state,&uav2,&gcs,&show_context,&verify_guard,msg))
        goto extra_done;
      step_ok=1;
extra_done:
      free(prefetch);free(manifest);revac_uav_clear(&extra);
      if(!require(step_ok,"连续撤销更新与完整 Show 失败")) goto cleanup;
    }
  }
'''
assert marker in s
s=s.replace(marker,extra+marker,1)
Path(__file__).with_name('sequence_show_demo.c').write_text(s)
