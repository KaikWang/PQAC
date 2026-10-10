#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "pack.h"
#include "proofsystem.h"
#include "labrador_core.h"
#include "labradoodle.h"
#include "labrador.h"
#include "labrador_tail.h"
#include "malloc.h"
#include "timing.h"
#include "revac_profile.h"

static const char *lab_round_name(const lab_params pp) {
  if (pp->compressed) return "Labradoodle";
  if (pp->tail) return "LaBRADOR-tail";
  return "LaBRADOR";
}

static void profile_lab_round(const char *mode, size_t round,
                              const lab_params pp, double elapsed) {
  if (!revac_labrador_profile_detail_enabled()) return;
  fprintf(stderr,
          "[Labrador 原生细分] %s 第%zu轮 %s：输入向量=%zu，输入多项式=%zu，最大分组=%zu，%.6f 秒\n",
          mode, round + 1, lab_round_name(pp), pp->r_old, pp->nn,
          pp->nmax, elapsed);
}

static int pack_debug_verify_enabled(void) {
  const char *value = getenv("REVAC_LABRADOR_DCH_VERIFY");
  return value != NULL && value[0] != '\0' && value[0] != '0';
}

/* Pack parameter generation historically continued a folding round only when
 * the combined proof/witness size improved by at least ten percent.  Expose
 * that engineering heuristic for controlled A/B experiments: every emitted
 * round remains a normal LaBRADOR/Labradoodle round, while a higher threshold
 * stops earlier and trades proof bytes for latency. */
static double pack_improve_threshold(void) {
  static int initialized = 0;
  static double threshold = 0.90;
  if (!initialized) {
    const char *value = getenv("REVAC_LABRADOR_PACK_IMPROVE_THRESHOLD");
    if (value != NULL && value[0] != '\0') {
      char *end = NULL;
      const double parsed = strtod(value, &end);
      if (end != value && *end == '\0' && parsed > 0.0 && parsed < 1.0)
        threshold = parsed;
    }
    initialized = 1;
  }
  return threshold;
}

static int pack_split_enabled(void) {
  const char *value = getenv("REVAC_LABRADOR_PACK_SPLIT");
  if (value == NULL || value[0] == '\0') return 1;
  return value[0] != '0';
}

static _Thread_local int last_prove_relation_valid = 1;

int pack_last_prove_relation_valid(void) {
  return last_prove_relation_valid;
}

static void pack_check_prove_relation(const char *kind, size_t round,
                                      const statement st,
                                      const witness wt) {
  if (pack_debug_verify_enabled()) {
    const int valid = verify(st, wt);
    fprintf(stderr, "[Labrador Pack 诊断] 第%zu轮 %s输出关系：%s\n",
            round + 1, kind, valid ? "通过" : "失败");
  }
}

void pack_proof_init(pack_proof pi, const pack_params pp){
  pi->np = pp->np;
  pi->p = pp->np ? _malloc(pi->np * sizeof(lab_proof)) : NULL;
  pi->zkp = pp->zkp != NULL ? _malloc(sizeof(lnp_proof)) : NULL;
}

void pack_proof_free(pack_proof pi){
  size_t i;
  for(i=0;i<pi->np;i++){
    lab_proof_free(pi->p[i]);
  }
  witness_free(pi->owt);
  free(pi->p);
  pi->np = 0;
  free(pi->zkp);
  pi->zkp = NULL;
}

void pack_params_gen(
  pack_params pp, 
  size_t *pibits, 
  const statement st, 
  int zk, 
  size_t iwtbits
)
{
  size_t maxnp = 20;
  size_t i, pisize[maxnp+1], wtsize[maxnp+1+1];
  int compressed, tail, iszk, split, global_norm_check, transition, improve;
  int ret;
  double slack_norm_check;
  const double improve_threshold = pack_improve_threshold();
  statement sts[maxnp+1+1];
  lab_params *cpp;

  pp->np = 0;
  pp->p = _malloc(maxnp * sizeof(lab_params));
  pp->zkp = zk ? _malloc(sizeof(lnp_params)) : NULL;
  pp->zkround = SIZE_MAX;

  // init statement with relevant info from input statement
  statement_init(sts[0], st->r, st->r);
  for(i=0;i<st->r;i++){
    sts[0]->n[i] = st->n[i];
    sts[0]->normsq[i] = st->normsq[i];
    sts[0]->normty[i] = st->normty[i];
  }
  sts[0]->zqcnst->nsigmam1 = st->zqcnst->nsigmam1;
  wtsize[0] = iwtbits;


  // labradoodle
  transition = 0;
  compressed = 1;
  tail = 0;
  iszk = zk;
  split = pack_split_enabled();

  do {

  while(pp->np < maxnp){
    if(!transition){
      global_norm_check = 0;
      slack_norm_check = JL_INF_SLACK;
    }
    else if(iszk){
      global_norm_check = 0;
      slack_norm_check = 2*JL_INF_SLACK;
    }
    else{
      global_norm_check = 1;
      slack_norm_check = JL_L2_SLACK;
    }
    cpp = &pp->p[pp->np];

    ret = lab_params_gen(*cpp, &pisize[pp->np], &wtsize[pp->np+1], sts[pp->np], 
                         compressed, tail, iszk, split, global_norm_check, 
                         slack_norm_check);
    improve= ((double)(pisize[pp->np] + wtsize[pp->np+1])) <
             improve_threshold*wtsize[pp->np];


    if(ret || (!improve && !transition && pp->np > 0 && pp->zkround != pp->np-1)){
      if(!ret){ 
        lab_params_free(*cpp);
      }
      if(pp->np == 0){
        if (zk){
          printf("ERROR: No Labradoodle params found, ZK not possible\n");
          exit(1);
        }
        break;
      }
      if(pp->zkround == pp->np-1){
        printf("ERROR: No Labradoodle round possible after LNP");
        exit(1);
      }
      transition = 1;
      sts[pp->np]->zqcnst->nsigmam1 = 0;
      statement_free(sts[pp->np]);
      lab_params_free(pp->p[pp->np-1]);
      pp->np--;
      continue;
    }

    pp->np++;

    statement_init(sts[pp->np], 5, 5);
    sts[pp->np]->n[0] = (*cpp)->nmax;
    sts[pp->np]->n[1] = (*cpp)->nmax;
    sts[pp->np]->n[2] = (*cpp)->len[LAB_INCOM] + (*cpp)->len[LAB_QUADG] 
                        + (*cpp)->len[LAB_LING] + (*cpp)->len[LAB_OUTCOM];
    sts[pp->np]->n[3] = (*cpp)->len[LAB_BIN];
    sts[pp->np]->n[4] = sts[pp->np]->n[3]; // include sigmam1
    sts[pp->np]->zqcnst->nsigmam1 = 1;
    for(i=0;i<4;i++){
      sts[pp->np]->normsq[i] = (*cpp)->normsq_new[i];
    }
    sts[pp->np]->normsq[4] = sts[pp->np]->normsq[3];

    if(transition){
      break;
    }
  }

  // lnp
  if (iszk) {    
    ret = lnp_params_gen (pp->zkp[0], &pisize[pp->np], &wtsize[pp->np+1], sts[pp->np]);
    if(ret && pp->np < 10){
      continue; // try another labradoodle round
    }
    else if (ret){
      printf("ERROR: No LNP params found, ZK not possible\n");
      exit(1);
    }

    memset (pp->p[pp->np], 0, sizeof(pp->p[pp->np])); // empty lab proof
    pp->zkround = pp->np;
    pp->np++;
    lnp_statement_init (sts[pp->np], sts[pp->np-1], pp->zkp[0]);
    sts[pp->np]->normsq[sts[pp->np]->r] = N * pp->zkp[0]->xbinlen;
    sts[pp->np]->normty[sts[pp->np]->r] = L2APPROX;
    //sts[pp->np]->rqcnst->nsparse = 1;
    //sts[pp->np]->rqcnst->ncom = 4 + 2;
    //sts[pp->np]->zqcnst->nsparse = LIFTS + 1;
    sts[pp->np]->zqcnst->nsigmam1 = 1;
    sts[pp->np]->r += 1;

    iszk = 0;
    transition = 0;
  } else {
    break;
  }
  
  } while (1);


  // labrador
  compressed = 0;
  tail = 0;
  iszk = 0;
  split = pack_split_enabled();
  global_norm_check = 1;
  slack_norm_check = JL_L2_SLACK;
  while(pp->np < maxnp){
    cpp = &pp->p[pp->np];

    ret = lab_params_gen(*cpp, &pisize[pp->np], &wtsize[pp->np+1], sts[pp->np], 
                         compressed, tail, iszk, split, global_norm_check, 
                         slack_norm_check);
    improve= ((double)(pisize[pp->np] + wtsize[pp->np+1])) <
             improve_threshold*wtsize[pp->np];

    if(ret || (!improve && pp->np > 0 && !pp->p[pp->np-1]->compressed)){
      if(!ret){ 
        lab_params_free(*cpp);
      }
      break;
    }

    pp->np++;

    statement_init(sts[pp->np], (*cpp)->fz + 1, (*cpp)->fz + 1);
    sts[pp->np]->n[0] = (*cpp)->nmax;
    sts[pp->np]->n[(*cpp)->fz - 1] = (*cpp)->nmax;
    sts[pp->np]->n[(*cpp)->fz] = (*cpp)->len[LAB_INCOM] + (*cpp)->len[LAB_QUADG] 
                                 + (*cpp)->len[LAB_LING];
    for(i=0;i<sts[pp->np]->r;i++){
      sts[pp->np]->normsq[i] = (*cpp)->normsq_new[i];
    }
  }

  // tail
  if(pp->np < maxnp){
    compressed = 0;
    tail = 1;
    iszk = 0;
    split = pack_split_enabled();
    global_norm_check = 0;
    slack_norm_check = 1;
    cpp = &pp->p[pp->np];

    ret = lab_params_gen(*cpp, &pisize[pp->np], &wtsize[pp->np+1], sts[pp->np], 
                          compressed, tail, iszk, split, global_norm_check, 
                          slack_norm_check);
    if(!ret)
      pp->np++;
  }

  *pibits = 0;
  for(i=0;i<pp->np;i++){
    *pibits += pisize[i];
  }
  *pibits += wtsize[pp->np];

  for(i=0;i<pp->np;i++){
    sts[i]->zqcnst->nsigmam1 = 0;
    statement_free(sts[i]);
  }
  if(pp->np == 0 || !pp->p[pp->np-1]->tail){
    sts[pp->np]->zqcnst->nsigmam1 = 0;
    statement_free(sts[pp->np]);
  }
}

void pack_params_print(const pack_params pp){
  size_t i;
  for(i=0;i<pp->np;i++){
    printf("\n*** Round %2zu ***\n", i+1);
    if(pp->zkp != NULL && i == pp->zkround){
      printf("LNP round\n");
      lnp_params_print(pp->zkp[0]);
    }
    else{
    lab_params_print(pp->p[i]);
    }
  }
  if(pp->np == 0){
    printf("No rounds: direct statement-witness check\n");
  }
}

void pack_params_free(pack_params pp){
  size_t i;
  for(i=0;i<pp->np;i++){
    lab_params_free(pp->p[i]);
  }
  free(pp->p);
  pp->np = 0;
  free (pp->zkp);
  pp->zkp = NULL;
  pp->zkround = 0;
}

static void lab_prove(
  lab_proof pi, 
  statement ost, 
  witness owt, 
  const statement ist,
  const witness iwt, 
  const lab_params pp
)
{
  timing time;

  if(pp->compressed){
    timing_start(&time, "Labradoodle Prover");
    ldd_prove(pi, ost, owt, ist, iwt, pp);
    if (pack_debug_verify_enabled()) {
      fprintf(stderr,
              "[Labrador Pack 诊断] Labradoodle compile_bincnst 前关系：%s\n",
              verify(ost, owt) ? "通过" : "失败");
    }
    compile_bincnst(ost, owt);
    if (pack_debug_verify_enabled()) {
      fprintf(stderr,
              "[Labrador Pack 诊断] Labradoodle compile_bincnst 后关系：%s\n",
              verify(ost, owt) ? "通过" : "失败");
    }
  }
  else if(!pp->tail){
    timing_start(&time, "Labrador Prover");
    ldr_prove(pi, ost, owt, ist, iwt, pp);
  }
  else{
    timing_start(&time, "Labrador Tail Prover");
    ldr_tail_prove(pi, ost, owt, ist, iwt, pp);
    if (!ldr_tail_last_prove_valid()) last_prove_relation_valid = 0;
  }
  timing_end(&time);
  timing_print(&time, 2);
}

static int lab_reduce(
  statement ost, 
  const statement ist, 
  const lab_proof pi,
  const lab_params pp
)
{
  int ret;
  timing time;

  if(pp->compressed){
    double compile_started = 0.0;
    timing_start(&time, "Labradoodle Verifier");
    ldd_reduce(ost, ist, pi, pp);
    if (revac_labrador_profile_detail_enabled()) {
      compile_started = revac_labrador_profile_seconds();
    }
    compile_bincnst(ost, NULL);
    if (revac_labrador_profile_detail_enabled()) {
      fprintf(stderr,
              "[Labrador 原生细分] 编译二进制度约束：%.6f 秒\n",
              revac_labrador_profile_seconds() - compile_started);
    }
    ret = 0;
  }
  else if(!pp->tail){
    timing_start(&time, "Labrador Verifier");
    ret = ldr_reduce(ost, ist, pi, pp);
  }
  else{
    timing_start(&time, "Labrador Tail Verifier");
    ret = ldr_tail_reduce(ost, ist, pi, pp);
  }
  timing_end(&time);
  timing_print(&time, 2);
  return ret;
}

void pack_prove(
  pack_proof pi,
  const statement ist,
  const witness iwt,
  const pack_params pp
)
{
  size_t i, j, limit;
  statement tst[2];
  witness twt[2];
  timing time;
  const int profile = revac_labrador_profile_enabled();
  double lab_pre = 0.0, lab_post = 0.0, lnp = 0.0, started;

  last_prove_relation_valid = 1;

  pack_proof_init(pi, pp);

  if(pp->np == 0){
    witness_copy(pi->owt, iwt);
    return;
  }
  if(pp->np == 1){
    if (profile) started = revac_labrador_profile_seconds();
    lab_prove(pi->p[0], tst[0], pi->owt, ist, iwt, pp->p[0]);
    pack_check_prove_relation(lab_round_name(pp->p[0]), 0, tst[0], pi->owt);
    if (profile) profile_lab_round("Prove", 0, pp->p[0],
                                   revac_labrador_profile_seconds() - started);
    if (profile) {
      fprintf(stderr, "[Labrador 原生计时] LaBRADOR 单轮：%.6f 秒\n",
              revac_labrador_profile_seconds() - started);
    }
    statement_free(tst[0]);
    return;
  }

  if (profile) started = revac_labrador_profile_seconds();
  lab_prove(pi->p[0], tst[0], twt[0], ist, iwt, pp->p[0]);
  pack_check_prove_relation(lab_round_name(pp->p[0]), 0, tst[0], twt[0]);
  if (profile) {
    const double elapsed = revac_labrador_profile_seconds() - started;
    lab_pre += elapsed;
    profile_lab_round("Prove", 0, pp->p[0], elapsed);
  }

  j = 0;
  limit = (pp->zkp != NULL) ? pp->zkround : pp->np - 1;
  for(i=1;i<limit;i++){
    if (profile) started = revac_labrador_profile_seconds();
    lab_prove(pi->p[i], tst[j^1], twt[j^1], tst[j], twt[j], pp->p[i]);
    pack_check_prove_relation(lab_round_name(pp->p[i]), i, tst[j^1],
                              twt[j^1]);
    if (profile) {
      const double elapsed = revac_labrador_profile_seconds() - started;
      lab_pre += elapsed;
      profile_lab_round("Prove", i, pp->p[i], elapsed);
    }
    statement_free(tst[j]);
    witness_free(twt[j]);
    j ^= 1;
  }
  if(pp->zkp != NULL){
    timing_start(&time,"LNP Prover");
    if (profile) started = revac_labrador_profile_seconds();
    lnp_prove (pi->zkp[0], tst[j^1], twt[j^1], tst[j], twt[j], pp->zkp[0]);
    if (profile) {
      lnp += revac_labrador_profile_seconds() - started;
      if (revac_labrador_profile_detail_enabled()) {
        fprintf(stderr,
                "[Labrador 原生细分] Prove 第%zu轮 LNP：%.6f 秒\n",
                i + 1, lnp);
      }
    }
    timing_end(&time);
    timing_print(&time, 2);

    compile_bincnst(tst[j^1], twt[j^1]);
    pack_check_prove_relation("LNP", i, tst[j^1], twt[j^1]);

    memset (pi->p[i], 0, sizeof(pi->p[i])); // empty lab proof
    statement_free(tst[j]);
    witness_free(twt[j]);
    j ^= 1;
    for(i=pp->zkround+1;i<pp->np-1;i++){
      if (profile) started = revac_labrador_profile_seconds();
      lab_prove(pi->p[i], tst[j^1], twt[j^1], tst[j], twt[j], pp->p[i]);
      pack_check_prove_relation(lab_round_name(pp->p[i]), i, tst[j^1],
                                twt[j^1]);
      if (profile) {
        const double elapsed = revac_labrador_profile_seconds() - started;
        lab_post += elapsed;
        profile_lab_round("Prove", i, pp->p[i], elapsed);
      }
      statement_free(tst[j]);
      witness_free(twt[j]);
      j ^= 1;
    }
    if(pp->zkround == pp->np-1){
      if (profile) {
        fprintf(stderr,
                "[Labrador 原生计时] LaBRADOR 折叠（LNP 前）：%.6f 秒；LNP 零知识阶段：%.6f 秒；LaBRADOR 折叠（LNP 后）：%.6f 秒\n",
                lab_pre, lnp, lab_post);
      }
      witness_copy(pi->owt, twt[j]);
      statement_free(tst[j]);
      witness_free(twt[j]);
      return;
    }
  }

  if (profile) started = revac_labrador_profile_seconds();
  lab_prove(pi->p[pp->np-1], tst[j^1], pi->owt, tst[j], twt[j], pp->p[pp->np-1]);
  pack_check_prove_relation(lab_round_name(pp->p[pp->np-1]), pp->np - 1,
                            tst[j^1], pi->owt);
  if (profile) {
    const double elapsed = revac_labrador_profile_seconds() - started;
    lab_post += elapsed;
    profile_lab_round("Prove", pp->np - 1, pp->p[pp->np-1], elapsed);
  }

  if (profile) {
    fprintf(stderr,
            "[Labrador 原生计时] LaBRADOR 折叠（LNP 前）：%.6f 秒；LNP 零知识阶段：%.6f 秒；LaBRADOR 折叠（LNP 后）：%.6f 秒\n",
            lab_pre, lnp, lab_post);
  }

  statement_free(tst[0]);
  statement_free(tst[1]);
  witness_free(twt[j]);
}

int pack_verify(
  const statement ist, 
  const pack_params pp, 
  const pack_proof pi
)
{
  size_t i, j, limit;
  int ret;
  statement tst[2];
  timing time;
  const int profile = revac_labrador_profile_enabled();
  double lab_pre = 0.0, lab_post = 0.0, lnp = 0.0, final = 0.0;
  double started = 0.0;

  if(pp->np == 0){
    if(verify(ist, pi->owt)){
      return 1;
    }
    else{
      fprintf(stderr, "pack_verify failed (verification of final witness)\n");
      return 0;
    }
  }

  if (profile) started = revac_labrador_profile_seconds();
  ret = lab_reduce(tst[0], ist, pi->p[0], pp->p[0]);
  if (profile) {
    const double elapsed = revac_labrador_profile_seconds() - started;
    lab_pre += elapsed;
    profile_lab_round("Verify", 0, pp->p[0], elapsed);
  }

  if(ret){
    fprintf(stderr, "pack_verify failed (reduction round 1)\n");
    return 0;
  }

  j = 0;
  limit = (pp->zkp != NULL) ? pp->zkround : pp->np;
  for(i=1;i<limit;i++){
    if (profile) started = revac_labrador_profile_seconds();
    ret = lab_reduce(tst[j^1], tst[j], pi->p[i], pp->p[i]);
    if (profile) {
      const double elapsed = revac_labrador_profile_seconds() - started;
      lab_pre += elapsed;
      profile_lab_round("Verify", i, pp->p[i], elapsed);
    }
    statement_free(tst[j]);
    if(ret){
      fprintf(stderr, "pack_verify failed (reduction round %zu)\n", i+1);
      return 0;
    }
    j ^= 1;
  }
  if(pp->zkp != NULL){
    timing_start(&time, "LNP Verifier");
    if (profile) started = revac_labrador_profile_seconds();
    ret = lnp_reduce(tst[j^1], tst[j], pi->zkp[0], pp->zkp[0]);
    if (profile) {
      lnp = revac_labrador_profile_seconds() - started;
      if (revac_labrador_profile_detail_enabled()) {
        fprintf(stderr,
                "[Labrador 原生细分] Verify 第%zu轮 LNP：%.6f 秒\n",
                i + 1, lnp);
      }
    }
    timing_end(&time);
    timing_print(&time, 2);

    compile_bincnst(tst[j^1], NULL);

    statement_free(tst[j]);
    if(ret){
      fprintf(stderr, "pack_verify failed (reduction round %zu)\n", i+1);
      return 0;
    }
    j ^= 1;
    for(i=pp->zkround+1;i<pp->np;i++){
      if (profile) started = revac_labrador_profile_seconds();
      ret = lab_reduce(tst[j^1], tst[j], pi->p[i], pp->p[i]);
      if (profile) {
        const double elapsed = revac_labrador_profile_seconds() - started;
        lab_post += elapsed;
        profile_lab_round("Verify", i, pp->p[i], elapsed);
      }
      statement_free(tst[j]);
      if(ret){
        fprintf(stderr, "pack_verify failed (reduction round %zu)\n", i+1);
        return 0;
      }
      j ^= 1;
    }
  }

  if (profile) started = revac_labrador_profile_seconds();
  if(verify(tst[j], pi->owt)){
    ret = 1;
  }
  else{
    fprintf(stderr, "pack_verify failed (verification of final witness)\n");
    ret = 0;
  }
  if (profile) {
    final = revac_labrador_profile_seconds() - started;
    fprintf(stderr,
            "[Labrador 原生计时] Verify 折叠（LNP 前）：%.6f 秒；LNP：%.6f 秒；折叠（LNP 后）：%.6f 秒；最终关系检查：%.6f 秒\n",
            lab_pre, lnp, lab_post, final);
  }

  statement_free(tst[j]);
  return ret;
}
