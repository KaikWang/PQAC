cmake_minimum_required(VERSION 3.11)

project(revac_entities_portable C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)
set(REVAC_LAZER_PYTHON "/home/wkk/miniconda3/envs/sage/bin/python"
  CACHE FILEPATH "Python interpreter that can import /home/wkk/pq_ac/lazer/python")
set(REVAC_SDE "/home/wkk/tools/sde-external-10.8.0-2026-03-15-lin/sde64"
  CACHE FILEPATH "Intel SDE executable for AVX-512 Labrador tests")
find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS C)
find_package(Threads REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
find_library(REVAC_MPFR_LIBRARY mpfr)
find_library(REVAC_GMP_LIBRARY gmp)
if(NOT REVAC_MPFR_LIBRARY OR NOT REVAC_GMP_LIBRARY)
  message(FATAL_ERROR "q96 requires MPFR and GMP")
endif()

set(ENTITIES_DIR "${CMAKE_CURRENT_LIST_DIR}")
set(REVAC_PROJECT_ROOT "${ENTITIES_DIR}/../..")
set(REVAC_PARAMETER_MANIFEST
  "${REVAC_PROJECT_ROOT}/params/revac128-candidate-v1.json"
  CACHE FILEPATH "Active revocable-credential parameter manifest")
set(COMMON_DIR "${ENTITIES_DIR}/common")
set(CORE_DIR "${COMMON_DIR}/core")

add_custom_target(revac_params_check
  COMMAND "${Python3_EXECUTABLE}"
          "${REVAC_PROJECT_ROOT}/scripts/params/generate_revac_profile.py"
          --manifest "${REVAC_PARAMETER_MANIFEST}" --check --quiet
  COMMAND "${Python3_EXECUTABLE}"
          "${REVAC_PROJECT_ROOT}/scripts/params/check_q88_candidate.py"
  WORKING_DIRECTORY "${REVAC_PROJECT_ROOT}"
  COMMENT "Check generated revocable-credential parameter profile"
  VERBATIM)

file(GLOB_RECURSE CORE_SOURCES CONFIGURE_DEPENDS
  "${CORE_DIR}/src/*.c")
list(FILTER CORE_SOURCES EXCLUDE REGEX ".*/CMakeFiles/.*")
list(FILTER CORE_SOURCES EXCLUDE REGEX ".*/build/.*")

add_library(revac_entities_core STATIC
  "${COMMON_DIR}/revocable_ac.c"
  "${COMMON_DIR}/revac_context.c"
  "${COMMON_DIR}/revac_sync.c"
  "${COMMON_DIR}/revac_labrador_bridge.c"
  "${COMMON_DIR}/revac_proof_codec.c"
  "${COMMON_DIR}/revac_cross_ring.c"
  "${COMMON_DIR}/revac_q96_link_wire.c"
  "${ENTITIES_DIR}/ta/revac_ta.c"
  "${ENTITIES_DIR}/uav/revac_uav.c"
  "${ENTITIES_DIR}/gcs/revac_gcs.c"
  ${CORE_SOURCES})
add_dependencies(revac_entities_core revac_params_check)

target_include_directories(revac_entities_core PUBLIC
  "${COMMON_DIR}"
  "${ENTITIES_DIR}/ta"
  "${ENTITIES_DIR}/uav"
  "${ENTITIES_DIR}/gcs"
  "${CORE_DIR}/include"
  "${CORE_DIR}/include/arith"
  "${CORE_DIR}/include/arith/arith_q"
  "${CORE_DIR}/include/arith/arith_qiss"
  "${CORE_DIR}/include/arith/arith_qshow"
  "${CORE_DIR}/include/arith/arith_real"
  "${CORE_DIR}/include/arith/arith_z"
  "${CORE_DIR}/flint_headers")

target_compile_options(revac_entities_core PRIVATE
  -Wall -Wextra -Wundef -Wshadow -Wpointer-arith
  -Wcast-align -Wstrict-prototypes -Wwrite-strings
  -Wcast-qual -Wswitch-default -Wswitch-enum -Wunreachable-code
  -march=native
  -O3)

# Resolve bridge executables and scripts when configuring each deployment.
# This keeps build artifacts portable between WSL and AutoDL without relying
# on developer-specific absolute paths compiled into the C binary.
target_compile_definitions(revac_entities_core PRIVATE
  REVAC_LABRADOR_DEFAULT_PYTHON="${REVAC_LAZER_PYTHON}"
  REVAC_LABRADOR_DEFAULT_SDE="${REVAC_SDE}"
  REVAC_LABRADOR_DEFAULT_SCRIPT="${ENTITIES_DIR}/zk_lazer/revac_acc_memver_import_labrador.py"
  REVAC_LABRADOR_DEFAULT_WORKER_SCRIPT="${ENTITIES_DIR}/zk_lazer/revac_acc_memver_worker.py"
  REVAC_LABRADOR_DEFAULT_Q96_SCRIPT="${ENTITIES_DIR}/zk_lazer/revac_acc96_memver_labrador.py"
  REVAC_LABRADOR_DEFAULT_Q96_WORKER_SCRIPT="${ENTITIES_DIR}/zk_lazer/revac_acc96_memver_worker.py")

target_link_libraries(revac_entities_core PUBLIC
  "${CORE_DIR}/lib/libflint.so"
  OpenMP::OpenMP_C
  Threads::Threads
  ZLIB::ZLIB
  "${REVAC_MPFR_LIBRARY}"
  "${REVAC_GMP_LIBRARY}"
  m)

set_target_properties(revac_entities_core PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(uav_demo
  "${ENTITIES_DIR}/uav/uav_demo.c")
target_link_libraries(uav_demo PRIVATE revac_entities_core)
set_target_properties(uav_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revocable_ac_demo
  "${ENTITIES_DIR}/test/revocable_ac_demo.c")
target_link_libraries(revocable_ac_demo PRIVATE revac_entities_core)
set_target_properties(revocable_ac_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revocable_ac_q19_legacy_demo
  "${ENTITIES_DIR}/test/revocable_ac_q19_legacy_demo.c")
target_link_libraries(revocable_ac_q19_legacy_demo PRIVATE revac_entities_core)
set_target_properties(revocable_ac_q19_legacy_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_acc_memver_export
  "${ENTITIES_DIR}/tools/revac_acc_memver_export.c")
target_link_libraries(revac_acc_memver_export PRIVATE revac_entities_core)
set_target_properties(revac_acc_memver_export PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_acc_eval_bench
  "${ENTITIES_DIR}/bench/revac_acc_eval_bench.c")
target_link_libraries(revac_acc_eval_bench PRIVATE revac_entities_core)
set_target_properties(revac_acc_eval_bench PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_acc_update_bound_bench
  "${ENTITIES_DIR}/bench/revac_acc_update_bound_bench.c")
target_link_libraries(revac_acc_update_bound_bench PRIVATE revac_entities_core)
set_target_properties(revac_acc_update_bound_bench PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_acc_lift_demo
  "${ENTITIES_DIR}/test/revac_acc_lift_demo.c")
target_link_libraries(revac_acc_lift_demo PRIVATE revac_entities_core)
set_target_properties(revac_acc_lift_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_acc_parallel_demo
  "${ENTITIES_DIR}/test/revac_acc_parallel_demo.c")
target_link_libraries(revac_acc_parallel_demo PRIVATE revac_entities_core)
set_target_properties(revac_acc_parallel_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_acc_base3_correctness
  "${ENTITIES_DIR}/test/revac_acc_base3_correctness.c")
target_link_libraries(revac_acc_base3_correctness PRIVATE revac_entities_core)
set_target_properties(revac_acc_base3_correctness PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_acc_balanced_odd_base_kat
  "${ENTITIES_DIR}/test/revac_acc_balanced_odd_base_kat.c")

add_executable(revac_acc96_crt_demo
  "${ENTITIES_DIR}/test/revac_acc80_crt_demo.c"
  "${CORE_DIR}/src/acc80_crt.c")
target_include_directories(revac_acc96_crt_demo PRIVATE
  "${CORE_DIR}/include")
target_compile_options(revac_acc96_crt_demo PRIVATE
  -Wall -Wextra -Werror -march=native -O3)

add_executable(revac_acc96_eval_demo
  "${ENTITIES_DIR}/test/revac_acc80_eval_demo.c"
  "${CORE_DIR}/src/acc80_crt.c"
  "${CORE_DIR}/src/acc80_eval.c")
target_include_directories(revac_acc96_eval_demo PRIVATE
  "${CORE_DIR}/include")
target_compile_options(revac_acc96_eval_demo PRIVATE
  -Wall -Wextra -Werror -march=native -O3)
target_link_libraries(revac_acc96_eval_demo PRIVATE OpenMP::OpenMP_C m)

# Compatibility build names used by earlier q80/q96 experiment instructions.
add_custom_target(revac_acc80_crt_demo DEPENDS revac_acc96_crt_demo)
add_custom_target(revac_acc80_eval_demo DEPENDS revac_acc96_eval_demo)

add_executable(revac_acc96_sampler_params_demo
  "${ENTITIES_DIR}/test/revac_acc96_sampler_params_demo.c"
  "${CORE_DIR}/src/acc80_crt.c")
target_include_directories(revac_acc96_sampler_params_demo PRIVATE
  "${CORE_DIR}/include")
target_compile_options(revac_acc96_sampler_params_demo PRIVATE
  -Wall -Wextra -Wpedantic -O3 -march=native)
target_link_libraries(revac_acc96_sampler_params_demo PRIVATE m)

add_executable(revac_acc96_relation_demo
  "${ENTITIES_DIR}/test/revac_acc96_relation_demo.c"
  "${CORE_DIR}/src/acc80_crt.c"
  "${CORE_DIR}/src/acc80_eval.c"
  "${CORE_DIR}/src/acc96_delta4_codec.c"
  "${CORE_DIR}/src/acc96_relation.c"
  "${CORE_DIR}/src/fips202.c")
target_include_directories(revac_acc96_relation_demo PRIVATE
  "${CORE_DIR}/include")
target_compile_options(revac_acc96_relation_demo PRIVATE
  -Wall -Wextra -Werror -march=native -O3)
target_link_libraries(revac_acc96_relation_demo PRIVATE OpenMP::OpenMP_C m ZLIB::ZLIB)

add_executable(revac_acc96_proof_layout_demo
  "${ENTITIES_DIR}/test/revac_acc96_proof_layout_demo.c"
  "${CORE_DIR}/src/acc80_crt.c"
  "${CORE_DIR}/src/acc96_proof_layout.c")
target_include_directories(revac_acc96_proof_layout_demo PRIVATE
  "${CORE_DIR}/include")
target_compile_options(revac_acc96_proof_layout_demo PRIVATE
  -Wall -Wextra -Werror -march=native -O3)

add_executable(revac_acc96_trapdoor_demo
  "${ENTITIES_DIR}/test/revac_acc96_trapdoor_demo.c"
  "${CORE_DIR}/src/acc80_crt.c"
  "${CORE_DIR}/src/acc80_eval.c"
  "${CORE_DIR}/src/acc96_relation.c"
  "${CORE_DIR}/src/acc96_trapdoor.c"
  "${CORE_DIR}/src/acc96_trapdoor_cert.c"
  "${CORE_DIR}/src/acc96_perturbation_cert.c"
  "${CORE_DIR}/src/fips202.c")
target_include_directories(revac_acc96_trapdoor_demo PRIVATE
  "${CORE_DIR}/include"
  "${CORE_DIR}/flint_headers")
target_compile_options(revac_acc96_trapdoor_demo PRIVATE
  -Wall -Wextra -Werror -march=native -O3)
target_link_libraries(revac_acc96_trapdoor_demo PRIVATE
  "${CORE_DIR}/lib/libflint.so" OpenMP::OpenMP_C m)
set_target_properties(revac_acc96_trapdoor_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_acc96_lifecycle_demo
  "${ENTITIES_DIR}/test/revac_acc96_lifecycle_demo.c"
  "${CORE_DIR}/src/acc80_crt.c"
  "${CORE_DIR}/src/acc80_eval.c"
  "${CORE_DIR}/src/acc96_relation.c"
  "${CORE_DIR}/src/acc96_trapdoor.c"
  "${CORE_DIR}/src/acc96_accumulator.c"
  "${CORE_DIR}/src/acc96_sampleleft_kat.c"
  "${CORE_DIR}/src/fips202.c")
target_include_directories(revac_acc96_lifecycle_demo PRIVATE
  "${CORE_DIR}/include")
target_compile_options(revac_acc96_lifecycle_demo PRIVATE
  -Wall -Wextra -Werror -march=native -O3)
target_link_libraries(revac_acc96_lifecycle_demo PRIVATE OpenMP::OpenMP_C m)

add_executable(revac_acc96_klein_exact_demo
  "${ENTITIES_DIR}/test/revac_acc96_klein_exact_demo.c"
  "${CORE_DIR}/src/acc80_crt.c"
  "${CORE_DIR}/src/acc96_klein_exact.c")
target_include_directories(revac_acc96_klein_exact_demo PRIVATE
  "${CORE_DIR}/include")
target_compile_options(revac_acc96_klein_exact_demo PRIVATE
  -Wall -Wextra -Werror -march=native -O3)
target_link_libraries(revac_acc96_klein_exact_demo PRIVATE
  "${CORE_DIR}/lib/libflint.so" gmp m)
set_target_properties(revac_acc96_klein_exact_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_q96_entities_demo
  "${ENTITIES_DIR}/test/revocable_ac_demo.c")
target_link_libraries(revac_q96_entities_demo PRIVATE revac_entities_core)
set_target_properties(revac_q96_entities_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_cross_ring_demo
  "${ENTITIES_DIR}/test/revac_cross_ring_demo.c")
target_link_libraries(revac_cross_ring_demo PRIVATE revac_entities_core)
set_target_properties(revac_cross_ring_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_acc96_gaussian_mpfr_demo
  "${ENTITIES_DIR}/test/revac_acc96_gaussian_mpfr_demo.c")
target_link_libraries(revac_acc96_gaussian_mpfr_demo PRIVATE revac_entities_core)
set_target_properties(revac_acc96_gaussian_mpfr_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_acc96_samplefz_demo
  "${ENTITIES_DIR}/test/revac_acc96_samplefz_demo.c")
target_link_libraries(revac_acc96_samplefz_demo PRIVATE revac_entities_core)
set_target_properties(revac_acc96_samplefz_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_acc96_sampleleft_demo
  "${ENTITIES_DIR}/test/revac_acc96_sampleleft_demo.c")
target_link_libraries(revac_acc96_sampleleft_demo PRIVATE revac_entities_core)
set_target_properties(revac_acc96_sampleleft_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_acc96_fresh_gate_demo
  "${ENTITIES_DIR}/test/revac_acc96_fresh_gate_demo.c")
target_link_libraries(revac_acc96_fresh_gate_demo PRIVATE revac_entities_core)
set_target_properties(revac_acc96_fresh_gate_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_acc96_fresh_sync_demo
  "${ENTITIES_DIR}/test/revac_acc96_fresh_sync_demo.c")
target_link_libraries(revac_acc96_fresh_sync_demo PRIVATE revac_entities_core)
set_target_properties(revac_acc96_fresh_sync_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_context_demo
  "${ENTITIES_DIR}/test/revac_context_demo.c")
target_link_libraries(revac_context_demo PRIVATE revac_entities_core)
set_target_properties(revac_context_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_acc_improvements_ablation
  "${ENTITIES_DIR}/test/revac_acc_improvements_ablation.c")
target_link_libraries(revac_acc_improvements_ablation PRIVATE revac_entities_core)
set_target_properties(revac_acc_improvements_ablation PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_executable(revac_proof_codec_demo
  "${ENTITIES_DIR}/test/revac_proof_codec_demo.c")
target_link_libraries(revac_proof_codec_demo PRIVATE revac_entities_core)
set_target_properties(revac_proof_codec_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")

add_custom_target(revac_same_x_zk_demo
  COMMAND "${REVAC_LAZER_PYTHON}" "${ENTITIES_DIR}/zk_lazer/revac_same_x_demo.py"
  WORKING_DIRECTORY "${ENTITIES_DIR}/zk_lazer"
  COMMENT "Run the LaZer same-x prototype for revocable anonymous show")

add_custom_target(revac_memver_quadratic_zk_demo
  COMMAND "${REVAC_SDE}" -icx -- "${REVAC_LAZER_PYTHON}" "${ENTITIES_DIR}/zk_lazer/revac_memver_quadratic_demo.py" --primesize 38
  WORKING_DIRECTORY "${ENTITIES_DIR}/zk_lazer"
  COMMENT "Run the Labrador quadratic prototype for revocable anonymous show")

enable_testing()
add_test(NAME revac_params_check
  COMMAND "${Python3_EXECUTABLE}"
          "${REVAC_PROJECT_ROOT}/scripts/params/generate_revac_profile.py"
          --manifest "${REVAC_PARAMETER_MANIFEST}" --check)
add_test(NAME revac_q88_params_check
  COMMAND "${Python3_EXECUTABLE}"
          "${REVAC_PROJECT_ROOT}/scripts/params/check_q88_candidate.py")
add_test(NAME revocable_ac_demo COMMAND revocable_ac_demo)
add_test(NAME revocable_ac_q19_legacy_demo COMMAND revocable_ac_q19_legacy_demo)
add_test(NAME revac_proof_codec_demo COMMAND revac_proof_codec_demo)
add_test(NAME revac_context_demo COMMAND revac_context_demo)
add_test(NAME revac_acc_improvements_ablation COMMAND revac_acc_improvements_ablation)
add_test(NAME revac_acc_lift_demo COMMAND revac_acc_lift_demo)
add_test(NAME revac_acc_parallel_demo COMMAND revac_acc_parallel_demo)
add_test(NAME revac_acc_base3_correctness COMMAND revac_acc_base3_correctness)
add_test(NAME revac_acc_balanced_odd_base_kat
  COMMAND revac_acc_balanced_odd_base_kat)
add_test(NAME revac_acc96_crt_demo COMMAND revac_acc96_crt_demo)
add_test(NAME revac_acc96_eval_demo COMMAND revac_acc96_eval_demo)
add_test(NAME revac_acc96_sampler_params_demo COMMAND revac_acc96_sampler_params_demo)
add_test(NAME revac_acc96_relation_demo COMMAND revac_acc96_relation_demo)
add_test(NAME revac_acc96_proof_layout_demo COMMAND revac_acc96_proof_layout_demo)
add_test(NAME revac_acc96_trapdoor_demo COMMAND revac_acc96_trapdoor_demo)
add_test(NAME revac_acc96_lifecycle_demo COMMAND revac_acc96_lifecycle_demo)
add_test(NAME revac_cross_ring_demo COMMAND revac_cross_ring_demo)
add_test(NAME revac_acc96_klein_exact_demo COMMAND revac_acc96_klein_exact_demo)
add_test(NAME revac_q96_entities_demo COMMAND revac_q96_entities_demo)
add_test(NAME revac_acc96_gaussian_mpfr_demo
  COMMAND revac_acc96_gaussian_mpfr_demo)
add_test(NAME revac_acc96_samplefz_demo COMMAND revac_acc96_samplefz_demo)
add_test(NAME revac_acc96_sampleleft_demo COMMAND revac_acc96_sampleleft_demo)
add_test(NAME revac_acc96_fresh_gate_demo COMMAND revac_acc96_fresh_gate_demo)
add_test(NAME revac_acc96_fresh_sync_demo COMMAND revac_acc96_fresh_sync_demo)
add_test(NAME revac_same_x_zk_demo
  COMMAND "${REVAC_LAZER_PYTHON}" "${ENTITIES_DIR}/zk_lazer/revac_same_x_demo.py")
add_test(NAME revac_memver_quadratic_zk_demo
  COMMAND "${REVAC_SDE}" -icx -- "${REVAC_LAZER_PYTHON}" "${ENTITIES_DIR}/zk_lazer/revac_memver_quadratic_demo.py" --primesize 38)
