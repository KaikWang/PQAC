cmake_minimum_required(VERSION 3.11)

project(revac_uav_portable C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)
find_package(ZLIB REQUIRED)

set(COMMON_DIR "${CMAKE_CURRENT_LIST_DIR}/../common")
set(CORE_DIR "${COMMON_DIR}/core")

file(GLOB_RECURSE CORE_SOURCES CONFIGURE_DEPENDS
  "${CORE_DIR}/src/*.c")
list(FILTER CORE_SOURCES EXCLUDE REGEX ".*/CMakeFiles/.*")
list(FILTER CORE_SOURCES EXCLUDE REGEX ".*/build/.*")

add_executable(uav_demo
  "uav_demo.c"
  "revac_uav.c"
  "${COMMON_DIR}/revocable_ac.c"
  "${COMMON_DIR}/revac_context.c"
  "${COMMON_DIR}/revac_sync.c"
  "${COMMON_DIR}/revac_labrador_bridge.c"
  "${COMMON_DIR}/revac_proof_codec.c"
  ${CORE_SOURCES})

target_include_directories(uav_demo PRIVATE
  "${CMAKE_CURRENT_LIST_DIR}"
  "${COMMON_DIR}"
  "${CORE_DIR}/include"
  "${CORE_DIR}/include/arith"
  "${CORE_DIR}/include/arith/arith_q"
  "${CORE_DIR}/include/arith/arith_qiss"
  "${CORE_DIR}/include/arith/arith_qshow"
  "${CORE_DIR}/include/arith/arith_real"
  "${CORE_DIR}/include/arith/arith_z"
  "${CORE_DIR}/flint_headers")

target_compile_options(uav_demo PRIVATE
  -Wall -Wextra -Wundef -Wshadow -Wpointer-arith
  -Wcast-align -Wstrict-prototypes -Wwrite-strings
  -Wcast-qual -Wswitch-default -Wswitch-enum -Wunreachable-code
  -march=native
  -O3)

target_link_libraries(uav_demo PRIVATE
  "${CORE_DIR}/lib/libflint.so"
  ZLIB::ZLIB
  m)

set_target_properties(uav_demo PROPERTIES
  BUILD_RPATH "${CORE_DIR}/lib")
