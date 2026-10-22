cmake_minimum_required(VERSION 3.15)
project(sample_size_calculator VERSION 0.1.0 LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(SSC_BUILD_TESTS "Build the C++ test executable" ON)
option(SSC_BUILD_CLI "Build the ssc_cli command-line tool" ON)
option(SSC_NATIVE_ARCH "Optimize for the build machine's CPU (-march=native); not for redistribution" OFF)
option(SSC_WASM_THREADS "WebAssembly: use threads (the page must be cross-origin isolated)" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
if(APPLE AND NOT CMAKE_OSX_DEPLOYMENT_TARGET)
  set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15" CACHE STRING "Minimum macOS version" FORCE)
endif()

find_package(Threads REQUIRED)
if(EMSCRIPTEN)
  # the engine reports errors with C++ exceptions: WebAssembly's own (every current browser and Node have them);
  # every object (the engine's and the tests') must agree, so these apply to the whole build
  add_compile_options(-fwasm-exceptions)
  add_link_options(-fwasm-exceptions)
  if(SSC_WASM_THREADS)
    add_compile_options(-pthread)
    add_link_options(-pthread)
  endif()
endif()

# ---------------------------------------------------------------------------
# Core sources (compiled once, linked into the shared library, CLI and tests)
# ---------------------------------------------------------------------------
file(GLOB_RECURSE SSC_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/core/src/*.cpp")

add_library(ssc_obj OBJECT ${SSC_SOURCES})
target_include_directories(ssc_obj PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/core/include"
                                          "${CMAKE_CURRENT_SOURCE_DIR}/core/src")
target_compile_definitions(ssc_obj PRIVATE SSC_BUILDING)
set_target_properties(ssc_obj PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON)

if(MSVC)
  # /utf-8: sources contain Greek letters in labels; /fp:precise keeps IEEE semantics.
  target_compile_options(ssc_obj PRIVATE /utf-8 /W3 /permissive- /fp:precise /bigobj /O2)
else()
  # Note: never use -ffast-math here (NaN/Inf tests and careful summations rely on IEEE semantics).
  # -fopenmp-simd: honor `#pragma omp simd` (loop vectorization only; no OpenMP runtime)
  target_compile_options(ssc_obj PRIVATE -Wall -Wextra -Wno-unused-parameter -fno-math-errno -fopenmp-simd)
  if(SSC_NATIVE_ARCH)
    target_compile_options(ssc_obj PRIVATE -march=native)
  endif()
endif()
target_link_libraries(ssc_obj PUBLIC Threads::Threads)

function(ssc_static_runtime target)
  if(MSVC)
    set_property(TARGET ${target} PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
  elseif(MINGW)
    target_link_options(${target} PRIVATE -static)
  endif()
endfunction()
if(MSVC)
  set_property(TARGET ssc_obj PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

if(EMSCRIPTEN)
  # ---------------------------------------------------------------------------
  # WebAssembly: the engine for a web page (web/), its C API called from JavaScript
  # ---------------------------------------------------------------------------
  set(SSC_WASM_EXPORTS
    _ssc_version _ssc_hardware_threads _ssc_catalog_json _ssc_procedure_json _ssc_search_json _ssc_run_json
    _ssc_free _ssc_token_new _ssc_token_cancel _ssc_token_progress _ssc_token_free _malloc _free)
  string(REPLACE ";" "," SSC_WASM_EXPORTS "${SSC_WASM_EXPORTS}")
  set(SSC_WASM_LINK -O3 --no-entry -sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=createEngine
      -sENVIRONMENT=web,worker,node -sALLOW_MEMORY_GROWTH=1 -sSTACK_SIZE=4MB
      "-sEXPORTED_FUNCTIONS=${SSC_WASM_EXPORTS}"
      -sEXPORTED_RUNTIME_METHODS=UTF8ToString,stringToUTF8,lengthBytesUTF8,HEAP8)
  if(SSC_WASM_THREADS)
    # growing shared memory slows JavaScript's access to it; the page only copies short strings in and out.
    # It starts at 128 MB, more than the largest answer needs (20,000 scenarios: 100 MB): Safari's engine (WebKit)
    # stopped the multi-core engine the first time its shared memory had to grow.  Operating systems back memory
    # only once it is used, so the unused part costs little.
    list(APPEND SSC_WASM_LINK -sPTHREAD_POOL_SIZE=navigator.hardwareConcurrency
         -sDEFAULT_PTHREAD_STACK_SIZE=2MB -sINITIAL_MEMORY=128MB -Wno-pthreads-mem-growth)
  endif()
  add_executable(ssc_web $<TARGET_OBJECTS:ssc_obj>)
  target_link_options(ssc_web PRIVATE ${SSC_WASM_LINK})
  set_target_properties(ssc_web PROPERTIES SUFFIX ".mjs" OUTPUT_NAME $<IF:$<BOOL:${SSC_WASM_THREADS}>,ssc_web_mt,ssc_web>)
  if(SSC_BUILD_TESTS AND SSC_WASM_THREADS)
    # the C++ tests (some test the thread pool), run under Node with the real file system for their fixtures:
    # node build-wasm-on/ssc_tests.js
    file(GLOB SSC_TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/core/tests/*.cpp")
    add_executable(ssc_tests ${SSC_TEST_SOURCES} $<TARGET_OBJECTS:ssc_obj>)
    target_include_directories(ssc_tests PRIVATE core/include core/src)
    target_compile_definitions(ssc_tests PRIVATE
      SSC_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/core/tests/fixtures")
    target_link_options(ssc_tests PRIVATE -O3 -sNODERAWFS=1 -sENVIRONMENT=node -sALLOW_MEMORY_GROWTH=1
                        -sSTACK_SIZE=8MB)
  endif()
  return()
endif()

# Shared library with a plain C ABI (loaded from Python through ctypes, so a
# single binary serves every Python 3 version).
add_library(ssc_core SHARED $<TARGET_OBJECTS:ssc_obj>)
target_include_directories(ssc_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/core/include")
target_link_libraries(ssc_core PRIVATE Threads::Threads)
set_target_properties(ssc_core PROPERTIES
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
  PREFIX "lib"            # libssc_core.{so,dylib,dll} on every platform
  WINDOWS_EXPORT_ALL_SYMBOLS OFF)
ssc_static_runtime(ssc_core)

if(SSC_BUILD_CLI)
  add_executable(ssc_cli core/cli/main.cpp $<TARGET_OBJECTS:ssc_obj>)
  target_include_directories(ssc_cli PRIVATE core/include core/src)
  # the core is linked in statically: no dllimport for the C API declarations
  target_compile_definitions(ssc_cli PRIVATE SSC_STATIC)
  target_link_libraries(ssc_cli PRIVATE Threads::Threads)
  if(MSVC)
    target_compile_options(ssc_cli PRIVATE /utf-8)
  endif()
  ssc_static_runtime(ssc_cli)
endif()

if(SSC_BUILD_TESTS)
  enable_testing()
  file(GLOB SSC_TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/core/tests/*.cpp")
  add_executable(ssc_tests ${SSC_TEST_SOURCES} $<TARGET_OBJECTS:ssc_obj>)
  target_include_directories(ssc_tests PRIVATE core/include core/src)
  target_link_libraries(ssc_tests PRIVATE Threads::Threads)
  target_compile_definitions(ssc_tests PRIVATE
    SSC_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/core/tests/fixtures")
  if(MSVC)
    target_compile_options(ssc_tests PRIVATE /utf-8)
  endif()
  ssc_static_runtime(ssc_tests)
  add_test(NAME ssc_tests COMMAND ssc_tests)
endif()

# ---------------------------------------------------------------------------
# Installation (used by scikit-build-core when building the Python wheel)
# ---------------------------------------------------------------------------
if(SKBUILD)
  install(TARGETS ssc_core
          LIBRARY DESTINATION ssc/_native
          RUNTIME DESTINATION ssc/_native)
else()
  install(TARGETS ssc_core LIBRARY DESTINATION lib RUNTIME DESTINATION bin ARCHIVE DESTINATION lib)
  install(FILES core/include/ssc/ssc.h DESTINATION include/ssc)
endif()
