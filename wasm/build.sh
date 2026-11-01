#!/usr/bin/env bash
# Compile primer3's thal.c plus the DimerDetective wrapper to WebAssembly.
# Requires Emscripten (emcc) on PATH. Output: web/thal.js (wasm embedded).
set -euo pipefail
cd "$(dirname "$0")"

emcc -O2 \
  dd_thal.c primer3/thal.c \
  -o ../web/thal.js \
  --embed-file primer3/primer3_config@/primer3_config/ \
  -sMODULARIZE=1 \
  -sEXPORT_NAME=createThal \
  -sENVIRONMENT=web,worker,node \
  -sSINGLE_FILE=1 \
  -sALLOW_MEMORY_GROWTH=1 \
  -sEXPORTED_FUNCTIONS=_dd_init,_dd_end_stability,_dd_last_structure,_dd_last_message,_dd_last_no_structure,_dd_last_tm,_dd_last_dh,_dd_last_ds \
  -sEXPORTED_RUNTIME_METHODS=ccall,cwrap,UTF8ToString

echo "built ../web/thal.js ($(wc -c < ../web/thal.js) bytes)"
