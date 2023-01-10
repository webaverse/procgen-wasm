mkdir -p bin
echo 'building...'
# m = 64*1024; s = 512 * 1024 * 1024; Math.floor(s/m)*m;
# emcc -s NO_EXIT_RUNTIME=1 -s TOTAL_MEMORY=52428800 -D__linux__ -s ALLOW_MEMORY_GROWTH=0 -sWASM_WORKERS=1 -g -s ASSERTIONS=1 -fexceptions \
# -g -s ASSERTIONS=1 -fsanitize=address

NUM_THREADS=8

emcc -std=c++20 -D NUM_THREADS=$NUM_THREADS -sNO_EXIT_RUNTIME=1 -s TOTAL_MEMORY=100MB -pthread -sPTHREAD_POOL_SIZE=$NUM_THREADS -sPTHREAD_POOL_SIZE_STRICT=$NUM_THREADS -s ALLOW_MEMORY_GROWTH=0 -sLLD_REPORT_UNDEFINED -O3 \
  binding.cc \
  procgen.cc \
  instance.cc \
  generation/heightfield-generator.cc \
  generation/noise.cc \
  generation/glsl.cc \
  polygonization/polygonizer.cc \
  polygonization/mesh.cc \
  task/octree.cc \
  task/result.cc \
  task/task.cc \
  task/sync.cc \
  task/promise.cc \
  task/tracker.cc \
  libs/FastNoise.cpp \
  libs/vectorMath.cc \
  libs/vector.cc \
  libs/MurmurHash3.cpp \
  libs/Worley.cpp \
  libs/NoiseTools.cpp \
  libs/NoiseBase.cpp \
  utils/util.cc \
  -I. \
  -o bin/pg.js

  #sed -Ei 's/var Module=typeof Module!="undefined"\?Module:\{\};/var Module = globalThis.Module??{};/g' bin/pg.js
  # sed -Ei 's/var asm=createWasm\(\);/asmLibraryArg.__cxa_atexit=()=>{};var asm=createWasm();/g' bin/pg.js
  sed -Ei 's/importScripts\(e.data.urlOrBlob\)/importScripts(e.data.urlOrBlob.replace(\/pg-worker\\.js.*$\/, "pg.js"))/g' bin/pg.worker.js
  echo 'let accept, reject;const p = new Promise((a, r) => {accept = a;  reject = r;});Module.postRun = () => {accept();};Module.waitForLoad = () => p;' >> bin/pg.js
  sed -Ei 's/scriptDirectory\+path/"\/"+path/g' bin/pg.js
  cp bin/pg.js bin/pg.module.js
  echo 'export default Module;' >>bin/pg.module.js
  cp -R bin/* ../app/public/
echo done
