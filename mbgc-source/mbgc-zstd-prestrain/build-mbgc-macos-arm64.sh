#!/bin/bash
# build-mbgc-macos-arm64.sh — 在 Apple Silicon 上原生构建集成 zstd 的 MBGC 内核
# 与 Linux 版 build-mbgc-zstd.sh 的差异：
#   1. macOS 无静态 libc，去掉 -static；
#   2. libdeflate 使用 lib/arm/cpu_features.c（而非 x86）；
#   3. 无 x86 专用 libaelf64.a，arm64 不需要汇编加速库；
#   4. OpenMP 使用 libomp（此处指向本机已有 conda 运行库），并把 rpath 设为 @executable_path；
#   5. zstd 静态库需在 macOS 侧重新 make libzstd.a。
# 用法: bash build-mbgc-macos-arm64.sh
set -e
cd "$(dirname "$0")/src"

CC=${CC:-clang}
CXX=${CXX:-clang++}
ZSTD=../vendor/zstd-1.5.6/lib
OMP_PREFIX=${OMP_PREFIX:-/opt/homebrew/opt/libomp}
BUILD=/tmp/mbgc-macos-build
mkdir -p "$BUILD"

if [ ! -f "$ZSTD/libzstd.a" ]; then
  echo "!! 缺少 $ZSTD/libzstd.a，请先执行： make -C ../vendor/zstd-1.5.6/lib libzstd.a CC=clang"
  exit 1
fi
if [ ! -f "$OMP_PREFIX/include/omp.h" ]; then
  echo "!! 缺少 OpenMP 头文件: $OMP_PREFIX/include/omp.h"
  exit 1
fi

echo "==> compiling libdeflate (C, arm64)"
$CC -O3 -c coders/libdeflate/lib/deflate_compress.c   -o "$BUILD/ld_compress.o"
$CC -O3 -c coders/libdeflate/lib/deflate_decompress.c -o "$BUILD/ld_decompress.o"
$CC -O3 -c coders/libdeflate/lib/utils.c             -o "$BUILD/ld_utils.o"
$CC -O3 -c coders/libdeflate/lib/arm/cpu_features.c  -o "$BUILD/ld_cpu.o"
$CC -O3 -c coders/libdeflate/lib/gzip_compress.c     -o "$BUILD/ld_gzcomp.o"
$CC -O3 -c coders/libdeflate/lib/gzip_decompress.c   -o "$BUILD/ld_gzdecomp.o"
$CC -O3 -c coders/libdeflate/lib/crc32.c             -o "$BUILD/ld_crc32.o"

echo "==> compiling lzma sdk (C)"
for f in LzmaEnc LzmaDec Alloc LzFindMt LzFind Ppmd7 Ppmd7Enc Ppmd7Dec Threads; do
  $CC -O3 -c coders/lzma/$f.c -o "$BUILD/lz_$f.o"
done

echo "==> compiling xxhash (C)"
$CC -O3 -c utils/hashes/xxhash.c -o "$BUILD/hash_xx.o"

echo "==> compiling C++ sources + linking (with zstd, libomp)"
$CXX -O3 -std=c++17 -Xpreprocessor -fopenmp -I"$OMP_PREFIX/include" \
  -I. -Imbgccoder -Icoders -Imatching -Iutils -Ilibs -Icoders/libdeflate -I$ZSTD \
  -o ../mbgc-macos-arm64 main.cpp \
  mbgccoder/MBGC_Decoder.cpp mbgccoder/MBGC_Encoder.cpp \
  coders/CodersLib.cpp coders/ContextAwareMismatchesCoder.cpp coders/LzmaCoder.cpp coders/PpmdCoder.cpp \
  coders/PropsLibrary.cpp coders/SymbolsPackingFacility.cpp coders/VarLenDNACoder.cpp \
  matching/MultipleGenomeMatchingProcessor.cpp matching/SimpleSequenceMatcher.cpp \
  matching/SlidingWindowSparseEMMatcher.cpp matching/TextMatchers.cpp matching/input_with_libdeflate_wrapper.cpp \
  matching/copmem/CopMEMMatcher.cpp \
  utils/helper.cpp utils/hashes/city.cpp \
  "$BUILD/ld_compress.o" "$BUILD/ld_decompress.o" "$BUILD/ld_utils.o" "$BUILD/ld_cpu.o" \
  "$BUILD/ld_gzcomp.o" "$BUILD/ld_gzdecomp.o" "$BUILD/ld_crc32.o" \
  "$BUILD/lz_LzmaEnc.o" "$BUILD/lz_LzmaDec.o" "$BUILD/lz_Alloc.o" "$BUILD/lz_LzFindMt.o" "$BUILD/lz_LzFind.o" \
  "$BUILD/lz_Ppmd7.o" "$BUILD/lz_Ppmd7Enc.o" "$BUILD/lz_Ppmd7Dec.o" "$BUILD/lz_Threads.o" \
  "$BUILD/hash_xx.o" \
  $ZSTD/libzstd.a -L"$OMP_PREFIX/lib" -lomp -Wl,-rpath,@executable_path -Wl,-rpath,/usr/lib -lpthread

echo "==> done: ../mbgc-macos-arm64"
ls -la ../mbgc-macos-arm64
