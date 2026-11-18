#!/bin/bash
# build-mbgc-zstd.sh — 集成 zstd 的 MBGC 内核构建（完全静态，默认 -march=native）
# 用法: bash build-mbgc-zstd.sh
#   MARCH=           可关闭 -march=native（通用 x86-64）
#   BUILD_DYNAMIC=1  改为动态链接（默认完全静态）
# 产物: ../mbgc.zstd
# 依赖: gcc, g++；静态版需要静态 glibc（libc6-dev/build-essential）
#       zstd 静态库来自 ../vendor/zstd-1.5.6/lib/libzstd.a
set -e
cd "$(dirname "$0")/src"

CXX=${CXX:-g++}
CC=${CC:-gcc}
MARCH=${MARCH:--march=native}
ZSTD=../vendor/zstd-1.5.6/lib

if [ ! -f "$ZSTD/libzstd.a" ]; then
  echo "!! vendor zstd 静态库缺失: $ZSTD/libzstd.a"
  echo "   请先在 vendor/zstd-1.5.6/lib 下执行: make libzstd.a"
  exit 1
fi

LINK_FLAGS=""
if [ "$BUILD_DYNAMIC" != "1" ]; then
  echo "==> STATIC build (default)"
  LINK_FLAGS="-static"
fi

echo "==> compiling libdeflate (C)"
$CC -O3 $MARCH -c coders/libdeflate/lib/deflate_compress.c   -o /tmp/zstd_ld_compress.o
$CC -O3 $MARCH -c coders/libdeflate/lib/deflate_decompress.c -o /tmp/zstd_ld_decompress.o
$CC -O3 $MARCH -c coders/libdeflate/lib/utils.c             -o /tmp/zstd_ld_utils.o
$CC -O3 $MARCH -c coders/libdeflate/lib/x86/cpu_features.c  -o /tmp/zstd_ld_cpu.o
$CC -O3 $MARCH -c coders/libdeflate/lib/gzip_compress.c     -o /tmp/zstd_ld_gzcomp.o
$CC -O3 $MARCH -c coders/libdeflate/lib/gzip_decompress.c   -o /tmp/zstd_ld_gzdecomp.o
$CC -O3 $MARCH -c coders/libdeflate/lib/crc32.c             -o /tmp/zstd_ld_crc32.o

echo "==> compiling lzma sdk (C)"
for f in LzmaEnc LzmaDec Alloc LzFindMt LzFind Ppmd7 Ppmd7Enc Ppmd7Dec Threads; do
  $CC -O3 $MARCH -c coders/lzma/$f.c -o /tmp/zstd_lz_$f.o
done

echo "==> compiling xxhash (C)"
$CC -O3 $MARCH -c utils/hashes/xxhash.c -o /tmp/zstd_hash_xx.o

echo "==> compiling C++ sources + linking (with zstd)"
$CXX -O3 $MARCH -fopenmp -std=c++17 $LINK_FLAGS -I. -Imbgccoder -Icoders -Imatching -Iutils -Ilibs -Icoders/libdeflate -I$ZSTD \
  -o ../mbgc.zstd main.cpp \
  mbgccoder/MBGC_Decoder.cpp mbgccoder/MBGC_Encoder.cpp \
  coders/CodersLib.cpp coders/ContextAwareMismatchesCoder.cpp coders/LzmaCoder.cpp coders/PpmdCoder.cpp \
  coders/PropsLibrary.cpp coders/SymbolsPackingFacility.cpp coders/VarLenDNACoder.cpp \
  matching/MultipleGenomeMatchingProcessor.cpp matching/SimpleSequenceMatcher.cpp \
  matching/SlidingWindowSparseEMMatcher.cpp matching/TextMatchers.cpp matching/input_with_libdeflate_wrapper.cpp \
  matching/copmem/CopMEMMatcher.cpp \
  utils/helper.cpp utils/hashes/city.cpp \
  /tmp/zstd_ld_compress.o /tmp/zstd_ld_decompress.o /tmp/zstd_ld_utils.o /tmp/zstd_ld_cpu.o \
  /tmp/zstd_ld_gzcomp.o /tmp/zstd_ld_gzdecomp.o /tmp/zstd_ld_crc32.o \
  /tmp/zstd_lz_LzmaEnc.o /tmp/zstd_lz_LzmaDec.o /tmp/zstd_lz_Alloc.o /tmp/zstd_lz_LzFindMt.o /tmp/zstd_lz_LzFind.o \
  /tmp/zstd_lz_Ppmd7.o /tmp/zstd_lz_Ppmd7Enc.o /tmp/zstd_lz_Ppmd7Dec.o /tmp/zstd_lz_Threads.o \
  /tmp/zstd_hash_xx.o \
  libs/libaelf64.a $ZSTD/libzstd.a -lpthread

echo "==> done: ../mbgc.zstd"
ls -la ../mbgc.zstd
