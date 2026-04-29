#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/param.h>
#include <string.h>
#ifndef __cplusplus
# include <stdbool.h>
#endif
#ifndef __APPLE__
#include <endian.h>
#endif
#if defined(unix) || defined(__unix)
# include <unistd.h>
#endif
#ifdef __has_include
#  if __has_include(<utility>)
#     include <utility>
#  endif
#  if __has_include(<memory>)
#     include <memory>
#  endif
# if __cplusplus >= 201103L
#  if __has_include(<type_traits>)
#     include <type_traits>
#  endif
# endif
#  if __has_include(<functional>)
#     include <functional>
#  endif
#  if __has_include(<string>)
#     include <string>
#  endif
#  if __cplusplus >= 201103L
#   if __has_include(<chrono>)
#     include <chrono>
#   endif
#  endif
#  if __has_include(<map>)
#     include <map>
#  endif
#  if __has_include(<set>)
#     include <set>
#  endif
#  if __has_include(<vector>)
#     include <vector>
#  endif
#  if __has_include(<iterator>)
#     include <iterator>
#  endif
#  if __has_include(<complex>)
#     include <complex>
#  endif
#  if __has_include(<iomanip>)
#     include <iomanip>
#  endif
#  if __cplusplus >= 201103L
#   if __has_include(<mutex>)
#     include <mutex>
#   endif
#  endif
#  if __cplusplus >= 201402L
#   if __has_include(<shared_mutex>)
#     include <shared_mutex>
#   endif
#  endif
#  if __has_include(<cstddef>)
#     include <cstddef>
#  endif
#  if __has_include(<new>)
#     include <new>
#  endif
#  if __has_include(<exception>)
#     include <exception>
#  endif
#  if __cplusplus >= 201103L
#   if __has_include(<tuple>)
#     include <tuple>
#   endif
#  endif
#  if __has_include(<optional>)
#     include <optional>
#  endif
#  if __has_include(<any>)
#     include <any>
#  endif
#  if __has_include(<memory_resource>)
#     include <memory_resource>
#  endif
#  if __has_include(<string_view>)
#     include <string_view>
#  endif
/*
// #  if __has_include(<execution>)
// #     include <execution>
// #  endif
*/
#  if __has_include(<forwardlist>)
#     include <forwardlist>
#  endif
#  if __has_include(<list>)
#     include <list>
#  endif
#  if __cplusplus >= 201103L
#   if __has_include(<unordered_map>)
#     include <unordered_map>
#   endif
#  endif
#  if __has_include(<numeric>)
#     include <numeric>
#  endif
#  if __has_include(<cmath>)
#     include <cmath>
#  endif
#  if __has_include(<filesystem>)
#     include <filesystem>
#  endif
#  if __cplusplus >= 201103L
#   if __has_include(<atomic>)
#     include <atomic>
#   endif
#  endif
#  if __cplusplus >= 201103L
#   if __has_include(<thread>)
#     include <thread>
#   endif
#  endif
#  if __has_include(<version>)
#     include <version>
#  endif
#  if __cplusplus >= 201103L
#   if __has_include(<atomic>)
#     include <atomic>
#   endif
#  endif
#  if __has_include(<barrier>)
#     include <barrier>
#  endif
#  if __cplusplus >= 201902L
#   if __has_include(<coroutine>)
#     include <coroutine>
#   endif
#  endif
#  if __has_include(<functional>)
#     include <functional>
#  endif
#  if __has_include(<bit>)
#     include <bit>
#  endif
#  if __cplusplus >= 201103L
#   if __has_include(<type_traits>)
#     include <type_traits>
#   endif
#  endif
#  if __has_include(<concepts>)
#     include <concepts>
#  endif
#  if __has_include(<format>)
#     include <format>
#  endif
#  if __cplusplus >= 201103L
#   if __has_include(<unordered_map>)
#     include <unordered_map>
#   endif
#   if __has_include(<unordered_set>)
#     include <unordered_set>
#   endif
#  endif
#  if __has_include(<flat_map>)
#     include <flat_map>
#  endif
#  if __has_include(<flat_set>)
#     include <flat_set>
#  endif
#  if __has_include(<latch>)
#     include <latch>
#  endif
#  if __has_include(<numbers>)
#     include <numbers>
#  endif
#  if __has_include(<new>)
#     include <new>
#  endif
#  if __cplusplus >= 201103L
#   if __has_include(<thread>)
#     include <thread>
#   endif
#  endif
#  if __cplusplus >= 201103L
#   if __has_include(<tuple>)
#     include <tuple>
#   endif
#  endif
#  if __has_include(<format>)
#     include <format>
#  endif
#  if __cplusplus >= 201902L
#   if __has_include(<execution>)
#     include <execution>
#   endif
#  endif
#  if __has_include(<utility>)
#     include <utility>
#  endif
#  if __has_include(<expected>)
#     include <expected>
#  endif
#  if __has_include(<mdspan>)
#     include <mdspan>
#  endif
#endif

#if ( defined __STDC__ && __STDC__ ) || defined __cplusplus
#define _XSTR_(s) _STR_(s) /**< Expand a macro first and then enclose in string */
#define _STR_(s) #s /**< Enclose in string without macro expansion. */
#else
#define _XSTR_(s) === /* No working conversion method known: let compilation fail. */
#define _STR_(s) "s" /**< This one also works with traditional C. */
#endif

#define SHOW(s) if ( strcmp(#s,_XSTR_(s)) != 0 ) { if ( strcmp("", _XSTR_(s)) != 0 ) printf( #s "=" _XSTR_(s) "\n" ); else printf( #s "\n" ); }

#define CPU_IS(s) if ( __builtin_cpu_is(s) ) printf("CPU is %s\n", s);
#define CPU_SUPPORTS(s) if ( __builtin_cpu_supports(s) ) printf("CPU supports %s\n", s);

int main(int argc, char **argv)
{
#if defined(unix) || defined(__unix)
printf("* System identification (using uname, /etc/release, ...; not from preprocessor):\n");
system("uname -a");
if ( access("/etc/release",R_OK) == 0 )
   system("cat /etc/release");
else if ( access("/etc/system-release",R_OK) == 0 )
   system("cat /etc/system-release");
else if ( access("/etc/redhat-release",R_OK) == 0 )
   system("cat /etc/redhat-release");
if ( access("/etc/lsb-release",R_OK) == 0 )
{
   printf("* LSB release information:\n");
   system("cat /etc/lsb-release");
}
if ( access("/proc/cpuinfo",R_OK) == 0 )
{
   printf("* CPU information:\n");
   system("cat /proc/cpuinfo | head -26");
   printf("...\n");
}
#endif

printf("\n* Preprocessor definitions:\n");
printf("\n* Example:\n*  #define TEST_ABC\n* shows up as: ");
#define TEST_ABC
SHOW(TEST_ABC);
printf("* and\n*  #define TEST_XYZ 123\n* shows up as: ");
#define TEST_XYZ 123
SHOW(TEST_XYZ);
printf("* For a macro assigned to itself:  #define TEST_TEST TEST_TEST\n* we get: ");
#define TEST_TEST TEST_TEST
SHOW(TEST_TEST);
printf("\n* (expected to get nothing)\n");
printf("* while for an undefined TEST_NONONO we get: ");
SHOW(TEST_NONONO);
printf("\n* (expected to get nothing)\n");

printf("\nCompiler:\n");
SHOW(__VERSION__);
SHOW(__GNUC__);
SHOW(__GNUC_VERSION__);
SHOW(__GNUC_MINOR__);
SHOW(__GNUC_PATCHLEVEL__);

SHOW(__MINGW32__);
SHOW(__MINGW32_MAJOR_VERSION);
SHOW(__MINGW32_MINOR_VERSION);
SHOW(__MINGW64__);
SHOW(__MINGW64_VERSION_MAJOR);
SHOW(__MINGW64_VERSION_MINOR);

SHOW(_MSC_VER);
SHOW(_MSC_FULL_VER);
SHOW(_MSC_BUILD);

SHOW(__llvm__);
SHOW(__clang__);
SHOW(__clang_major__);
SHOW(__clang_minor__);
SHOW(__clang_patchlevel__);
SHOW(__clang_version__);
SHOW(__apple_build_version__);

SHOW(__BORLANDC__);
SHOW(__CODEGEARC__);

SHOW(__DECC);
SHOW(__DECC_VER);
SHOW(__DECCXX);
SHOW(__DECCXX_VER);
SHOW(__VAXC);
SHOW(VAXC);

SHOW(__COVERITY__);

SHOW(_CRAYC);
SHOW(_RELEASE);
SHOW(_RELEASE_MINOR);

SHOW(__EDG__);
SHOW(__EDG_VERSION__);

SHOW(__xlc__);
SHOW(__xlC__);
SHOW(__xlC_ver__);
SHOW(__IBMC__);
SHOW(__IBMCPP__);
SHOW(__COMPILER_VER__);

SHOW(__ibmxl__);
SHOW(__ibmxl_vrm__);
SHOW(__ibmxl_version__);
SHOW(__ibmxl_release__);
SHOW(__ibmxl_modification__);
SHOW(__ibmxl_ptf_fix_level__);

SHOW(__ICL);
SHOW(__ICC);
SHOW(__ECC);
SHOW(__INTEL_COMPILER);
SHOW(__INTEL_COMPILER_BUILD_DATE);

SHOW(__KCC);
SHOW(__KCC_VERSION);

SHOW(__MWERKS__);
SHOW(__CWCC__);

SHOW(__SUNPRO_C);
SHOW(__SUNPRO_CC);

SHOW(__PGI);
SHOW(__PGIC__);
SHOW(__PGIC_MINOR__);
SHOW(__PGIC_PATCHLEVEL__);

SHOW(SDCC);

SHOW(__SC__);

SHOW(THINKC3);
SHOW(THINKC4);
SHOW(__TINYC__);
SHOW(__TURBOC__);
SHOW(__WATCOMC__);
SHOW(__ZTC__);

SHOW(__CC_ARM);
SHOW(__ARMCC_VERSION);
SHOW(__REVISION__);

SHOW(__BUILD__);
SHOW(__sgi);
SHOW(sgi);
SHOW(_COMPILER_VERSION);
SHOW(_SGI_COMPILER_VERSION);

SHOW(__TI_COMPILER_VERSION__);
SHOW(_TMS320C6X);

SHOW(__APPLE_CC__);
SHOW(__SUNPRO_CC);
SHOW(__SUNPRO_CC_COMPAT);

SHOW(__OPEN64__);
SHOW(__OPENCC__);
SHOW(__OPENCC_MINOR__);
SHOW(__OPENCC_PATCHLEVEL__);

printf("\nLibrary:\n");
SHOW(__GNU_LIBRARY__);
SHOW(__GLIBC__);
SHOW(__GLIBC_MINOR__);
SHOW(__GLIBC_MAJOR__);
SHOW(__GXX_ABI_VERSION);

printf("\nLanguage:\n");
SHOW(__STDC__);
SHOW(__STDC_VERSION__);
SHOW(__STDC_HOSTED__);
SHOW(ANSI_C);
SHOW(_ANSI_C);
SHOW(__STRICT_ANSI__);
SHOW(__cplusplus);
SHOW(__cplusplus_cli);
SHOW(__embedded_cplusplus);

printf("\nSelected:\n");
SHOW(__STRICT_ANSI__);
SHOW(_ISOC90_SOURCE);
SHOW(_ISOC95_SOURCE);
SHOW(_ISOC99_SOURCE);
SHOW(_ISOC11_SOURCE);
SHOW(_ISOC17_SOURCE);
SHOW(_ISOC18_SOURCE);
SHOW(_ISOC20_SOURCE);
SHOW(_ISOC2X_SOURCE);
SHOW(_ISOCXX11_SOURCE);
SHOW(_ISOCXX14_SOURCE);
SHOW(_ISOCXX17_SOURCE);
SHOW(_ISOCXX20_SOURCE);
SHOW(_ISOCXX23_SOURCE);
SHOW(__STDC_WANT_LIB_EXT2__);
SHOW(__STDC_WANT_IEC_60559_BFP_EXT__);
SHOW(__STDC_WANT_IEC_60559_FUNCS_EXT__);
SHOW(__STDC_WANT_IEC_60559_TYPES_EXT__);
SHOW(__STDC_WANT_IEC_60559_EXT__);
SHOW(__STDC_IEC_60559_TYPES__);
SHOW(__STDC_IEC_60559_BFP__);
SHOW(__STDC_IEC_60559_DFP__);
SHOW(_ANSI_SOURCE);
SHOW(_POSIX_SOURCE);
SHOW(_POSIX_C_SOURCE);
SHOW(_POSIX_C_LEVEL);
SHOW(_DARWIN_C_SOURCE);
SHOW(_DARWIN_C_LEVEL);
SHOW(_XOPEN_SOURCE);
SHOW(_XOPEN_SOURCE_EXTENDED);
SHOW(_LARGEFILE_SOURCE);
SHOW(_LARGEFILE64_SOURCE);
SHOW(_FILE_OFFSET_BITS);
SHOW(_TIME_BITS);
SHOW(_BSD_SOURCE);
SHOW(_SVID_SOURCE);
SHOW(_ATFILE_SOURCE);
SHOW(_GNU_SOURCE);
SHOW(_DEFAULT_SOURCE);
SHOW(_REENTRANT);
SHOW(_THREAD_SAFE);
SHOW(_FORTIFY_SOURCE);
SHOW(_DYNAMIC_STACK_SIZE_SOURCE);

printf("\nFeatures:\n");
SHOW(__USE_ATFILE);
SHOW(__USE_BSD);
SHOW(__USE_FILE_OFFSET64);
SHOW(__USE_FORTIFY_LEVEL);
SHOW(__USE_GNU);
SHOW(__USE_ISOC11);
SHOW(__USE_ISOC17);
SHOW(__USE_ISOC18);
SHOW(__USE_ISOC20);
SHOW(__USE_ISOC23);
SHOW(__USE_ISOC89);
SHOW(__USE_ISOC90);
SHOW(__USE_ISOC95);
SHOW(__USE_ISOC99);
SHOW(__USE_ISOC9X);
SHOW(__USE_ISOCXX11);
SHOW(__USE_ISOCXX14);
SHOW(__USE_ISOCXX17);
SHOW(__USE_ISOCXX20);
SHOW(__USE_ISOCXX23);
SHOW(__USE_ISOCXX98);
SHOW(__USE_LARGEFILE);
SHOW(__USE_LARGEFILE64);
SHOW(__USE_MISC);
SHOW(__USE_POSIX);
SHOW(__USE_POSIX199309);
SHOW(__USE_POSIX199506);
SHOW(__USE_POSIX2);
SHOW(__USE_REENTRANT);
SHOW(__USE_SVID);
SHOW(__USE_UNIX98);
SHOW(__USE_XOPEN);
SHOW(__USE_XOPEN_EXTENDED);
SHOW(__USE_XOPEN2K);
SHOW(__USE_XOPEN2K8);
SHOW(__USE_XOPEN2K8XSI);
SHOW(__USE_XOPEN2KXSI);
SHOW(__USE_XOPEN_EXTENDED);

SHOW(__GLIBC_USE);
#if defined __GLIBC_USE
SHOW(__GLIBC_USE (LIB_EXT2));
SHOW(__GLIBC_USE (IEC_60559_TYPES_EXT));
SHOW(__GLIBC_USE (IEC_60559_BFP_EXT));
SHOW(__GLIBC_USE (DEPRECATED_SCANF));
SHOW(__GLIBC_USE (DEPRECATED_GETS));
#endif
SHOW(__FAVOR_BSD);
SHOW(__LSB_VERSION__);
SHOW(__KERNEL_STRICT_NAMES);

SHOW(__STDC_MB_MIGHT_NEQ_WC__);
SHOW(__STDC_ISO_10646__);
SHOW(__STDCPP_STRICT_POINTER_SAFETY__);
SHOW(__STDCPP_THREADS__);
SHOW(__PIC__);
SHOW(__pic__);
SHOW(_MT);
SHOW(_OPENMP);

printf("\nHardware and byte order:\n");
/* See, among others: https://sourceforge.net/p/predef/wiki/Architectures/ */
SHOW(__IA32__);
SHOW(__I86__);
SHOW(__INTEL__);
SHOW(_M_I86);
SHOW(_M_IX86);
SHOW(__X86__);
SHOW(_X86_);
SHOW(__THW_INTEL__);
SHOW(__x86_64);
SHOW(__x86_64__);
SHOW(__x86_64);
SHOW(__amd64);
SHOW(__amd64__);
SHOW(_M_X64);
SHOW(_M_AMD64);
SHOW(__i386__);
SHOW(__i386);
SHOW(i386);
SHOW(__i486__);
SHOW(__i586__);
SHOW(__i686__);
SHOW(__pentium4);
SHOW(__pentium4__ );

SHOW(__alpha);
SHOW(__alpha__);
SHOW(_M_ALPHA);
SHOW(__alpha_ev4__);
SHOW(__alpha_ev5__);
SHOW(__alpha_ev6__);

SHOW(__arm);
SHOW(__arm__);
SHOW(__thumb__);
SHOW(_ARM);
SHOW(_M_ARM);
SHOW(_M_ARMT);
SHOW(__ARM_ARCH_2__);
SHOW(__ARM_ARCH_3__);
SHOW(__ARM_ARCH_3M__);
SHOW(__ARM_ARCH_4T__);
SHOW(__TARGET_ARM_4T);
SHOW(__ARM_ARCH_5__);
SHOW(__ARM_ARCH_5E__);
SHOW(__ARM_ARCH_5T__);
SHOW(__ARM_ARCH_5TE__);
SHOW(__ARM_ARCH_5TEJ__);
SHOW(__ARM_ARCH_6__);
SHOW(__ARM_ARCH_6J__);
SHOW(__ARM_ARCH_6K__);
SHOW(__ARM_ARCH_6Z__);
SHOW(__ARM_ARCH_6ZK__);
SHOW(__ARM_ARCH_6T2__);
SHOW(__ARM_ARCH_7__);
SHOW(__ARM_ARCH_7A__);
SHOW(__ARM_ARCH_7R__);
SHOW(__ARM_ARCH_7M__);
SHOW(__ARM_ARCH_7S__);

SHOW(__aarch64__);

SHOW(__bfin);
SHOW(__BFIN__);

SHOW(__convex__);
SHOW(__convex_c1__);
SHOW(__convex_c2__);
SHOW(__convex_c32__);
SHOW(__convex_c34__);
SHOW(__convex_c38__);

SHOW(__epiphany__);

SHOW(__hppa);
SHOW(__hppa__);
SHOW(__HPPA__);
SHOW(_PA_RISC1_0);
SHOW(_PA_RISC1_1);
SHOW(__HPPA11__);
SHOW(__PA7100__);
SHOW(_PA_RISC2_0);
SHOW(__RISC2_0__);
SHOW(__HPPA20__);
SHOW(__PA8000__);

SHOW(__ia64);
SHOW(__ia64__);
SHOW(_IA64);
SHOW(__IA64__);
SHOW(_M_IA64);
SHOW(__itanium__);

SHOW(__m68k__);
SHOW(__MC68K__);
SHOW(M68000);
SHOW(__mc68000__);
SHOW(__MC68000__);
SHOW(__mc68010__);
SHOW(__mc68020__);
SHOW(__MC68020__);
SHOW(__mc68030__);
SHOW(__MC68030__);
SHOW(__mc68040__);
SHOW(__mc68060__);

SHOW(__mips__);
SHOW(__mips);
SHOW(mips);
SHOW(__MIPS__);
SHOW(_R3000);
SHOW(_R4000);
SHOW(_R5900);
SHOW(_MIPS_ISA_MIPS1);
SHOW(_MIPS_ISA_MIPS2);
SHOW(_MIPS_ISA_MIPS3);
SHOW(_MIPS_ISA_MIPS4);
SHOW(__MIPS_ISA2__);
SHOW(__MIPS_ISA3__);
SHOW(__MIPS_ISA4__);

SHOW(__powerpc__);
SHOW(__powerpc);
SHOW(__powerpc64__);
SHOW(__POWERPC__);
SHOW(__ppc__);
SHOW(__ppc);
SHOW(__ppc64__);
SHOW(__PPC__);
SHOW(__PPC64__);
SHOW(_ARCH_PPC);
SHOW(_ARCH_PPC64);
SHOW(_M_PPC);
SHOW(__PPCGECKO__);
SHOW(__PPCBROADWAY__);
SHOW(_XENON);
SHOW(__ppc601__);
SHOW(__ppc603__);
SHOW(__ppc604__);
SHOW(_ARCH_440);
SHOW(_ARCH_450);
SHOW(_ARCH_601);
SHOW(_ARCH_603);
SHOW(_ARCH_604);
SHOW(TARGET_CPU_PPC);

SHOW(pyr);

SHOW(__THW_RS6000);
SHOW(_IBMR2);
SHOW(_POWER);
SHOW(_ARCH_PWR);
SHOW(_ARCH_PWR2);
SHOW(_ARCH_PWR3);
SHOW(_ARCH_PWR4);

SHOW(__sparc__);
SHOW(__sparc);
SHOW(__sparc_v8__);
SHOW(__sparc_v9__);
SHOW(__sparcv8);
SHOW(__sparcv9);

SHOW(__sh__);
SHOW(__sh1__);
SHOW(__sh2__);
SHOW(__sh3__);
SHOW(__SH3__);
SHOW(__SH4__);
SHOW(__SH5__);

SHOW(__370__);
SHOW(__THW_370__);
SHOW(__s390__);
SHOW(__s390x__);
SHOW(__zarch__);
SHOW(__SYSC_ZARCH__);

SHOW(_TMS320C2XX);
SHOW(__TMS320C2000__);
SHOW(_TMS320C5X);
SHOW(__TMS320C55X__);
SHOW(_TMS320C6X);
SHOW(__TMS320C6X__);
SHOW(_TMS320C28X);
SHOW(_TMS320C5XX);
SHOW(__TMS320C55X__);
SHOW(_TMS320C6200);
SHOW(_TMS320C6400);
SHOW(_TMS320C6400_PLUS);
SHOW(_TMS320C6600);
SHOW(_TMS320C6700);
SHOW(_TMS320C6700_PLUS);
SHOW(_TMS320C6740);

SHOW(__TMS470__);

SHOW(__VOS__);

SHOW(__WORDSIZE);
SHOW(__WORDSIZE32_PTRDIFF_LONG);
SHOW(__WORDSIZE32_SIZE_ULONG);
SHOW(sun);
SHOW(__sun);
SHOW(sparc);
SHOW(__sparc);
SHOW(__ARMEB__);
SHOW(__ARMEL__);
SHOW(__THUMBEB__);
SHOW(__THUMBEL__);
SHOW(__AARCH64EB__);
SHOW(__AARCH64EL__);
SHOW(_MIPSEB);
SHOW(_MIPSEL);
SHOW(__MMX__);
SHOW(__SSE__);
SHOW(__SSE2__);
SHOW(__SSE3__);
SHOW(__SSE_MATH__);
SHOW(__SSE2_MATH__);
SHOW(__SSSE3__);
SHOW(__SSE4_1__);
SHOW(__SSE4_2__);
SHOW(__AVX__);
SHOW(__AVX2__);
SHOW(__BIG_ENDIAN__);
SHOW(__LITTLE_ENDIAN__);
SHOW(__PDP_ENDIAN);
SHOW(__BYTE_ORDER__);
SHOW(__ORDER_LITTLE_ENDIAN__);
SHOW(__ORDER_BIG_ENDIAN__);
SHOW(__ORDER_PDP_ENDIAN__);
SHOW(__ARMEB__);
SHOW(__ARMEL__);
SHOW(__THUMBEB__);
SHOW(__THUMBEL__);
SHOW(__AARCH64EB__);
SHOW(__AARCH64EL__);
SHOW(_MIPSEB);
SHOW(_MIPSEL);
SHOW(__MIPSEB);
SHOW(__MIPSEL);
SHOW(__MIPSEB__);
SHOW(__MIPSEL__);
SHOW(__FLOAT_WORD_ORDER__);

#ifdef __GNUC__
/* Built-in functions are compiler-dependent */
printf("\nCPU type and supported features:\n");
__builtin_cpu_init ();
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 8 )
/* No exactly sure at which compiler version these were known. */
/* Right now I don't have a gcc version before 4.8.5 to test this. */
/* Since clang/clang++ identifies itself as a 4.2.x that is left out. */
CPU_IS("amd");
CPU_IS("intel");
CPU_IS("atom");
CPU_IS("core2");
CPU_IS("corei7")
CPU_IS("nehalem");
CPU_IS("westmere");
CPU_IS("sandybridge");
CPU_IS("amdfam10h");
CPU_IS("barcelona");
CPU_IS("shanghai");
CPU_IS("istanbul");
CPU_IS("amdfam15h");
CPU_IS("bdver1");
CPU_IS("bdver2");
CPU_IS("bdver3");
#endif
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 8 )
CPU_IS("slm");
CPU_IS("ivybridge");
CPU_IS("haswell");
CPU_IS("broadwell");
CPU_IS("bonnell");
CPU_IS("silvermont");
CPU_IS("bdver4");
CPU_IS("btver1");
CPU_IS("btver2");
#endif
#if __GNUC__ > 4
CPU_IS("knl");
#endif
#if __GNUC__ > 5
CPU_IS("skylake");
CPU_IS("skylake-avx512");
CPU_IS("amdfam17h");
CPU_IS("znver1");
#endif
#if __GNUC__ > 7
CPU_IS("cannonlake");
CPU_IS("icelake-client");
CPU_IS("icelake-server");
CPU_IS("knm");
#endif
#if __GNUC__ > 8
CPU_IS("cascadelake");
CPU_IS("goldmont");
CPU_IS("goldmont-plus");
CPU_IS("tremont");
CPU_IS("znver2");
#endif
#if __GNUC__ > 9
CPU_IS("tigerlake");
CPU_IS("cooperlake");
#endif
#if __GNUC__ > 10
CPU_IS("sapphirerapids");
CPU_IS("alderlake");
CPU_IS("rocketlake");
CPU_IS("amdfam19h");
CPU_IS("znver3");
#endif
#if __GNUC__ > 12
CPU_IS("graniterapids");
CPU_IS("graniterapids-d");
CPU_IS("sierraforest");
CPU_IS("grandridge");
CPU_IS("lujiazui");
CPU_IS("znver4");
#endif
#if __GNUC__ > 13
CPU_IS("arrowlake");
CPU_IS("arrowlake-s");
CPU_IS("pantherlake");
CPU_IS("clearwaterforest");
CPU_IS("yongfeng");
#endif

#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 8 )
/* No exactly sure at which compiler version these were known */
CPU_SUPPORTS("cmov");
CPU_SUPPORTS("mmx");
CPU_SUPPORTS("popcnt");
CPU_SUPPORTS("sse");
CPU_SUPPORTS("sse2");
CPU_SUPPORTS("sse3");
CPU_SUPPORTS("ssse3");
CPU_SUPPORTS("sse4.1");
CPU_SUPPORTS("sse4.2");
CPU_SUPPORTS("avx");
CPU_SUPPORTS("avx2");
#endif
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 8 )
CPU_SUPPORTS("sse4a");
CPU_SUPPORTS("fma4");
CPU_SUPPORTS("xop");
CPU_SUPPORTS("fma");
#endif
#if __GNUC__ > 4
CPU_SUPPORTS("avx512f");
CPU_SUPPORTS("bmi");
CPU_SUPPORTS("bmi2");
#endif
#if __GNUC__ > 5
CPU_SUPPORTS("aes");
CPU_SUPPORTS("pclmul");
CPU_SUPPORTS("avx512vl");
CPU_SUPPORTS("avx512bw");
CPU_SUPPORTS("avx512dq");
CPU_SUPPORTS("avx512cd");
CPU_SUPPORTS("avx512er");
CPU_SUPPORTS("avx512pf");
CPU_SUPPORTS("avx512vbmi");
CPU_SUPPORTS("avx512ifma");
#endif
#if __GNUC__ > 6
CPU_SUPPORTS("avx5124vnniw");
CPU_SUPPORTS("avx5124fmaps");
CPU_SUPPORTS("avx512vpopcntdq");
#endif
#if __GNUC__ > 7
CPU_SUPPORTS("avx512vbmi2");
CPU_SUPPORTS("gfni");
CPU_SUPPORTS("vpclmulqdq");
CPU_SUPPORTS("avx512vnni");
CPU_SUPPORTS("avx512bitalg");
#endif
#if __GNUC__ > 11
CPU_SUPPORTS("x86-64");
CPU_SUPPORTS("x86-64-v2");
CPU_SUPPORTS("x86-64-v3");
CPU_SUPPORTS("x86-64-v4");
#endif
#endif

printf("\nSystem:\n");
SHOW(MSDOS);
SHOW(__MSDOS__);
SHOW(vax);
SHOW(vms);
SHOW(linux);
SHOW(__linux__);
SHOW(__gnu_linux__);
SHOW(ultrix);
SHOW(lynx);
SHOW(__lynx__);
SHOW(AIX);
SHOW(__osf__);
SHOW(__alpha);
SHOW(__alpha__);
SHOW(__APPLE__);
SHOW(__MACH__);
SHOW(__ENVIRONMENT_MAC_OS_X_VERSION_MIN_REQUIRED__);
SHOW(OSK);
SHOW(M68K);
SHOW(__hpux);
SHOW(__hppa);
SHOW(unix);
SHOW(__unix);
SHOW(BSD);
SHOW(_Windows_);
SHOW(__NT__);
SHOW(_WIN32);
SHOW(_WIN64);
SHOW(WIN32_PLATFORM_HPCPRO);
SHOW(_WIN32_WCE);

printf("\nFilesystem limits:\n");
printf("pathconf(\"/\",_PC_PATH_MAX) = %ld\n", pathconf("/",_PC_PATH_MAX));
printf("pathconf(\".\",_PC_PATH_MAX) = %ld\n", pathconf(".",_PC_PATH_MAX));
SHOW(PATH_MAX);
SHOW(_POSIX_PATH_MAX);
SHOW(MAXPATH);
SHOW(MAXPATHLEN);
printf("pathconf(\"/\",_PC_NAME_MAX}) = %ld\n", pathconf("/",_PC_NAME_MAX));
printf("pathconf(\".\",_PC_NAME_MAX}) = %ld\n", pathconf(".",_PC_NAME_MAX));
SHOW(FILENAME_MAX);
SHOW(NAME_MAX);
SHOW(_POSIX_NAME_MAX);

printf("\nDefinition of variable types:\n");
SHOW(__SIZE_TYPE__);
SHOW(__PTRDIFF_TYPE__);
SHOW(__BOOL_TYPE__);
SHOW(__CHAR_TYPE__);
SHOW(__WCHAR_TYPE__);
SHOW(__WINT_TYPE__);
SHOW(__INTMAX_TYPE__);
SHOW(__UINTMAX_TYPE__);
SHOW(__SIG_ATOMIC_TYPE__);
SHOW(__INT8_TYPE__);
SHOW(__INT16_TYPE__);
SHOW(__INT32_TYPE__);
SHOW(__INT64_TYPE__);
SHOW(__UINT8_TYPE__);
SHOW(__UINT16_TYPE__);
SHOW(__UINT32_TYPE__);
SHOW(__UINT64_TYPE__);
SHOW(__INT_LEAST8_TYPE__);
SHOW(__INT_LEAST16_TYPE__);
SHOW(__INT_LEAST32_TYPE__);
SHOW(__INT_LEAST64_TYPE__);
SHOW(__UINT_LEAST8_TYPE__);
SHOW(__UINT_LEAST16_TYPE__);
SHOW(__UINT_LEAST32_TYPE__);
SHOW(__UINT_LEAST64_TYPE__);
SHOW(__INT_FAST8_TYPE__);
SHOW(__INT_FAST16_TYPE__);
SHOW(__INT_FAST32_TYPE__);
SHOW(__INT_FAST64_TYPE__);
SHOW(__UINT_FAST8_TYPE__);
SHOW(__UINT_FAST16_TYPE__);
SHOW(__UINT_FAST32_TYPE__);
SHOW(__UINT_FAST64_TYPE__);
SHOW(__INTPTR_TYPE__);
SHOW(__UINTPTR_TYPE__);

SHOW(__DEV_T_TYPE);
SHOW(__UID_T_TYPE);
SHOW(__GID_T_TYPE);
SHOW(__INO_T_TYPE);
SHOW(__INO64_T_TYPE);
SHOW(__MODE_T_TYPE);
SHOW(__NLINK_T_TYPE);
SHOW(__OFF_T_TYPE);
SHOW(__OFF64_T_TYPE);
SHOW(__PID_T_TYPE);
SHOW(__FSID_T_TYPE);
SHOW(__CLOCK_T_TYPE);
SHOW(__ID_T_TYPE);
SHOW(__TIME_T_TYPE);
SHOW(__USECONDS_T_TYPE);
SHOW(__SUSECONDS_T_TYPE);
SHOW(__SUSECONDS64_T_TYPE);

printf("\nAvailable floating point sizes:\n");
SHOW(__HAVE_FLOAT16);
SHOW(__HAVE_FLOAT32);
SHOW(__HAVE_FLOAT64);
SHOW(__HAVE_FLOAT80);
SHOW(__HAVE_FLOAT128);
SHOW(__HAVE_FLOAT32X);
SHOW(__HAVE_FLOAT64X);
SHOW(__HAVE_FLOAT80X);
SHOW(__HAVE_FLOAT128X);
SHOW(__HAVE_DISTINCT_FLOAT16);
SHOW(__HAVE_DISTINCT_FLOAT128);
SHOW(__LDBL_COMPAT);

printf("\nSize of variable types:\n");
printf("sizeof(bool) = %lu\n", sizeof(bool));
printf("sizeof(char) = %lu\n", sizeof(char));
printf("sizeof(wchar_t) = %lu\n", sizeof(wchar_t));
printf("sizeof(short) = %lu\n", sizeof(short));
printf("sizeof(int) = %lu\n", sizeof(int));
printf("sizeof(long) = %lu\n", sizeof(long));
printf("sizeof(long long) = %lu\n", sizeof(long long));

printf("sizeof(int_least8_t) = %lu\n", sizeof(int_least8_t));
printf("sizeof(int_least16_t) = %lu\n", sizeof(int_fast16_t));
printf("sizeof(int_least32_t) = %lu\n", sizeof(int_fast32_t));
printf("sizeof(int_least64_t) = %lu\n", sizeof(int_fast64_t));
printf("sizeof(int_fast8_t) = %lu\n", sizeof(int_fast8_t));
printf("sizeof(int_fast16_t) = %lu\n", sizeof(int_fast16_t));
printf("sizeof(int_fast32_t) = %lu\n", sizeof(int_fast32_t));
printf("sizeof(int_fast64_t) = %lu\n", sizeof(int_fast64_t));
printf("sizeof(uint_least8_t) = %lu\n", sizeof(uint_least8_t));
printf("sizeof(uint_least16_t) = %lu\n", sizeof(uint_fast16_t));
printf("sizeof(uint_least32_t) = %lu\n", sizeof(uint_fast32_t));
printf("sizeof(uint_least64_t) = %lu\n", sizeof(uint_fast64_t));
printf("sizeof(uint_fast8_t) = %lu\n", sizeof(uint_fast8_t));
printf("sizeof(uint_fast16_t) = %lu\n", sizeof(uint_fast16_t));
printf("sizeof(uint_fast32_t) = %lu\n", sizeof(uint_fast32_t));
printf("sizeof(uint_fast64_t) = %lu\n", sizeof(uint_fast64_t));

printf("sizeof(float) = %lu\n", sizeof(float));
printf("sizeof(double) = %lu\n", sizeof(double));
printf("sizeof(long double) = %lu\n", sizeof(long double));

SHOW(__CHAR_BIT__);
SHOW(__SIZEOF_INT__);
SHOW(__SIZEOF_LONG__);
SHOW(__SIZEOF_LONG_LONG__);
SHOW(__SIZEOF_SHORT__);
SHOW(__SIZEOF_POINTER__);
SHOW(__SIZEOF_FLOAT__);
SHOW(__SIZEOF_DOUBLE__);
SHOW(__SIZEOF_LONG_DOUBLE__);
SHOW(__SIZEOF_SIZE_T__);
SHOW(__SIZEOF_WCHAR_T__);
SHOW(__SIZEOF_WINT_T__);
SHOW(__SIZEOF_POINTER__);
SHOW(__SIZEOF_PTRDIFF_T__);
SHOW(__SIZEOF_INT128__);

SHOW(__SIZEOF_SEM_T);
SHOW(__SIZEOF_PTHREAD_MUTEX_T);
SHOW(__SIZEOF_PTHREAD_MUTEXATTR_T);
SHOW(__SIZEOF_PTHREAD_ATTR_T);
SHOW(__SIZEOF_PTHREAD_COND_T);
SHOW(__SIZEOF_PTHREAD_CONDATTR_T);
SHOW(__SIZEOF_PTHREAD_RWLOCK_T);
SHOW(__SIZEOF_PTHREAD_RWLOCKATTR_T);
SHOW(__SIZEOF_PTHREAD_BARRIER_T);
SHOW(__SIZEOF_PTHREAD_BARRIERATTR_T);

printf("\nLimits for variable types:\n");
SHOW(__SCHAR_MAX__);
SHOW(__WCHAR_MAX__);
SHOW(__SHRT_MAX__);
SHOW(__INT_MAX__);
SHOW(__LONG_MAX__);
SHOW(__LONG_LONG_MAX__);
SHOW(__WINT_MAX__);
SHOW(__SIZE_MAX__);
SHOW(__PTRDIFF_MAX__);
SHOW(__INTMAX_MAX__);
SHOW(__UINTMAX_MAX__);
SHOW(__SIG_ATOMIC_MAX__);
SHOW(__INT8_MAX__);
SHOW(__INT16_MAX__);
SHOW(__INT32_MAX__);
SHOW(__INT64_MAX__);
SHOW(__UINT8_MAX__);
SHOW(__UINT16_MAX__);
SHOW(__UINT32_MAX__);
SHOW(__UINT64_MAX__);
SHOW(__INT_LEAST8_MAX__);
SHOW(__INT_LEAST16_MAX__);
SHOW(__INT_LEAST32_MAX__);
SHOW(__INT_LEAST64_MAX__);
SHOW(__UINT_LEAST8_MAX__);
SHOW(__UINT_LEAST16_MAX__);
SHOW(__UINT_LEAST32_MAX__);
SHOW(__UINT_LEAST64_MAX__);
SHOW(__INT_FAST8_MAX__);
SHOW(__INT_FAST16_MAX__);
SHOW(__INT_FAST32_MAX__);
SHOW(__INT_FAST64_MAX__);
SHOW(__UINT_FAST8_MAX__);
SHOW(__UINT_FAST16_MAX__);
SHOW(__UINT_FAST32_MAX__);
SHOW(__UINT_FAST64_MAX__);
SHOW(__INTPTR_MAX__);
SHOW(__UINTPTR_MAX__);
SHOW(__WCHAR_MIN__);
SHOW(__WINT_MIN__);
SHOW(__SIG_ATOMIC_MIN__);

printf("\nOther:\n");
SHOW(lint);
SHOW(__INTMAX_TYPE__);
SHOW(__UINTMAX_TYPE__);
SHOW(INTMAX_MAX);
SHOW(UINTMAX_MAX);
SHOW(__CHAR_UNSIGNED__);
SHOW(__WCHAR_UNSIGNED__);
SHOW(__EXCEPTIONS);
SHOW(M_PI);
SHOW(M_PI_2);
SHOW(M_PI_4);
SHOW(M_PIl);
SHOW(M_PI_2l);
SHOW(M_PI_4l);
SHOW(M_E);
SHOW(M_LN2);
SHOW(M_LN10);

printf("\nCompiler features:\n");
#ifdef __has_include
   printf("   __has_include\n");
#endif
#ifdef __has_builtin
   printf("   __has_builtin\n");
#endif
#ifdef __has_cpp_attribute
   printf("   __has_cpp_attribute\n");
#else
# define __has_cpp_attribute(a) 0
#endif
#ifdef __cpp_attributes
   printf("   __cpp_attributes\n");
#endif
#ifdef __has_c_attribute
   printf("   __has_c_attribute\n");
#endif
#ifdef __has_extension
   printf("   __has_extension\n");
#endif
#if defined __clang__ && defined __has_feature
   printf("   __has_feature\n");
# if !defined(has_feature)
/*
//   #define has_feature(ftname) __has_feature(cxx_ ## ftname) || __has_feature(cpp_ ## ftname) || 
//           __has_feature(__cpp_ ## ftname) || defined(__cpp_ ## ftname)
*/
#   define has_feature(ftname)  __has_feature(__cpp_ ## ftname)
#   define has_feature2(ftname) __has_feature(cxx_ ## ftname)
#   define has_feature3(ftname) __has_feature(cpp_ ## ftname)
# endif
# if !defined(has_lib_feature)
/*
//   #define has_lib_feature(ftname) __has_feature(__cpp_lib_ ## ftname) ||  __has_feature(cxx_ ## ftname) || 
//            defined(__cpp_lib_ ## ftname)
*/
#   define has_lib_feature(ftname) __has_feature(__cpp_lib_ ## ftname)
#   define has_lib_feature2(ftname) __has_feature(cxx_lib ## ftname) 
# endif
   printf("   Features (*) might get shown with the standard proposed '__cpp_...' macros\n"
          "   instead of actuall clang-specific variant of '__has_feature(cxx_...)' or similar.\n");
#else
/* No C++ and/or no __has_feature */
#  define has_feature(ftname) __cpp_ ## ftname
#  define has_feature2(ftname) 0
#  define has_feature3(ftname) 0
#  define has_lib_feature(ftname) __cpp_lib_ ## ftname
#  define has_lib_feature2(ftname) 0
#  define __has_feature(ftname) __cpp_ ## ftname
#endif

#ifndef has_feature
#  define has_feature(ftname) 0
#endif
#ifndef has_feature2
#  define has_feature2(ftname) 0
#endif
#ifndef has_feature3
#  define has_feature3(ftname) 0
#endif
#ifndef has_lib_feature
#  define has_lib_feature(ftname) 0
#endif
#ifndef has_lib_feature2
#  define has_lib_feature2(ftname) 0
#endif

#define has_feature_x(ftname) __has_feature(cxx_ ## ftname) || __has_feature(cpp_ ## ftname) || \
           __has_feature(__cpp_ ## ftname) || has_feature(ftname)

/* For feature tests see: https://en.cppreference.com/w/cpp/feature_test */

   /* ----------- C++98 features ------------ */

   printf("  C++98 features available:\n");

#  if defined (__cpp_rtti)
      printf("   __cpp_rtti = %ld\n",
          (long) __cpp_rtti);
#  elif __has_feature(__cpp_rtti)
      printf("   __has_feature(__cpp_rtti) = %ld\n",
          (long) __has_feature(__cpp_rtti));
#  elif __has_feature(cpp_rtti)
      printf("   __has_feature(cpp_rtti) = %ld\n",
          (long) __has_feature(cpp_rtti));
#  elif __has_feature(cxx_rtti)
      printf("   __has_feature(cxx_rtti) = %ld\n",
          (long) __has_feature(cxx_rtti));
#  elif has_feature(rtti)
      printf("   has_feature(rtti) = %ld (*)\n",
          (long) has_feature(rtti));
#  endif

#  if defined (__cpp_exceptions)
      printf("   __cpp_exceptions = %ld\n",
          (long) __cpp_exceptions);
#  elif __has_feature(__cpp_exceptions)
      printf("   __has_feature(__cpp_exceptions) = %ld\n",
          (long) __has_feature(__cpp_exceptions));
#  elif __has_feature(cpp_exceptions)
      printf("   __has_feature(cpp_exceptions) = %ld\n",
          (long) __has_feature(cpp_exceptions));
#  elif __has_feature(cxx_exceptions)
      printf("   __has_feature(cxx_exceptions) = %ld\n",
          (long) __has_feature(cxx_exceptions));
#  elif has_feature(exceptions)
      printf("   has_feature(exceptions) = %ld (*)\n",
          (long) has_feature(exceptions));
#  endif

   /* ----------- C++11 features ------------ */

   printf("  C++11 features available:\n");

#  if defined (__cpp_unicode_characters)
      printf("   __cpp_unicode_characters = %ld\n",
          (long) __cpp_unicode_characters);
#  elif has_feature(unicode_characters)
      printf("   has_feature(unicode_characters) = %ld (*)\n",
          (long) has_feature(unicode_characters));
#  endif

#  if defined (__cpp_raw_strings)
      printf("   __cpp_raw_strings = %ld\n",
          (long) __cpp_raw_strings);
#  elif __has_feature(raw_strings)
      printf("   __has_feature(raw_strings) = %ld\n",
          (long) __has_feature(raw_strings));
#  elif has_feature(raw_strings)
      printf("   has_feature(raw_strings) = %ld\n",
          (long) has_feature(raw_strings));
#  endif

#  if defined (__cpp_unicode_literals)
      printf("   __cpp_unicode_literals = %ld\n",
          (long) __cpp_unicode_literals);
#  elif __has_feature(unicode_literals)
      printf("   __has_feature(unicode_literals) = %ld\n",
          (long) __has_feature(unicode_literals));
#  elif has_feature(unicode_literals)
      printf("   has_feature(unicode_literals) = %ld\n",
          (long) has_feature(unicode_literals));
#  endif

#  if defined(__cpp_user_defined_literals)
      printf("   __cpp_user_defined_literals = %ld\n",
          (long) __cpp_user_defined_literals);
#  elif __has_feature(user_defined_literals)
      printf("   __has_feature(user_defined_literals) = %ld\n",
         (long)  __has_feature(user_defined_literals));
#  elif has_feature(user_defined_literals)
      printf("   has_feature(user_defined_literals) = %ld\n",
          (long) has_feature(user_defined_literals));
#  endif

#  if defined(__cpp_threadsafe_static_init)
      printf("   __cpp_threadsafe_static_init = %ld\n",
          (long) __cpp_threadsafe_static_init);
#  endif
#  if defined(__cpp_lambdas)
      printf("   __cpp_lambdas = %ld\n",
          (long) __cpp_lambdas);
#  endif
#  if defined(__cpp_constexpr)
#    if ( __cpp_constexpr >= 202110L )
      printf("   __cpp_constexpr = %ld (C++23)\n",
          (long) __cpp_constexpr);
#    elif ( __cpp_constexpr >= 201806L )
      printf("   __cpp_constexpr = %ld (C++20)\n",
          (long) __cpp_constexpr);
#    elif ( __cpp_constexpr >= 201603L )
      printf("   __cpp_constexpr = %ld (C++17)\n",
          (long) __cpp_constexpr);
#    elif ( __cpp_constexpr >= 201304L )
      printf("   __cpp_constexpr = %ld (C++14)\n",
          (long) __cpp_constexpr);
#    else
      printf("   __cpp_constexpr = %ld (C++11)\n",
          (long) __cpp_constexpr);
#    endif
#  endif
#  if defined(__cpp_range_based_for)
#     if __cpp_range_based_for >= 202211L 
         printf("   __cpp_range_based_for = %ld (C++23)\n",
            (long) __cpp_range_based_for);
#     elif __cpp_range_based_for >= 201603L
         printf("   __cpp_range_based_for = %ld (C++17)\n",
            (long) __cpp_range_based_for);
#     else
         printf("   __cpp_range_based_for = %ld (C++11)\n",
            (long) __cpp_range_based_for);
#     endif
#  endif
#  if defined(__cpp_static_assert)
      printf("   __cpp_static_assert = %ld\n",
          (long) __cpp_static_assert);
#  endif
#  if defined(__cpp_decltype)
      printf("   __cpp_decltype = %ld\n",
          (long) __cpp_decltype);
#  endif
#  if defined(__cpp_attributes)
      printf("   __cpp_attributes = %ld\n",
          (long) __cpp_attributes);
#     if __has_cpp_attribute(noreturn)
      printf("   __has_cpp_attribute(noreturn) = %ld\n",
          (long) __has_cpp_attribute(noreturn));
#     endif
#     if __has_cpp_attribute(carries_dependency)
      printf("   __has_cpp_attribute(carries_dependency) = %ld\n",
          (long) __has_cpp_attribute(carries_dependency));
#     endif
#  endif
#  if defined(__cpp_rvalue_references)
      printf("   __cpp_rvalue_references = %ld\n",
          (long) __cpp_rvalue_references);
#  endif
#  if defined(__cpp_variadic_templates)
      printf("   __cpp_variadic_templates = %ld\n",
          (long) __cpp_variadic_templates);
#  endif
#  if defined(__cpp_initializer_lists)
      printf("   __cpp_initializer_lists = %ld\n",
          (long) __cpp_initializer_lists);
#  endif
#  if defined(__cpp_delegating_constructors)
      printf("   __cpp_delegating_constructors = %ld\n",
          (long) __cpp_delegating_constructors);
#  endif
#  if defined(__cpp_nsdmi)
      printf("   __cpp_nsdmi = %ld\n", (long) __cpp_nsdmi);
#  endif
#  if defined(__cpp_inheriting_constructors)
      printf("   __cpp_inheriting_constructors = %ld\n",
          (long) __cpp_inheriting_constructors);
#  endif
#  if defined(__cpp_ref_qualifiers)
      printf("   __cpp_ref_qualifiers = %ld\n",
          (long) __cpp_ref_qualifiers);
#  endif
#  if defined(__cpp_alias_templates)
      printf("   __cpp_alias_templates = %ld\n",
          (long) __cpp_alias_templates);
#  endif

   /* ----------- C++14 features ------------ */

   printf("  C++14 features available:\n");
#  if defined(__cpp_binary_literals)
      printf("   __cpp_binary_literals = %ld\n",
          (long) __cpp_binary_literals);
#  endif
#  if defined(__cpp_init_captures)
      printf("   __cpp_init_captures = %ld\n",
          (long) __cpp_init_captures);
#  endif
#  if defined(__cpp_generic_lambdas)
      printf("   __cpp_generic_lambdas = %ld\n",
          (long) __cpp_generic_lambdas);
#  endif
#  if defined(__cpp_sized_deallocation)
      printf("   __cpp_sized_deallocation = %ld\n",
          (long) __cpp_sized_deallocation);
#  endif
#  if defined(__cpp_constexpr)
     if ( has_feature(constexpr) > 201200 )
      printf("   __cpp_constexpr = %ld\n",
          (long) __cpp_constexpr);
#  endif
#  if defined(__cpp_decltype_auto)
      printf("   __cpp_decltype_auto = %ld\n",
          (long) __cpp_decltype_auto);
#  endif
#  if defined(__cpp_return_type_deduction)
      printf("   __cpp_return_type_deduction = %ld\n",
          (long) __cpp_return_type_deduction);
#  endif
#if __has_cpp_attribute(deprecated)
      printf("   __has_cpp_attribute(deprecated) = %ld\n",
          (long) __has_cpp_attribute(deprecated));
#endif
#  if defined(__cpp_aggregate_nsdmi)
      printf("   __cpp_aggregate_nsdmi = %ld\n",
          (long) __cpp_aggregate_nsdmi);
#  endif
#  if defined(__cpp_variable_templates)
      printf("   __cpp_variable_templates = %ld\n",
          (long) __cpp_variable_templates);
#  endif
#  if defined(__cpp_lib_integer_sequence)
      printf("   __cpp_lib_integer_sequence = %ld (depends on <utility>)\n",
          (long) __cpp_lib_integer_sequence); // <utility>
#  endif
#  if defined(__cpp_lib_exchange_function)
      printf("   __cpp_lib_exchange_function = %ld (depends on <utility>)\n",
          (long) __cpp_lib_exchange_function); // <utility>
#  endif
#  if defined(__cpp_lib_tuples_by_type)
      printf("   __cpp_lib_tuples_by_type = %ld (depends on <utility>)\n",
          (long) __cpp_lib_tuples_by_type); // <utility>
#  endif
#  if defined(__cpp_lib_tuple_element_t)
      printf("   __cpp_lib_tuple_element_t = %ld (depends on <utility>)\n",
          (long) __cpp_lib_tuple_element_t); // <utility>
#  endif
#  if defined(__cpp_lib_make_unique)
      printf("   __cpp_lib_make_unique = %ld (depends on <memory>)\n",
          (long) __cpp_lib_make_unique); // <memory>
#  endif
#  if defined(__cpp_lib_transparent_operators)
      printf("   __cpp_lib_transparent_operators = %ld (depends on <functional>)\n",
          (long) __cpp_lib_transparent_operators); // <functional>
#  endif
#  if defined(__cpp_lib_integral_constant_callable)
      printf("   __cpp_lib_integral_constant_callable = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_integral_constant_callable); // <type_traits>
#  endif
#  if defined(__cpp_lib_transformation_trait_aliases)
      printf("   __cpp_lib_transformation_trait_aliases = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_transformation_trait_aliases); // <type_traits>
#  endif
#  if defined(__cpp_lib_result_of_sfinae)
      printf("   __cpp_lib_result_of_sfinae = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_result_of_sfinae); // <functional> <type_traits>
#  endif
#  if defined(__cpp_lib_is_final)
      printf("   __cpp_lib_is_final = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_is_final); // <type_traits>
#  endif
#  if defined(__cpp_lib_is_null_pointer)
      printf("   __cpp_lib_is_null_pointer = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_is_null_pointer); // <type_traits>
#  endif
#  if defined(__cpp_lib_chrono_udls)
      printf("   __cpp_lib_chrono_udls = %ld (depends on <chrono>)\n",
          (long) __cpp_lib_chrono_udls); // <chrono>
#  endif
#  if defined(__cpp_lib_string_udls)
      printf("   __cpp_lib_string_udls = %ld (depends on <string>)\n",
          (long) __cpp_lib_string_udls); // <string>
#  endif
#  if defined(__cpp_lib_generic_associative_lookup)
      printf("   __cpp_lib_generic_associative_lookup = %ld (depends on  <map>, <set>)\n",
          (long) __cpp_lib_generic_associative_lookup); // <map> <set>
#  endif
#  if defined(__cpp_lib_null_iterators)
      printf("   __cpp_lib_null_iterators = %ld (depends on <iterator>)\n",
          (long) __cpp_lib_null_iterators); // <iterator>
#  endif
#  if defined(__cpp_lib_make_reverse_iterator)
      printf("   __cpp_lib_make_reverse_iterator = %ld (depends on <iterator>)\n",
          (long) __cpp_lib_make_reverse_iterator); // <iterator>
#  endif
#  if defined(__cpp_lib_robust_nonmodifying_seq_ops)
      printf("   __cpp_lib_robust_nonmodifying_seq_ops = %ld (depends on <algorithm>)\n",
          (long) __cpp_lib_robust_nonmodifying_seq_ops); // <algorithm>
#  endif
#  if defined(__cpp_lib_complex_udls)
      printf("   __cpp_lib_complex_udls = %ld (depends on <complex>)\n",
          (long) __cpp_lib_complex_udls); // <complex>
#  endif
#  if defined(__cpp_lib_)
      printf("   __cpp_lib_quoted_string_io = %ld (depends on <iomanip>)\n",
          (long) __cpp_lib_); // <iomanip>
#  endif
#  if defined(__cpp_lib_shared_timed_mutex)
      printf("   __cpp_lib_shared_timed_mutex = %ld (depends on <shared_mutex>)\n",
          (long) __cpp_lib_shared_timed_mutex); // <shared_mutex>
#  endif

   /* ----------- C++17 features ------------ */

   printf("  C++17 features available:\n");

#  if defined(__cpp_hex_float)
      printf("   __cpp_hex_float = %ld\n",
          (long) __cpp_hex_float);
#  elif __has_feature(hex_float)
      printf("   __has_feature(hex_float) = %ld\n",
          (long) __has_feature(hex_float));
#  elif has_feature(hex_float)
      printf("   has_feature(hex_float) = %ld\n",
          (long) has_feature(hex_float));
#  elif defined(cxx_hex_float)
      printf("   cxx_hex_float) = %ld\n",
          (long) cxx_hex_float);
#  elif defined(cpp_hex_float)
      printf("   cpp_hex_float = %ld\n",
          (long) cpp_hex_float);
#  endif

#  if defined(__cpp_inline_variables)
      printf("   __cpp_inline_variables = %ld\n",
          (long) __cpp_inline_variables);
#  elif has_feature(inline_variables)
      printf("   has_feature(inline_variables) = %ld\n",
          (long) has_feature(inline_variables));
#  elif has_feature_x(inline_variables)
      printf("   has feature: inline_variables\n",
          (long) has_feature_x(inline_variables));
#  endif

#  if defined(__cpp_aligned_new)
      printf("   __cpp_aligned_new = %ld\n",
          (long) __cpp_aligned_new);
#  endif
#  if defined(__cpp_guaranteed_copy_elision)
      printf("   __cpp_guaranteed_copy_elision = %ld\n",
          (long) __cpp_guaranteed_copy_elision);
#  endif
#  if defined(__cpp_noexcept_function_type)
      printf("   __cpp_noexcept_function_type = %ld\n",
          (long) __cpp_noexcept_function_type);
#  endif
#  if defined(__cpp_fold_expressions)
      printf("   __cpp_fold_expressions = %ld\n",
          (long) __cpp_fold_expressions);
#  endif
#  if defined(__cpp_capture_star_this)
      printf("   __cpp_capture_star_this = %ld\n",
          (long) __cpp_capture_star_this);
#  endif
#  if defined(__cpp_constexpr)
     if ( __cpp_constexpr >= 201411 )
      printf("   __cpp_constexpr = %ld\n",
          (long) __cpp_constexpr);
#  endif
#  if defined(__cpp_if_constexpr)
      printf("   __cpp_if_constexpr = %ld\n",
          (long) __cpp_if_constexpr);
#  endif
#  if defined(__cpp_range_based_for)
     if ( __cpp_range_based_for >= 201411 )
      printf("   __cpp_range_based_for = %ld\n",
          (long) __cpp_range_based_for);
#  endif
#  if defined(__cpp_static_assert)
     if ( __cpp_static_assert >= 201411 )
      printf("   __cpp_static_assert = %ld\n",
          (long) __cpp_static_assert);
#  endif
#  if defined(__cpp_deduction_guides)
      printf("   __cpp_deduction_guides = %ld\n",
          (long) __cpp_deduction_guides);
#  endif
#  if defined(__cpp_template_auto)
      printf("   __cpp_template_auto = %ld\n",  
          (long) __cpp_template_auto);
#  endif
#  if defined(__cpp_nontype_template_parameter_auto)
      printf("   __cpp_nontype_template_parameter_auto = %ld\n",
          (long) __cpp_nontype_template_parameter_auto);
#  endif
#  if defined(__cpp_structured_bindings)
      printf("   __cpp_structured_bindings = %ld\n",  
          (long) __cpp_structured_bindings);
#  endif
#  if defined(__cpp_namespace_attributes)
      printf("   __cpp_namespace_attributes = %ld\n",
          (long) __cpp_namespace_attributes);
#  endif
#if defined(__cpp_nested_namespace_definitions)
      printf("   __cpp_nested_namespace_definitions = %ld\n",  
          (long) __cpp_nested_namespace_definitions);
#endif
#  if defined(_cpp_enumerator_attributes)
      printf("   __cpp_enumerator_attributes = %ld\n",
          (long) __cpp_enumerator_attributes);
#  endif
#  if defined(__cpp_inheriting_constructors)
      printf("   __cpp_inheriting_constructors = %ld\n",
          (long) __cpp_inheriting_constructors);
#  endif
#  if defined(__cpp_variadic_using)
      printf("   __cpp_variadic_using = %ld\n",
          (long) __cpp_variadic_using);
#  endif
#  if __has_cpp_attribute(fallthrough)
      printf("   __has_cpp_attribute(fallthrough) = %ld\n",
          (long) __has_cpp_attribute(fallthrough));
#  endif
#  if __has_cpp_attribute(nodiscard)
      printf("   __has_cpp_attribute(nodiscard) = %ld\n",
          (long) __has_cpp_attribute(nodiscard));
#  endif
#  if __has_cpp_attribute(maybe_unused)
      printf("   __has_cpp_attribute(maybe_unused) = %ld\n",
          (long) __has_cpp_attribute(maybe_unused));
#  endif
#  if defined(__cpp_structured_bindings)
      printf("   __cpp_structured_bindings = %ld\n",
          (long) __cpp_structured_bindings);
#  endif
#  if defined(__cpp_aggregate_bases)
      printf("   __cpp_aggregate_bases = %ld\n",
          (long) __cpp_aggregate_bases);
#  endif
#  if defined(__cpp_nontype_template_args)
      printf("   __cpp_nontype_template_args = %ld\n",
          (long) __cpp_nontype_template_args);
#  endif
#  if defined(__cpp_template_template_args)
      printf("   __cpp_template_template_args = %ld\n",
          (long) __cpp_template_template_args);
#  endif
#  if defined(__cpp_fold_expressions)
      printf("   __cpp_fold_expressions = %ld\n",
          (long) __cpp_fold_expressions);
#  endif
#  if defined(__cpp_lib_byte)
      printf("   __cpp_lib_byte = %ld (depends on <cstddef>)\n",
          (long) __cpp_lib_byte);
#  endif
#  if defined(__cpp_lib_hardware_interference_size)
      printf("   __cpp_lib_hardware_interference_size = %ld (depends on <new>)\n",
          (long) __cpp_lib_hardware_interference_size);
#  endif
#  if defined(__cpp_lib_launder)
      printf("   __cpp_lib_launder = %ld (depends on <new>)\n",
          (long) __cpp_lib_launder);
#  endif
#  if defined(__cpp_lib_uncaught_exceptions)
      printf("   __cpp_lib_uncaught_exceptions = %ld (depends on <exception>)\n",
          (long) __cpp_lib_uncaught_exceptions);
#  endif
#  if defined(__cpp_lib_as_const)
      printf("   __cpp_lib_as_const = %ld (depends on <utility>)\n",
          (long) __cpp_lib_as_const);
#  endif
#  if defined(__cpp_lib_make_from_tuple)
      printf("   __cpp_lib_make_from_tuple = %ld (depends on <tuple>)\n",
          (long) __cpp_lib_make_from_tuple);
#  endif
#  if defined(__cpp_lib_apply)
      printf("   __cpp_lib_apply = %ld (depends on <tuple>)\n",
          (long) __cpp_lib_apply);
#  endif
#  if defined(__cpp_lib_optional)
#     if __cpp_lib_optional >= 202110L
         printf("   __cpp_lib_optional = %ld (C++23, depends on <optional>)\n",
            (long) __cpp_lib_optional);
#     elif __cpp_lib_optional >= 202106L
         printf("   __cpp_lib_optional = %ld (C++20 DR, depends on <optional>)\n",
            (long) __cpp_lib_optional);
#     else
         printf("   __cpp_lib_optional = %ld (C++17, depends on <optional>)\n",
            (long) __cpp_lib_optional);
#     endif
#  endif
#  if defined(__cpp_lib_any)
      printf("   __cpp_lib_any = %ld (depends on <any>)\n",
          (long) __cpp_lib_any);
#  endif
#  if defined(__cpp_lib_memory_resource)
      printf("   __cpp_lib_memory_resource = %ld (depends on <memory_resource>)\n",
          (long) __cpp_lib_memory_resource);
#  endif
#  if defined(__cpp_lib_boyer_moore_searcher)
      printf("   __cpp_lib_boyer_moore_searcher = %ld (depends on <functional>)\n",
          (long) __cpp_lib_boyer_moore_searcher);
#  endif
#  if defined(__cpp_lib_string_view)
      printf("   __cpp_lib_string_view = %ld (depends on <string_view>)\n",
          (long) __cpp_lib_string_view);
#  endif
#  if defined(__cpp_lib_sample)
      printf("   __cpp_lib_sample = %ld (depends on <algorithm>)\n",
          (long) __cpp_lib_sample);
#  endif
#  if defined(__cpp_lib_variant)
      printf("   __cpp_lib_variant = %ld (depends on <variant >)\n",
          (long) __cpp_lib_variant);
#  endif
#  if defined(__cpp_lib_addressof_constexpr)
      printf("   __cpp_lib_addressof_constexpr = %ld (depends on <memory>)\n",
          (long) __cpp_lib_addressof_constexpr);
#  endif
#  if defined(__cpp_lib_raw_memory_algorithms)
      printf("   __cpp_lib_raw_memory_algorithms = %ld (depends on <memory>)\n",
          (long) __cpp_lib_raw_memory_algorithms);
#  endif
#  if defined(__cpp_lib_transparent_operators)
      printf("   __cpp_lib_transparent_operators = %ld (depends on <memory>, <functional>)\n",
          (long) __cpp_lib_transparent_operators);
#  endif
#  if defined(__cpp_lib_enable_shared_from_this)
      printf("   __cpp_lib_enable_shared_from_this = %ld (depends on <memory>)\n",
          (long) __cpp_lib_enable_shared_from_this);
#  endif
#  if defined(shared_ptr_weak_type)
      printf("   __cpp_lib_shared_ptr_weak_type = %ld (depends on <memory>)\n",
          (long) __cpp_lib_shared_ptr_weak_type);
#  endif
#  if defined(__cpp_lib_shared_ptr_arrays)
      printf("   __cpp_lib_shared_ptr_arrays = %ld (depends on <memory>)\n",
          (long) __cpp_lib_shared_ptr_arrays);
#  endif
#  if defined(__cpp_lib_invoke)
      printf("   __cpp_lib_invoke = %ld (depends on <functional>)\n",
          (long) __cpp_lib_invoke);
#  endif
#  if defined(__cpp_lib_not_fn)
      printf("   __cpp_lib_not_fn = %ld (depends on <functional>)\n",
          (long) __cpp_lib_not_fn);
#  endif
#  if defined(__cpp_lib_void_t)
      printf("   __cpp_lib_void_t = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_void_t);
#  endif
#  if defined(__cpp_lib_bool_constant)
      printf("   __cpp_lib_bool_constant = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_bool_constant);
#  endif
#  if defined(__cpp_lib_type_trait_variable_templates)
      printf("   __cpp_lib_type_trait_variable_templates = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_type_trait_variable_templates);
#  endif
#  if defined(__cpp_lib_logical_traits)
      printf("   __cpp_lib_logical_traits = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_logical_traits);
#  endif
#  if defined(__cpp_lib_is_swappable)
      printf("   __cpp_lib_is_swappable = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_is_swappable);
#  endif
#  if defined(__cpp_lib_is_invocable)
      printf("   __cpp_lib_is_invocable = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_is_invocable);
#  endif
#  if defined(__cpp_lib_has_unique_object_representations)
      printf("   __cpp_lib_has_unique_object_representations = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_has_unique_object_representations);
#  endif
#  if defined(__cpp_lib_is_aggregate)
      printf("   __cpp_lib_is_aggregate = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_is_aggregate);
#  endif
#  if defined(__cpp_lib_chrono)
      printf("   __cpp_lib_chrono = %ld (depends on <chrono>)\n",
          (long) __cpp_lib_chrono);
#  endif
#  if defined(__cpp_lib_execution)
      printf("   __cpp_lib_execution = %ld (depends on <execution>)\n",
          (long) __cpp_lib_execution);
#  endif
#  if defined(__cpp_lib_parallel_algorithm)
      printf("   __cpp_lib_parallel_algorithm = %ld (depends on <algorithm>, <numeric>)\n",
          (long) __cpp_lib_parallel_algorithm);
#  endif
#  if defined(__cpp_lib_to_chars)
      printf("   __cpp_lib_to_chars = %ld (depends on <utility>)\n",
          (long) __cpp_lib_to_chars);
#  endif
#  if defined(__cpp_lib_string_view)
      printf("   __cpp_lib_string_view = %ld (depends on <string>, <string_view>)\n",
          (long) __cpp_lib_string_view);
#  endif
#  if defined(__cpp_lib_allocator_traits_is_always_equal)
      printf("   __cpp_lib_allocator_traits_is_always_equal = %ld (depends on <memory>, ...)\n",
          (long) __cpp_lib_allocator_traits_is_always_equal);
#  endif
#  if defined(__cpp_lib_incomplete_container_elements)
      printf("   __cpp_lib_incomplete_container_elements = %ld (depends on <forwardlist>, <list>, <vector>)\n",
          (long) __cpp_lib_incomplete_container_elements);
#  endif
#  if defined(__cpp_lib_map_try_emplace)
      printf("   __cpp_lib_map_try_emplace = %ld (depends on <map>)\n",
          (long) __cpp_lib_map_try_emplace);
#  endif
#  if defined(__cpp_lib_unordered_map_try_emplace)
      printf("   __cpp_lib_unordered_map_try_emplace = %ld (depends on <unordered_map>)\n",
          (long) __cpp_lib_unordered_map_try_emplace);
#  endif
#  if defined(__cpp_lib_node_extract)
      printf("   __cpp_lib_node_extract = %ld (depends on <map>, <set>, ...)\n",
          (long) __cpp_lib_node_extract);
#  endif
#  if defined(__cpp_lib_array_constexpr)
      printf("   __cpp_lib_array_constexpr = %ld (depends on <iterator>, <array>)\n",
          (long) __cpp_lib_array_constexpr);
#  endif
#  if defined(__cpp_lib_nonmember_container_access)
      printf("   __cpp_lib_nonmember_container_access = %ld (depends on <iterator>, ...)\n",
          (long) __cpp_lib_nonmember_container_access);
#  endif
#  if defined(__cpp_lib_clamp)
      printf("   __cpp_lib_clamp = %ld (depends on <algorithm>)\n",
          (long) __cpp_lib_clamp);
#  endif
#  if defined(__cpp_lib_gcd_lcm)
      printf("   __cpp_lib_gcd_lcm = %ld (depends on <numeric>)\n",
          (long) __cpp_lib_gcd_lcm);
#  endif
#  if defined(__cpp_lib_hypot)
      printf("   __cpp_lib_hypot = %ld (depends on <cmath>)\n",
          (long) __cpp_lib_hypot);
#  endif
#  if defined(__cpp_lib_math_special_functions)
      printf("   __cpp_lib_math_special_functions = %ld (depends on <cmath>)\n",
          (long) __cpp_lib_math_special_functions);
#  endif
#  if defined(__cpp_lib_filesystem)
      printf("   __cpp_lib_filesystem = %ld (depends on <filesystem>)\n",
          (long) __cpp_lib_filesystem);
#  endif
#  if defined(__cpp_lib_atomic_is_always_lock_free)
      printf("   __cpp_lib_atomic_is_always_lock_free = %ld (depends on <atomic>)\n",
          (long) __cpp_lib_atomic_is_always_lock_free);
#  endif
#  if defined(__cpp_lib_shared_mutex)
      printf("   __cpp_lib_shared_mutex = %ld (depends on <shared_mutex>)\n",
          (long) __cpp_lib_shared_mutex);
#  endif
#  if defined(__cpp_lib_lock_guard_variadic)
      printf("   __cpp_lib_lock_guard_variadic = %ld (depends on <thread>)\n",
          (long) __cpp_lib_lock_guard_variadic);
#  endif
#  if defined(__cpp_lib_scoped_lock)
      printf("   __cpp_lib_scoped_lock = %ld (depends on <mutex>)\n",
          (long) __cpp_lib_scoped_lock);
#  endif

   /* ----------- C++20 features ------------ */

   printf("  C++20 features available:\n");
#if defined(__cpp_modules)
      printf("   __cpp_modules = %ld\n",
          (long) __cpp_modules);
#endif
#if defined(__cpp_concepts)
      printf("   __cpp_concepts = %ld (%s)\n",
          (long) __cpp_concepts,
                 __cpp_concepts >= 201907 ? "C++20" : "pre-C++20");
#endif
/*
// #if defined(__cpp_deduction_guides)
// # if __cpp_deduction_guides > 201611
//       printf("   __cpp_deduction_guides = %ld\n",  (long) __cpp_deduction_guides);
// # endif
// #endif
*/
#if defined(__cpp_impl_three_way_comparison)
      printf("   __cpp_impl_three_way_comparison = %ld\n",
          (long) __cpp_impl_three_way_comparison);
#endif
#if defined(__cpp_constexpr)
# if __cpp_constexpr > 201603
      printf("   __cpp_constexpr = %ld\n",
         (long)  __cpp_constexpr);
# endif
#endif
#if defined(__cpp_constexpr_in_decltype)
# if __cpp_constexpr_in_decltype >= 201711
      printf("   __cpp_constexpr_in_decltype = %ld\n",
         (long)  __cpp_constexpr_in_decltype);
# endif
#endif
#if defined(__cpp_constexpr_dynamic_alloc)
      printf("   __cpp_constexpr_dynamic_alloc = %ld\n",
          (long) __cpp_constexpr_dynamic_alloc);
#endif
#if defined(__cpp_constinit)
      printf("   __cpp_constinit = %ld\n",
          (long) __cpp_constinit);
#endif
#if defined(__cpp_consteval)
      printf("   __cpp_consteval = %ld\n",  
          (long) __cpp_consteval);
#endif
#if defined(__cpp_coroutines)
      printf("   __cpp_coroutines = %ld\n",
         (long)  __cpp_coroutines);
#endif
#if defined(__cpp_aggregate_paren_init)
      printf("   __cpp_aggregate_paren_init = %ld\n",
          (long) __cpp_aggregate_paren_init);
#endif
#if defined(__cpp_char8_t)
      printf("   __cpp_char8_t = %ld\n",
          (long) __cpp_char8_t);
#endif
#if defined(__cpp_conditional_explicit)
      printf("   __cpp_conditional_explicit = %ld\n",
          (long) __cpp_conditional_explicit);
#endif
#if defined(__cpp_impl_destroying_delete)
      printf("   __cpp_impl_destroying_delete = %ld\n",
          (long) __cpp_impl_destroying_delete);
#endif
#if defined(__cpp_nontype_template_parameter_class)
      printf("   __cpp_nontype_template_parameter_class = %ld\n",
          (long) __cpp_nontype_template_parameter_class);
#endif
#if defined(__cpp_using_enum)
      printf("   __cpp_using_enum = %ld\n",
          (long) __cpp_using_enum);
#endif
#if __has_cpp_attribute(nodiscard)
# if __has_cpp_attribute(nodiscard) >= 201907
     printf("   __has_cpp_attribute(nodiscard) = %ld\n",
         (long) __has_cpp_attribute(nodiscard));
# endif
#endif
#if __has_cpp_attribute(likely)
      printf("   __has_cpp_attribute(likely) = %ld\n",
          (long) __has_cpp_attribute(likely));
#endif
#if __has_cpp_attribute(unlikely)
      printf("   __has_cpp_attribute(unlikely) = %ld\n",
          (long) __has_cpp_attribute(unlikely));
#endif
#if __has_cpp_attribute(no_unique_address)
      printf("   __has_cpp_attribute(no_unique_address) = %ld\n",
          (long) __has_cpp_attribute(no_unique_address));
#endif
#if defined(__cpp_designated_initializers)
      printf("   __cpp_designated_initializers = %ld\n",
          (long) __cpp_designated_initializers);
#endif
#if defined(__cpp_generic_lambdas)
#  if __cpp_generic_lambdas >= 201707
      printf("   __cpp_generic_lambdas = %ld (C++20)\n",
          (long) __cpp_generic_lambdas);
#  endif
#endif
#if defined(__cpp_init_captures)
#   if __cpp_init_captures >= 201803 
      printf("   __cpp_init_captures = %ld\n",
          (long) __cpp_init_captures);
#   endif
#endif
#if defined(__cpp_nontype_template_parameter_class)
      printf("   __cpp_nontype_template_parameter_class = %ld\n",  
          (long) __cpp_nontype_template_parameter_class);
#endif
#if defined(__cpp_impl_destroying_delete)
      printf("   __cpp_impl_destroying_delete = %ld\n",  
          (long) __cpp_impl_destroying_delete);
#endif
#if defined(__cpp_aggregate_paren_init)
      printf("   __cpp_aggregate_paren_init = %ld\n",  
          (long) __cpp_aggregate_paren_init);
#endif


/* __cpp_lib_....: */
#if defined(__cpp_lib_array_constexpr)
# if __cpp_lib_array_constexpr >= 201803
      printf("   __cpp_lib_array_constexpr = %ld\n",
          (long) __cpp_lib_array_constexpr);
# endif
#endif
#if defined(__cpp_lib_atomic_flag_test)
      printf("   __cpp_lib_atomic_flag_test = %ld (depends on <atomic>)\n",
          (long) __cpp_lib_atomic_flag_test);
#endif
#if defined(__cpp_lib_atomic_ref)
      printf("   __cpp_lib_atomic_ref = %ld (depends on <atomic>)\n",
          (long) __cpp_lib_atomic_ref);
#endif
#if defined(__cpp_lib_atomic_wait)
      printf("   __cpp_lib_atomic_wait = %ld (depends on <atomic>)\n",
          (long) __cpp_lib_atomic_wait);
#endif
#if defined(__cpp_lib_barrier)
      printf("   __cpp_lib_barrier = %ld (depends on <barrier>\n",
          (long) __cpp_lib_barrier);
#endif
#if defined(__cpp_lib_bind_front)
      printf("   __cpp_lib_bind_front = %ld (depends on <functional>)\n",
         (long)  __cpp_lib_bind_front);
#endif
#if defined(__cpp_lib_bit_cast)
      printf("   __cpp_lib_bit_cast = %ld (depends on <bit>)\n",
          (long) __cpp_lib_bit_cast);
#endif
#if defined(__cpp_lib_bitops)
      printf("   __cpp_lib_bitops = %ld (depends on <bit>)\n",
           (long) __cpp_lib_bitops);
#endif
#if defined(__cpp_lib_bounded_array_traits)
      printf("   __cpp_lib_bounded_array_traits = %ld (depends on <type_traits>)\n",
           (long) __cpp_lib_bounded_array_traits);
#endif
#if defined(__cpp_lib_concepts)
   #if __cpp_lib_concepts >= 202207L
      printf("   __cpp_lib_concepts = %ld (C++23, depends on <compare> <concepts>)\n",
          (long) __cpp_lib_concepts);
   #else
      printf("   __cpp_lib_concepts = %ld (C++20, depends on <concepts>)\n",
          (long) __cpp_lib_concepts);
# endif
#endif
#if defined(__cpp_lib_constexpr)
      printf("   __cpp_lib_constexpr = %ld\n",
          (long) __cpp_lib_constexpr);
#endif
#if defined(__cpp_lib_constexpr_algorithms)
      printf("   __cpp_lib_constexpr_algorithms = %ld (depends on <algorithm>)\n",
          (long) __cpp_lib_constexpr_algorithms);
#endif
#if defined(__cpp_lib_constexpr_dynamic_alloc)
      printf("   __cpp_lib_constexpr_dynamic_alloc = %ld (depends on <memory>)\n",
          (long) __cpp_lib_constexpr_dynamic_alloc);
#endif
#if defined(__cpp_lib_constexpr_invoke)
      printf("   __cpp_lib_constexpr_invoke = %ld (depends on <functional>)\n",
          (long) __cpp_lib_constexpr_invoke);
#endif
#if defined(__cpp_lib_constexpr_string)
      printf("   __cpp_lib_constexpr_string = %ld (depends on <string>)\n",
          (long) __cpp_lib_constexpr_string);
#endif
#if defined(__cpp_lib_constexpr_vector)
      printf("   __cpp_lib_constexpr_vector = %ld (depends on <vector>)\n",
          (long) __cpp_lib_constexpr_vector);
#endif
#if defined(__cpp_lib_destroying_delete)
      printf("   __cpp_lib_destroying_delete = %ld (depends on <new>)\n",
          (long) __cpp_lib_destroying_delete);
#endif
#if defined(__cpp_lib_endian)
      printf("   __cpp_lib_endian = %ld (depends on <bit>)\n",
          (long) __cpp_lib_endian);
#endif
#if defined(__cpp_lib_erase_if)
      printf("   __cpp_lib_erase_if = %ld (depends on <deque>, ...)\n",
          (long) __cpp_lib_erase_if);
#endif
#if defined(__cpp_lib_execution)
# if __cpp_lib_execution >= 201902
      printf("   __cpp_lib_execution = %ld (depends on <execution>)\n",
          (long) __cpp_lib_execution);
# endif
#endif
/*
// #if defined(__cpp_lib_filesystem)
// # if __cpp_lib_filesystem >= 201703
//       printf("   __cpp_lib_filesystem = %ld (depends on <filesystem>)\n", (long) __cpp_lib_filesystem);
// # endif
// #endif
*/
#if defined(__cpp_lib_format)
#   if __cpp_lib_format >= 202207L
      printf("   __cpp_lib_format = %ld (C++23, depends on <format>)\n",
          (long) __cpp_lib_format);
#   elif __cpp_lib_format >= 202106L
      printf("   __cpp_lib_format = %ld (C++20 DR, depends on <format>)\n",
          (long) __cpp_lib_format);
#   else
      printf("   __cpp_lib_format = %ld (C++20, depends on <format>)\n",
          (long) __cpp_lib_format);
#   endif
#endif
#if defined(__cpp_lib_generic_unordered_lookup)
      printf("   __cpp_lib_generic_unordered_lookup = %ld (depends on <unordered_map>, <unordered_set>)\n",
          (long) __cpp_lib_generic_unordered_lookup);
#endif
#if defined(__cpp_lib_interpolate)
      printf("   __cpp_lib_interpolate = %ld (depends on <cmath>, <numeric>)\n",
          (long) __cpp_lib_interpolate);
#endif
#if defined(__cpp_lib_is_constant_evaluated)
      printf("   __cpp_lib_is_constant_evaluated = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_is_constant_evaluated);
#endif
#if defined(__cpp_lib_is_layout_compatible)
      printf("   __cpp_lib_is_layout_compatible = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_is_layout_compatible);
#endif
#if defined(__cpp_lib_is_pointer_interconvertible)
      printf("   __cpp_lib_is_pointer_interconvertible = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_is_pointer_interconvertible);
#endif
#if defined(__cpp_lib_jthread)
      printf("   __cpp_lib_jthread = %ld (depends on <stop_token> <thread>)\n",
          (long) __cpp_lib_jthread);
#endif
#if defined(__cpp_lib_latch)
      printf("   __cpp_lib_latch = %ld (depends on <latch>)\n",
          (long) __cpp_lib_latch);
#endif
#if defined(__cpp_lib_list_remove_return_type)
      printf("   __cpp_lib_list_remove_return_type = %ld (depends on <forward_list> <list>)\n",
          (long) __cpp_lib_list_remove_return_type);
#endif
#if defined(__cpp_lib_math_constants)
      printf("   __cpp_lib_math_constants = %ld (depends on <numbers>)\n",
          (long) __cpp_lib_math_constants);
#endif
#if defined(__cpp_lib_ranges)
#  if __cpp_lib_ranges >= 202202L
      printf("   __cpp_lib_ranges = %ld (C++23, depends on <algorithm> ... <ranges>)\n",
          (long) __cpp_lib_ranges);
#  elif __cpp_lib_ranges >= 202106L
      printf("   __cpp_lib_ranges = %ld (C++20 DR, depends on <algorithm> ... <ranges>)\n",
          (long) __cpp_lib_ranges);
#  else
      printf("   __cpp_lib_ranges = %ld (C++20, depends on <algorithm> ... <ranges>)\n",
          (long) __cpp_lib_ranges);
#  endif
#endif
#if defined(__cpp_lib_semaphore)
      printf("   __cpp_lib_semaphore = %ld (depends on <semaphore>)\n",
          (long) __cpp_lib_semaphore);
#endif
#if defined(__cpp_lib_source_location)
      printf("   __cpp_lib_source_location = %ld (depends on <source_location>)\n",
          (long) __cpp_lib_source_location);
#endif
#if defined(__cpp_lib_string_view)
# if __cpp_lib_string_view > 201606
      printf("   __cpp_lib_string_view = %ld (depends on )\n",
          (long) __cpp_lib_string_view);
# endif
#endif
#if defined(__cpp_lib_spaceship)
      printf("   __cpp_lib_spaceship = %ld (depends on <compare>)\n",
          (long) __cpp_lib_spaceship);
#endif
#if defined(__cpp_lib_three_way_comparison)
      printf("   __cpp_lib_three_way_comparison = %ld (depends on <compare>)\n",
          (long) __cpp_lib_three_way_comparison);
#endif
#if defined(__cpp_lib_to_array)
      printf("   __cpp_lib_to_array = %ld (depends on <array>)\n",
          (long) __cpp_lib_to_array);
#endif
#if defined(__cpp_lib_constexpr_complex)
      printf("   __cpp_lib_constexpr_complex = %ld (depends on <complex>)\n",
          (long) __cpp_lib_constexpr_complex);
#endif

   /* ----------- C++23 features ------------ */

   printf("  C++23 features available:\n");

#  if __has_cpp_attribute(assume)
      printf("   __has_cpp_attribute(assume) = %ld\n",
          (long) __has_cpp_attribute(assume));
#  endif
#if defined(__cpp_auto_cast)
      printf("   __cpp_auto_cast = %ld\n",
          (long) __cpp_auto_cast);
#endif
#if defined(__cpp_explicit_this_parameter)
      printf("   __cpp_explicit_this_parameter = %ld\n",
          (long) __cpp_explicit_this_parameter);
#endif
#if defined(__cpp_if_consteval)
      printf("   __cpp_if_consteval = %ld\n",
          (long) __cpp_if_consteval);
#endif
#if defined(__cpp_implicit_move)
      printf("   __cpp_implicit_move = %ld\n",
          (long) __cpp_implicit_move);
#endif
#if defined(__cpp_multidimensional_subscript)
      printf("   __cpp_multidimensional_subscript = %ld\n",
          (long) __cpp_multidimensional_subscript);
#endif
#if defined(__cpp_named_character_escapes)
      printf("   __cpp_named_character_escapes = %ld\n",
          (long) __cpp_named_character_escapes);
#endif
#if defined(__cpp_size_t_suffix)
      printf("   __cpp_size_t_suffix = %ld\n",
          (long) __cpp_size_t_suffix);
#endif
#if defined(__cpp_static_call_operator)
      printf("   __cpp_static_call_operator = %ld\n",
          (long) __cpp_static_call_operator);
#endif

#if defined(__cpp_lib_adaptor_iterator_pair_constructor)
      printf("   __cpp_lib_adaptor_iterator_pair_constructor  = %ld (depends on <queue> <stack>)\n",
          (long) __cpp_lib_adaptor_iterator_pair_constructor );
#endif
#if defined( __cpp_lib_algorithm_iterator_requirements )
      printf("   __cpp_lib_algorithm_iterator_requirements  = %ld (depends on <algorithm> <memory> <numeric>)\n",
          (long) __cpp_lib_algorithm_iterator_requirements );
#endif
#if defined(__cpp_lib_allocate_at_least)
      printf("   __cpp_lib_allocate_at_least = %ld (depends on <memory>)\n",
          (long) __cpp_lib_allocate_at_least);
#endif
#if defined(__cpp_lib_associative_heterogeneous_erasure)
      printf("   __cpp_lib_associative_heterogeneous_erasure = %ld (depends on <map> <set> <unordered_map> <unordered_set>)\n",
          (long) __cpp_lib_associative_heterogeneous_erasure);
#endif
#if defined(__cpp_lib_barrier)
      printf("   __cpp_lib_barrier = %ld (depends on <barrier>)\n",
          (long) __cpp_lib_barrier);
#endif
#if defined(__cpp_lib_bind_back)
      printf("   __cpp_lib_bind_back = %ld (depends on <functional>)\n",
          (long) __cpp_lib_bind_back);
#endif
#if defined(__cpp_lib_byteswap)
      printf("   __cpp_lib_byteswap = %ld (depends on <bit>)\n",
          (long) __cpp_lib_byteswap);
#endif
#if defined(__cpp_lib_common_reference)
      printf("   __cpp_lib_common_reference = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_common_reference);
#endif
#if defined(__cpp_lib_common_reference_wrapper)
      printf("   __cpp_lib_common_reference_wrapper = %ld (depends on <functional>)\n",
          (long) __cpp_lib_common_reference_wrapper);
#endif
#if defined(__cpp_lib_constexpr_bitset)
      printf("   __cpp_lib_constexpr_bitset = %ld (depends on <bitset>)\n",
          (long) __cpp_lib_constexpr_bitset);
#endif
#if defined(__cpp_lib_constexpr_charconv)
      printf("   __cpp_lib_constexpr_charconv = %ld (depends on <charconv>)\n",
          (long) __cpp_lib_constexpr_charconv);
#endif
#if defined(__cpp_lib_constexpr_cmath)
      printf("   __cpp_lib_constexpr_cmath = %ld (depends on <cmath> <cstdlib>)\n",
          (long) __cpp_lib_constexpr_cmath);
#endif
#if defined(__cpp_lib_constexpr_memory)
      printf("   __cpp_lib_constexpr_memory = %ld (depends on <memory>)\n",
          (long) __cpp_lib_constexpr_memory);
#endif
#if defined(__cpp_lib_containers_ranges)
      printf("   __cpp_lib_containers_ranges = %ld (depends on <vector> ...)\n",
          (long) __cpp_lib_containers_ranges);
#endif
#if defined(__cpp_lib_expected)
      printf("   __cpp_lib_expected = %ld (depends on <expected>)\n",
          (long) __cpp_lib_expected);
#endif
#if defined(__cpp_lib_flat_map)
      printf("   __cpp_lib_flat_map = %ld (depends on <flat_map>)\n",
          (long) __cpp_lib_flat_map);
#endif
#if defined(__cpp_lib_flat_set)
      printf("   __cpp_lib_flat_set = %ld (depends on <flat_set>)\n",
          (long) __cpp_lib_flat_set);
#endif
#if defined(__cpp_lib_format_ranges)
      printf("   __cpp_lib_format_ranges = %ld (depends on <format>)\n",
          (long) __cpp_lib_format_ranges);
#endif
#if defined(__cpp_lib_formatters)
      printf("   __cpp_lib_formatters = %ld (depends on <stacktrace> <thread>)\n",
          (long) __cpp_lib_formatters);
#endif
#if defined(__cpp_lib_forward_like)
      printf("   __cpp_lib_forward_like = %ld (depends on <utility>)\n",
          (long) __cpp_lib_forward_like);
#endif
#if defined(__cpp_lib_generator)
      printf("   __cpp_lib_generator = %ld (depends on <generator>)\n",
          (long) __cpp_lib_generator);
#endif
#if defined(__cpp_lib_invoke_r)
      printf("   __cpp_lib_invoke_r = %ld (depends on <functional>)\n",
          (long) __cpp_lib_invoke_r);
#endif
#if defined(__cpp_lib_ios_noreplace)
      printf("   __cpp_lib_ios_noreplace = %ld (depends on <ios>)\n",
          (long) __cpp_lib_ios_noreplace);
#endif
#if defined(__cpp_lib_is_implicit_lifetime)
      printf("   __cpp_lib_is_implicit_lifetime = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_is_implicit_lifetime);
#endif
#if defined(__cpp_lib_is_scoped_enum)
      printf("   __cpp_lib_is_scoped_enum = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_is_scoped_enum);
#endif
#if defined(__cpp_lib_mdspan)
      printf("   __cpp_lib_mdspan = %ld (depends on <mdspan>)\n",
          (long) __cpp_lib_mdspan);
#endif
#if defined(__cpp_lib_modules)
      printf("   __cpp_lib_modules = %ld (standard library modules std and std.compat)\n",
          (long) __cpp_lib_modules);
#endif
#if defined(__cpp_lib_move_iterator_concept)
      printf("   __cpp_lib_move_iterator_concept = %ld (depends on <iterator>)\n",
          (long) __cpp_lib_move_iterator_concept);
#endif
#if defined(__cpp_lib_move_only_function)
      printf("   __cpp_lib_move_only_function = %ld (depends on <functional>)\n",
          (long) __cpp_lib_move_only_function);
#endif
#if defined(__cpp_lib_out_ptr)
      printf("   __cpp_lib_out_ptr = %ld (depends on <memory>)\n",
          (long) __cpp_lib_out_ptr);
#endif
#if defined(__cpp_lib_print)
      printf("   __cpp_lib_print = %ld (depends on <ostream> <print>)\n",
          (long) __cpp_lib_print);
#endif
#if defined(__cpp_lib_ranges_as_const)
      printf("   __cpp_lib_ranges_as_const = %ld (depends on <ranges>)\n",
          (long) __cpp_lib_ranges_as_const);
#endif
#if defined(__cpp_lib_ranges_as_rvalue)
      printf("   __cpp_lib_ranges_as_rvalue = %ld (depends on <ranges>)\n",
          (long) __cpp_lib_ranges_as_rvalue);
#endif
#if defined(__cpp_lib_ranges_cartesian_product)
      printf("   __cpp_lib_ranges_cartesian_product = %ld (depends on <ranges>)\n",
          (long) __cpp_lib_ranges_cartesian_product);
#endif
#if defined(__cpp_lib_ranges_chunk)
      printf("   __cpp_lib_ranges_chunk = %ld (depends on <ranges>)\n",
          (long) __cpp_lib_ranges_chunk);
#endif
#if defined(__cpp_lib_ranges_chunk_by)
      printf("   __cpp_lib_ranges_chunk_by = %ld (depends on <ranges>)\n",
          (long) __cpp_lib_ranges_chunk_by);
#endif
#if defined(__cpp_lib_ranges_contains)
      printf("   __cpp_lib_ranges_contains = %ld (depends on <algorithm>)\n",
          (long) __cpp_lib_ranges_contains);
#endif
#if defined(__cpp_lib_ranges_enumerate)
      printf("   __cpp_lib_ranges_enumerate = %ld (depends on <ranges>)\n",
          (long) __cpp_lib_ranges_enumerate);
#endif
#if defined(__cpp_lib_ranges_find_last)
      printf("   __cpp_lib_ranges_find_last = %ld (depends on <algorithm>)\n",
          (long) __cpp_lib_ranges_find_last);
#endif
#if defined(__cpp_lib_ranges_fold)
      printf("   __cpp_lib_ranges_fold = %ld (depends on <algorithm>)\n",
          (long) __cpp_lib_ranges_fold);
#endif
#if defined(__cpp_lib_ranges_iota)
      printf("   __cpp_lib_ranges_iota = %ld (depends on <numeric>)\n",
          (long) __cpp_lib_ranges_iota);
#endif
#if defined(__cpp_lib_ranges_join_with)
      printf("   __cpp_lib_ranges_join_with = %ld (depends on <ranges>)\n",
          (long) __cpp_lib_ranges_join_with);
#endif
#if defined(__cpp_lib_ranges_repeat)
      printf("   __cpp_lib_ranges_repeat = %ld (depends on <ranges>)\n",
          (long) __cpp_lib_ranges_repeat);
#endif
#if defined(__cpp_lib_ranges_slide)
      printf("   __cpp_lib_ranges_slide = %ld (depends on <ranges>)\n",
          (long) __cpp_lib_ranges_slide);
#endif
#if defined(__cpp_lib_ranges_starts_ends_with)
      printf("   __cpp_lib_ranges_starts_ends_with = %ld (depends on <algorithm>)\n",
          (long) __cpp_lib_ranges_starts_ends_with);
#endif
#if defined(__cpp_lib_ranges_stride)
      printf("   __cpp_lib_ranges_stride = %ld (depends on <ranges>)\n",
          (long) __cpp_lib_ranges_stride);
#endif
#if defined(__cpp_lib_ranges_to_container)
      printf("   __cpp_lib_ranges_to_container = %ld (depends on <ranges>)\n",
          (long) __cpp_lib_ranges_to_container);
#endif
#if defined(__cpp_lib_ranges_zip)
      printf("   __cpp_lib_ranges_zip = %ld (depends on <ranges> <tuple> <utility>)\n",
          (long) __cpp_lib_ranges_zip);
#endif
#if defined(__cpp_lib_reference_from_temporary)
      printf("   __cpp_lib_reference_from_temporary = %ld (depends on <type_traits>)\n",
          (long) __cpp_lib_reference_from_temporary);
#endif
#if defined(__cpp_lib_shift)
      printf("   __cpp_lib_shift = %ld (depends on <algorithm>)\n",
          (long) __cpp_lib_shift);
#endif
#if defined(__cpp_lib_spanstream)
      printf("   __cpp_lib_spanstream = %ld (depends on <spanstream>)\n",
          (long) __cpp_lib_spanstream);
#endif
#if defined(__cpp_lib_stacktrace)
      printf("   __cpp_lib_stacktrace = %ld (depends on <stacktrace>)\n",
          (long) __cpp_lib_stacktrace);
#endif
#if defined(__cpp_lib_start_lifetime_as)
      printf("   __cpp_lib_start_lifetime_as = %ld (depends on <memory>)\n",
          (long) __cpp_lib_start_lifetime_as);
#endif
#if defined(__cpp_lib_stdatomic_h)
      printf("   __cpp_lib_stdatomic_h = %ld (depends on <stdatomic.h>)\n",
          (long) __cpp_lib_stdatomic_h);
#endif
#if defined(__cpp_lib_string_contains)
      printf("   __cpp_lib_string_contains = %ld (depends on <string> <string_view>)\n",
          (long) __cpp_lib_string_contains);
#endif
#if defined(__cpp_lib_string_resize_and_overwrite)
      printf("   __cpp_lib_string_resize_and_overwrite = %ld (depends on <string>)\n",
          (long) __cpp_lib_string_resize_and_overwrite);
#endif
#if defined(__cpp_lib_to_underlying)
      printf("   __cpp_lib_to_underlying = %ld (depends on <utility>)\n",
          (long) __cpp_lib_to_underlying);
#endif
#if defined(__cpp_lib_tuple_like)
      printf("   __cpp_lib_tuple_like = %ld (depends on <map> <tuple> <unordered_map> <utility>)\n",
          (long) __cpp_lib_tuple_like);
#endif
#if defined(__cpp_lib_unreachable)
      printf("   __cpp_lib_unreachable = %ld (depends on <utility>)\n",
          (long) __cpp_lib_unreachable);
#endif

 
#ifndef __cplusplus
   printf("\nCompiled by C compiler instead of C++ compiler.\n");
# ifndef NO_C99_CHECK
#  ifndef ASSUME_C99
#   define NO_C99 1
#  endif
#  ifdef __STDC_VERSION__
#   if __STDC_VERSION__ >= 199901L
#    undef NO_C99
#   endif
#  else
#   ifdef __USE_ISOC95
#    define GNU90
#   endif
#  endif
#  ifdef GNU90
   // C++ style comments known to work with -std=gnu90
   printf("C compiler apparently working in -std=gnu90 mode.\n");
#  endif
#  ifndef NO_C99
   // Test C++ style comments also when compiling by C compiler
   // #if Test with C++ style comment hiding preprocessor line
   /* Check for loop variables declared in the loop: */
   printf("Testing C99 mode with 'for ( int i=0; i<4; i++ )': ");
   for ( int i=0; i<4; i++ )
      printf("%d ",i);
   printf(" (OK)\n");
#  else
   printf("C compiler not claiming to be C99 compatible.\n");
#  endif
# endif
#endif

   return 0;
}
