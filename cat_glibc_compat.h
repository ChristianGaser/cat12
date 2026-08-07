/* cat_glibc_compat.h
 *
 * Keep the Linux mex-files loadable on an older glibc than the build machine's.
 *
 * An ELF object records, for every symbol it imports, the glibc symbol version
 * it was linked against, and the dynamic loader treats that as a hard minimum.
 * The mex-files are built on the oldest runner GitHub still offers, Ubuntu
 * 22.04 with glibc 2.35, so without this header MATLAB rejects them on any
 * older system with
 *
 *   Invalid MEX-file '.../cat_sanlm.mexa64':
 *   /lib64/libc.so.6: version `GLIBC_2.34' not found
 *
 * which covers RHEL/Rocky 8, Ubuntu 20.04 and Debian 10/11 -- still ordinary
 * cluster installations. Nothing in CAT needs anything new; the whole
 * incompatibility comes from two changes inside glibc itself:
 *
 *   pthread_create, pthread_join
 *       glibc 2.34 merged libpthread into libc, and symbols that had lived only
 *       in libpthread were given a fresh 2.34 version. The functions did not
 *       change. (pthread_mutex_* were always in libc and still bind at 2.2.5,
 *       which is why they are not listed here. pthread_exit binds at the
 *       baseline today too, but it came through the same merge and is pinned
 *       along with the other two rather than left to whichever default the
 *       linker happens to pick.)
 *
 *   exp, log, pow
 *       glibc 2.29 added more accurate implementations under a new version.
 *
 * glibc keeps the previous implementations as compat symbols indefinitely, so
 * asking for them by name costs nothing and is not a downgrade of anything CAT
 * relies on: these are the very functions every binary built before 2019 used.
 *
 * .symver rewrites the *reference*, so this header has to be seen by each
 * translation unit that makes one of these calls. That is mostly not the mex
 * entry points but the helper sources compiled alongside them -- Amap.c,
 * MrfPrior.c, vollib.c, ornlm_float.c, sanlm_float.c -- which is where the
 * arithmetic and the threading actually live. It is included from the sources
 * rather than injected through compile.m, so that mex and mkoctfile behave the
 * same and no flag change can quietly drop it.
 *
 * A source that starts calling one of these and does not include this header
 * raises the floor again for the whole mex-file. That is what the workflow
 * check is for: it is meant to be found there and not by a user.
 *
 * The floor drops from glibc 2.34 (2021) to glibc 2.14 (2011) -- what remains is
 * memcpy@GLIBC_2.14 -- which is below the glibc any supported MATLAB requires of
 * its own accord. .github/workflows/compile_mex.yml enforces that floor on every
 * build, because it cannot be noticed on the build machine.
 *
 * The failure mode is benign: a version named here that did not exist would fail
 * the link on the build machine. It cannot quietly emit a binary that is broken
 * somewhere else.
 *
 * Define CAT_NO_GLIBC_COMPAT to switch this off.
 */

#ifndef CAT_GLIBC_COMPAT_H
#define CAT_GLIBC_COMPAT_H

#if defined(__linux__) && !defined(CAT_NO_GLIBC_COMPAT)

/* __GLIBC__ is only defined once a libc header has been seen. */
#include <limits.h>

#if defined(__GLIBC__) && defined(__ELF__)

/* The baseline version node is the one that architecture's ABI started at, so
   it exists in every glibc able to run the binary at all. It is not the same
   everywhere: x86-64 entered at 2.2.5, aarch64 only at 2.17. On anything else
   the header does nothing, which is the safe direction -- guessing a baseline
   would turn into a link error rather than a fallback. */
#if defined(__x86_64__) && !defined(__ILP32__)
#define CAT_GLIBC_BASE "GLIBC_2.2.5"
#elif defined(__aarch64__)
#define CAT_GLIBC_BASE "GLIBC_2.17"
#endif

#ifdef CAT_GLIBC_BASE

/* A .symver for a symbol the translation unit never calls emits nothing at all,
   so this list costs a file only for the functions it really uses. */
__asm__(".symver pthread_create,pthread_create@" CAT_GLIBC_BASE);
__asm__(".symver pthread_join,pthread_join@" CAT_GLIBC_BASE);
__asm__(".symver pthread_exit,pthread_exit@" CAT_GLIBC_BASE);
__asm__(".symver exp,exp@" CAT_GLIBC_BASE);
__asm__(".symver log,log@" CAT_GLIBC_BASE);
__asm__(".symver pow,pow@" CAT_GLIBC_BASE);

#endif /* CAT_GLIBC_BASE */

#endif /* __GLIBC__ && __ELF__ */

#endif /* __linux__ && !CAT_NO_GLIBC_COMPAT */

#endif /* CAT_GLIBC_COMPAT_H */
