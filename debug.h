// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#ifndef _DEBUG_H
#define _DEBUG_H

// Debug-only checks of internal invariants, compiled into the builds that define RN_DEBUG_CHECKS (CMakeLists.txt: the Debug
// configuration, and rn_obj_search_stats, whose times are never used). A failing check prints its expression and where it failed,
// then aborts. In the Release rn_obj_search and rn_obj_search_par it expands to nothing and its expression is not evaluated. Unlike assert,
// it does not depend on NDEBUG
#if defined(RN_DEBUG_CHECKS)
#include <cstdio>
#include <cstdlib>
namespace debug_checks
{
    [[noreturn]] inline void failed(const char* expression, const char* file, int line, const char* function)
    {
        std::fprintf(stderr, "%s:%d: %s: Debug check `%s' failed.\n", file, line, function, expression);
        std::abort();
    }
}
#define RN_DEBUG_ASSERT(...) ((__VA_ARGS__) ? static_cast<void>(0) : debug_checks::failed(#__VA_ARGS__, __FILE__, __LINE__, __PRETTY_FUNCTION__))
#else
#define RN_DEBUG_ASSERT(...) static_cast<void>(0)
#endif

#endif // _DEBUG_H
