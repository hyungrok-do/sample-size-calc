/*
 * Sample Size Calculator — C interface to the C++ computational core.
 *
 * The interface is deliberately small and JSON based so that any language
 * with a C FFI (Python ctypes, R, Julia, C#, ...) can drive the engine, and so
 * that one compiled library works with every Python 3 version.
 *
 * Memory: every `char*` returned by this library must be released with
 * ssc_free().  Strings passed in are UTF-8 and only borrowed for the call.
 * ssc_catalog_json, ssc_procedure_json and ssc_search_json return NULL when
 * there is not enough memory for their answer (ssc_procedure_json also for an
 * unknown id); ssc_run_json answers that it ran out of memory, and returns
 * NULL only when there is not enough memory even for that short answer.
 *
 * Threading: all functions are thread-safe.  Long computations can be
 * monitored and cancelled from another thread through an ssc_token.
 */
#ifndef SSC_SSC_H
#define SSC_SSC_H

#if defined(SSC_STATIC)
/* linked directly into an executable (the bundled CLI and tests) */
#  define SSC_API
#elif defined(_WIN32) || defined(__CYGWIN__)
#  ifdef SSC_BUILDING
#    define SSC_API __declspec(dllexport)
#  else
#    define SSC_API __declspec(dllimport)
#  endif
#else
#  define SSC_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ssc_token ssc_token;

/** Library version string (static storage, do not free). */
SSC_API const char* ssc_version(void);

/** Number of hardware threads available (>= 1). */
SSC_API int ssc_hardware_threads(void);

/** JSON array describing every procedure (id, name, category, ...). */
SSC_API char* ssc_catalog_json(void);

/** Full JSON schema of one procedure (parameters, solve-for options, ...).
 *  Returns NULL if the id is unknown, or if there is not enough memory for
 *  the schema. */
SSC_API char* ssc_procedure_json(const char* procedure_id);

/**
 * Search the catalog: procedures ranked for a free-text query, tolerating typos,
 * plurals, hyphenation, abbreviations and synonyms.  Returns
 *   {"query": "...", "corrected": "...", "partial": false,
 *    "hits": [{"id": "...", "name": "...", "category": "...", "score": 12.3, "terms": [...]}]}
 * where `corrected` is the query with misspellings fixed ("" if none) and `partial`
 * is true when no procedure matches every word.  `limit` <= 0 returns every match.
 */
SSC_API char* ssc_search_json(const char* query, int limit);

/**
 * Run a calculation.  `request_json` has the form
 *   {"procedure": "...", "solve_for": "...", "inputs": {"id": "text", ...},
 *    "options": {"threads": 0, "power_curve": true, ...}}
 * and the result is a JSON object (see docs/api.md).  `token` may be NULL.
 * A calculation that runs out of memory answers
 *   {"ok": false, "error": "There is not enough memory for this calculation: ..."};
 * NULL is returned only when there is not enough memory even for a short
 * answer such as that one.
 */
SSC_API char* ssc_run_json(const char* request_json, ssc_token* token);

/** DAG simulation and estimand-linked power. JSON operations: simulate, plan,
 * population, calibrate, analyze, nullGraph, power. Uses the calculator's RNG,
 * statistical core and thread pool. Returns {ok:true,result:...} or
 * {ok:false,error:...}; release the returned string with ssc_free(). */
SSC_API char* ssc_dag_json(const char* request_json, ssc_token* token);

/** Release a string returned by this library. */
SSC_API void ssc_free(char* p);

/** Cancellation / progress token. */
SSC_API ssc_token* ssc_token_new(void);
SSC_API void ssc_token_cancel(ssc_token* token);
/** Fraction of work completed in [0, 1]. */
SSC_API double ssc_token_progress(const ssc_token* token);
SSC_API void ssc_token_free(ssc_token* token);

/**
 * Evaluate a distribution / special function by name, e.g.
 *   "pnorm"(x, lower), "qt"(p, df, lower), "pnt"(t, df, ncp, lower).
 * Returns 0 on success (result in *out), non-zero for an unknown name or a
 * wrong number of arguments.  Mainly intended for testing and scripting.
 */
SSC_API int ssc_eval(const char* name, const double* args, int nargs, double* out);

#ifdef __cplusplus
}
#endif

#endif /* SSC_SSC_H */
