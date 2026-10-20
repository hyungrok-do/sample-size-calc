#include "util/memory.hpp"

#if defined(__GLIBC__)
#include <malloc.h>
#elif defined(__APPLE__)
#include <malloc/malloc.h>
#elif defined(_WIN32)
#include <malloc.h>
#endif

namespace ssc::memory {

void release_free_memory() {
#if defined(__GLIBC__)
  malloc_trim(0);
#elif defined(__APPLE__)
  malloc_zone_pressure_relief(nullptr, 0);
#elif defined(_WIN32)
  _heapmin();
#endif
}

}  // namespace ssc::memory
