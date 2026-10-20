// Returning freed memory to the operating system.
//
// After a simulation the allocator keeps the memory it handed out (replicates, model fits, worker threads'
// arenas) for reuse, so the process stays as large as at its peak.  release_free_memory() asks the C library
// to give the free parts back: malloc_trim on glibc, the zone pressure relief on macOS, _heapmin on Windows.
#pragma once

namespace ssc::memory {

void release_free_memory();

}  // namespace ssc::memory
