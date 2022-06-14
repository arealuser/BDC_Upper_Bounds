#include "cached_transition_probs.h"
#include "cache_io.h"
#include <cassert>

std::array<std::array<std::vector<Int>, MAX_BIT_CACHE_SIZE>, MAX_BIT_CACHE_SIZE> cached_transition_probs;

void load_transition_cache(size_t n, size_t k){
	assert(n < MAX_BIT_CACHE_SIZE);
	assert(k < MAX_BIT_CACHE_SIZE);

	cached_transition_probs[n][k] = load_data_from_cache_file(n, k);
	// printf("Loaded cache %lu, %lu. It is now of length %lu.\n", n, k, cached_transition_probs[n][k].size());
}