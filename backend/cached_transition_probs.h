#pragma once
#include <array>
#include "utils.h"


constexpr size_t MAX_BIT_CACHE_SIZE = 40;
extern std::array<std::array<std::vector<Int>, MAX_BIT_CACHE_SIZE>, MAX_BIT_CACHE_SIZE> cached_transition_probs;

/*
Loads the transition cache from words of length n to works of length k.
*/
void load_transition_cache(size_t n, size_t k);

/*
Returns the number of transitions from a word of length n and bitwise representation trans_num to
	a word of length k and bitwise representation rec_num
*/
inline size_t get_transition_count_cache(size_t n, size_t trans_num, size_t k, size_t rec_num){
	// printf("get_transition_count_cache(%lu, %lu, %lu, %lu)\n", n, trans_num, k, rec_num); fflush(stdout);
	const std::vector<Int>& cache = cached_transition_probs[n][k];
	// printf("cache.size() = %lu\n", cache.size()); fflush(stdout);
	size_t index = (trans_num << k) ^ rec_num;
	return cache[index];
}