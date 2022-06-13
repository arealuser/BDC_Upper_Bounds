#pragma once
#include "bit_channel.h"
#include "transition_counts/cached_transition_count.h"



extern std::vector<std::vector<Float> > _normalization_factors;
extern Float _deletion_prob;

inline size_t get_num_transition_possibilities(const BitCodeWord& transmitted, const BitCodeWord& recieved, bool verbose){
	size_t st = transmitted.size() + 1;
	size_t sr = recieved.size() + 1;

	std::vector<size_t> dynamic_programming_state;
	dynamic_programming_state.resize(st * sr);

	for (size_t it = 0; it < st; ++it)
	{
		dynamic_programming_state[it] = 1;
	}
	for (size_t ir = 1; ir < sr; ++ir)
	{
		for (size_t it = 0; it < st; ++it)
		{
			if (it == 0)
			{
				continue;
			}
			if (it + (sr - ir) < sr)
			{
				continue;
			}

			dynamic_programming_state[(ir*st)+it] = dynamic_programming_state[(ir*st)+it-1];
			if (transmitted[it - 1] == recieved[ir - 1])
			{
				dynamic_programming_state[(ir * st) + it] += dynamic_programming_state[((ir-1) * st) + (it-1)];
			}
		}
	}
	if (verbose)
	{
		printf("dynamic_programming_state = [\n");
		for (size_t i = 0; i < st; ++i)
		{
			printf("\t[");
			for (size_t j = 0; j < sr; ++j)
			{
				printf("%lu\t", dynamic_programming_state[(i*st)+j]);
			}
			printf("]\n");
		}
		printf("]\n");
	}

	return dynamic_programming_state[(sr * st) - 1];
}


inline Float get_bit_transition_prob(const BitCodeWord& transmitted, const BitCodeWord& recieved, bool verbose, bool use_cache){
	size_t st = transmitted.size() + 1;
	size_t sr = recieved.size() + 1;

	Float base_prob = _normalization_factors[st - 1][sr - 1];
	
	size_t count;
	if (verbose)
	{
		printf("use_cache=%d\n", use_cache);
	}
	if (use_cache)
	{
		count = get_num_transition_possibilities_using_cache(transmitted, recieved, verbose);
	} else {
		count = get_num_transition_possibilities(transmitted, recieved, verbose);
	}
	if (verbose)
	{
		printf("%lu, %lu\n", transmitted.size(), recieved.size());
		printf("transmitted = [");
		for (size_t i = 0; i < transmitted.size(); ++i)
		{
			printf("%d, ", transmitted[i]);
		}
		printf("]\n");

		printf("recieved = [");
		for (size_t i = 0; i < recieved.size(); ++i)
		{
			printf("%d, ", recieved[i]);
		}
		printf("]\n");
		printf("base_prob = %f, count = %lu, final_prob = %f\n", 
			base_prob, count, count * base_prob);
	}
	return base_prob * count;
}


inline Float get_bit_transition_prob_fast(const EfficientBitCodeWord& transmitted, const EfficientBitCodeWord& recieved, bool verbose){
	size_t st = transmitted.len + 1;
	size_t sr = recieved.len + 1;

	Float base_prob = _normalization_factors[st - 1][sr - 1];
	
	size_t count;
	count = get_num_transition_possibilities_using_cache_fast(transmitted, recieved, verbose);
	if (verbose)
	{
		printf("%lu, %lu\n", transmitted.len, recieved.len);
		printf("transmitted = %016lx\n", num_to_idx(transmitted.num, transmitted.len));

		printf("recieved = %016lx\n", num_to_idx(recieved.num, recieved.len));
		printf("base_prob = %f, count = %lu, final_prob = %f\n", 
			base_prob, count, count * base_prob);
	}
	return base_prob * count;
}


inline size_t get_num_transition_possibilities_using_cache(const BitCodeWord& transmitted, const BitCodeWord& recieved, bool verbose){
	size_t n = transmitted.size();
	size_t k = recieved.size();
	size_t n1 = (n+1) / 2;
	size_t n2 = n - n1;
	size_t total = 0;
	size_t trans_num1 = btc_to_num(transmitted.begin(), transmitted.begin() + n1);
	size_t trans_num2 = btc_to_num(transmitted.begin() + n1, transmitted.end());
	for (size_t i = 0; i <= k; ++i)
	{
		size_t k1 = i;
		size_t k2 = k - i;
		if ((k1 > n1) or (k2 > n2))
		{
			if (verbose)
			{
				printf("Skipping n1=%lu\tk1=%lu\tn2=%lu\tk2=%lu\n", n1, k1, n2, k2);
			}
			continue;
		}
		size_t rec_num1 = btc_to_num(recieved.begin(), recieved.begin() + k1);
		size_t rec_num2 = btc_to_num(recieved.begin() + k1, recieved.end());
		if (verbose)
		{
			printf("n1=%lu\t k1=%lu\t n2=%lu\t k2=%lu\t trans_num1=%lu\t trans_num2=%lu\t rec_num1=%lu\t rec_num2=%lu\n", n1, k1, n2, k2,
				trans_num1, trans_num2, rec_num1, rec_num2);
			fflush(stdout);
		}
		size_t count1 = get_transition_count_cache(n1, trans_num1, k1, rec_num1);
		size_t count2 = get_transition_count_cache(n2, trans_num2, k2, rec_num2);
		total += count1*count2;
		if (verbose)
		{
			printf("n=%lu\tk=%lu\tcount1=%lu\tcount2=%lutotal=%lu\n", n, k, count1, count2, total);
		}
	}

	return total;
}

inline size_t get_num_transition_possibilities_using_cache_fast(const EfficientBitCodeWord& transmitted, const EfficientBitCodeWord& recieved, bool verbose){
	size_t n = transmitted.len;
	size_t k = recieved.len;
	size_t n1 = (n+1) / 2;
	size_t n2 = n - n1;
	size_t total = 0;
	size_t trans_num1 = transmitted.num >> n2;
	size_t trans_num2 = transmitted.num  & ((1 << n2) - 1);
	for (size_t i = 0; i <= k; ++i)
	{
		size_t k1 = i;
		size_t k2 = k - i;
		if ((k1 > n1) or (k2 > n2))
		{
			if (verbose)
			{
				printf("Skipping n1=%lu\tk1=%lu\tn2=%lu\tk2=%lu\n", n1, k1, n2, k2);
			}
			continue;
		}
		size_t rec_num1 = recieved.num >> k2;
		size_t rec_num2 = recieved.num & ((1 << k2) - 1);
		if (verbose)
		{
			printf("n1=%lu\t k1=%lu\t n2=%lu\t k2=%lu\t trans_num1=%lu\t trans_num2=%lu\t rec_num1=%lu\t rec_num2=%lu\n", n1, k1, n2, k2,
				trans_num1, trans_num2, rec_num1, rec_num2);
			fflush(stdout);
		}
		size_t count1 = get_transition_count_cache(n1, trans_num1, k1, rec_num1);
		size_t count2 = get_transition_count_cache(n2, trans_num2, k2, rec_num2);
		total += count1*count2;
		if (verbose)
		{
			printf("count1=%lu\tcount2=%lu\ttotal=%lu\n", count1, count2, total);
		}
	}

	return total;
}

inline uint64_t btc_to_num(const BitCodeWord& codeword){
	uint64_t res = 0;
	for (size_t i = 0; i < codeword.size(); ++i)
	{
		res <<= 1;
		res ^= codeword[i];
	}
	return res;
}
inline BitCodeWord num_to_btc(uint64_t num, size_t len){
	BitCodeWord res;
	res.reserve(len);
	for (size_t i = 0; i < len; ++i)
	{
		res.push_back((num >> i) & 0x01);
	}
	return res;
}

inline uint64_t btc_to_idx(const BitCodeWord& codeword){
	return num_to_idx(btc_to_num(codeword), codeword.size());
}

inline uint64_t btc_to_idx(const EfficientBitCodeWord& codeword){
	return num_to_idx(codeword.num, codeword.len);
}

inline uint64_t num_to_idx(uint64_t num, size_t len){
	return num ^ (1ULL << len);
}


template <typename _InputIter>
inline uint64_t btc_to_num(_InputIter first, _InputIter last){
	if (first == last)
	{
		return 0;
	}
	size_t res = 0;
	for(; first != last; ++first){
		res <<= 1;
		res ^= *first;
	}
	return res;
}

inline bool operator< (const EfficientBitCodeWord& a, const EfficientBitCodeWord& b){
	if (a.len != b.len)
	{
		return a.len < b.len;
	}
	uint64_t n1 = std::min(a.num, (uint64_t)((1ULL << a.len) - 1) ^ a.num);
	uint64_t n2 = std::min(b.num, (uint64_t)((1ULL << b.len) - 1) ^ b.num);
	if (n1 != n2)
	{
		return n1 < n2;
	}
	return a.num < b.num;
}


inline EfficientBitCodeWord operator~ (const EfficientBitCodeWord& a){
	return EfficientBitCodeWord(~a.num, a.len);
}
