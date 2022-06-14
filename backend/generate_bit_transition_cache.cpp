#include "bit_channel.h"
#include "cache_io.h"

/*
Given an input length n, an output length k and an output file, produces the transition counts of all of the codewords
	of length n (2 ** n in total) to all the codewords of length at most k ((2 ** (k+1)) - 1 in total) and saves the
	results to the given output as a CPP code.
*/
int main(int argc, char const *argv[])
{
	if (argc != 3)
	{
		fprintf(stderr, "Usage %s n k\n", argv[0]);
		return 1;
	}

	size_t n = atol(argv[1]);
	size_t k = atol(argv[2]);

	auto transmitted_codewords = get_all_bit_codewords(n, false);
	std::sort(transmitted_codewords.begin(), transmitted_codewords.end(),
		[](const BitCodeWord& word1, BitCodeWord& word2){return btc_to_idx(word1) < btc_to_idx(word2);});

	auto receveived_codewords = get_all_bit_codewords(k, false);
	std::sort(receveived_codewords.begin(), receveived_codewords.end(),
		[](const BitCodeWord& word1, BitCodeWord& word2){return btc_to_idx(word1) < btc_to_idx(word2);});

	std::vector<Int> data; data.reserve(transmitted_codewords.size() * receveived_codewords.size());

	for (int i = 0; i < transmitted_codewords.size(); ++i)
	{
		for (int j = 0; j < receveived_codewords.size(); ++j)
		{
			data.push_back(get_num_transition_possibilities(transmitted_codewords[i], receveived_codewords[j], false));
		}
	}
	save_data_to_cache_file(n, k, data);

	return 0;
}