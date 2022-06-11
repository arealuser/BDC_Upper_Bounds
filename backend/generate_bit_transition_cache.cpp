#include "bit_channel.h"

/*
Given an input length n, an output length k and an output file, produces the transition counts of all of the codewords
	of length n (2 ** n in total) to all the codewords of length at most k ((2 ** (k+1)) - 1 in total) and saves the
	results to the given output as a CPP code.
*/
int main(int argc, char const *argv[])
{
	if (argc != 6)
	{
		fprintf(stderr, "Usage %s n k output.cc output_path.h output_filename.h\n", argv[0]);
		return 1;
	}

	size_t n = atol(argv[1]);
	size_t k = atol(argv[2]);
	const char* cpp_filename = argv[3];
	const char* header_path = argv[4];
	const char* header_filename = argv[5];
	FILE* cpp_file = try_to_open_file(cpp_filename, "wb");
	FILE* header_file = try_to_open_file(header_path, "wb");

	auto transmitted_codewords = get_all_bit_codewords(n, false);
	std::sort(transmitted_codewords.begin(), transmitted_codewords.end(),
		[](const BitCodeWord& word1, BitCodeWord& word2){return btc_to_idx(word1) < btc_to_idx(word2);});

	auto receveived_codewords = get_all_bit_codewords(k, false);
	std::sort(receveived_codewords.begin(), receveived_codewords.end(),
		[](const BitCodeWord& word1, BitCodeWord& word2){return btc_to_idx(word1) < btc_to_idx(word2);});
	

	fprintf(cpp_file, "#include \"%s\"\nstd::array<std::array<size_t, %lu>, %lu> bit_possibility_cache_%lu_%lu = {", header_filename, 
		receveived_codewords.size(), transmitted_codewords.size(), n, k);
	for (size_t i = 0; i < transmitted_codewords.size(); ++i)
	{
		fprintf(cpp_file, "\n\t");
		for (size_t j = 0; j < receveived_codewords.size(); ++j)
		{
			size_t count = get_num_transition_possibilities(transmitted_codewords[i], receveived_codewords[j], false);
			fprintf(cpp_file, "0x%lxLL", count);
			if ((i < transmitted_codewords.size()-1) or (j < receveived_codewords.size()-1))
			{
				fprintf(cpp_file, ",");
			}
		}
	}
	fprintf(cpp_file, "\n};\n");

	fprintf(header_file, "#pragma once\n#include<cstdlib>\n#include <array>\nextern std::array<std::array<size_t, %lu>, %lu> bit_possibility_cache_%lu_%lu;\n",
		receveived_codewords.size(), transmitted_codewords.size(), n, k);

	fprintf(header_file, "inline size_t get_transition_count_%lu_%lu(size_t trans_num, size_t rec_num){\n", n, k);
	fprintf(header_file, "\treturn bit_possibility_cache_%lu_%lu[trans_num][rec_num];\n", n, k);
	fprintf(header_file, "}\n");

	fclose(cpp_file);
	fclose(header_file);

	return 0;
}