#include "bit_channel.h"
#include <cmath>


Float get_bit_transition_prob(const BitCodeWord& transmitted, const BitCodeWord& recieved, Float deletion_prob, bool verbose){
	size_t st = transmitted.size() + 1;
	size_t sr = recieved.size() + 1;
	// printf("\n");

	std::vector<size_t> dynamic_programming_state;
	dynamic_programming_state.resize(st * sr);

	// printf("sr = %lu, st = %lu\n", sr, st);

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
				// dynamic_programming_state[(ir*st)+it] = 1;
				continue;
			}
			if (it + (sr - ir) < sr)
			{
				continue;
			}

			dynamic_programming_state[(ir*st)+it] = dynamic_programming_state[(ir*st)+it-1];
			if (transmitted[it - 1] == recieved[ir - 1])
			{
				// printf("ir = %lu, it = %lu\n", 
				// 	ir, it); fflush(stdout);
				dynamic_programming_state[(ir * st) + it] += dynamic_programming_state[((ir-1) * st) + (it-1)];
			}
		}
		// for (size_t i = 0; i < st; ++i)
		// {
		// 	printf("%lu\t", dynamic_programming_state[(ir * st) + i]);
		// }
		// printf("\n");
	}


	Float base_prob = pow(deletion_prob, transmitted.size() - recieved.size()) * pow(1 - deletion_prob, recieved.size());
	size_t count = dynamic_programming_state[(sr * st) - 1];
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
		printf("base_prob = %f * %f = %f, count = %lu, final_prob = %f\n", 
			pow(deletion_prob, transmitted.size() - recieved.size()),
			pow(1 - deletion_prob, recieved.size()), 
			base_prob, count, count * base_prob);
	}
	return base_prob * count;
}


CodeWord convert_to_run_word(const BitCodeWord& bit_code){
	std::vector<Run> res;
	if (bit_code.size() == 0)
	{
		return res;
	}

	res.push_back(Run(*bit_code.begin(), 0));
	for(auto iter = bit_code.begin(); iter != bit_code.end(); ++iter){
		if (*iter == res.rbegin() -> value)
		{
			res.rbegin() -> length++;
		} else{
			res.push_back(Run(*iter, 1));
		}
	}
	return res;
}


BitCodeWord convert_to_bit_word(const CodeWord& run_word){
	BitCodeWord res;
	for(auto iter = run_word.begin(); iter != run_word.end(); ++iter){
		for (size_t i = 0; i < iter -> length; ++i)
		{
			res.push_back(iter -> value);
		}
	}
	return res;
}
