#include "channel.h"
#include <cmath>

std::array<std::array<Float, MAX_RUN_LENGTH>, MAX_RUN_LENGTH> RUN_TO_RUN_PROB;
std::array<std::array<Float, MAX_RUN_LENGTH>, MAX_RUN_LENGTH> RUN_TO_RUN_PROB_TRANSPOSE;
std::array<std::array<std::array<Float, MAX_RUN_LENGTH>, MAX_RUN_LENGTH>, MAX_RUN_LENGTH> RUN_TO_RUN_PROB_CONV;

void initialize_channel(Float deletion_probability){
	for (size_t transmitted_run = 0; transmitted_run < MAX_RUN_LENGTH; ++transmitted_run)
	{
		for (size_t received_run = 0; received_run < MAX_RUN_LENGTH; ++received_run)
		{
			if (received_run > transmitted_run)
			{
				RUN_TO_RUN_PROB[transmitted_run][received_run] = 0.0;
				continue;
			}
			Float probability = pow(1 - deletion_probability, received_run) * pow(deletion_probability, transmitted_run - received_run);
			for (size_t i = 0; i < received_run; ++i)
			{
				probability *= (((Float) transmitted_run - i) / (i + 1));
			}
			RUN_TO_RUN_PROB[transmitted_run][received_run] = probability;
		}
	}

	for (size_t transmitted_run = 0; transmitted_run < MAX_RUN_LENGTH; ++transmitted_run)
	{
		for (size_t received_run = 0; received_run < MAX_RUN_LENGTH; ++received_run)
		{
			RUN_TO_RUN_PROB_TRANSPOSE[received_run][transmitted_run] = RUN_TO_RUN_PROB[transmitted_run][received_run];
		}
	}

	for (size_t transmitted_run = 0; transmitted_run < MAX_RUN_LENGTH; ++transmitted_run)
	{
		Float ndel_prob = (1 - len_to_len_transition_prob(transmitted_run, 0));
		std::array<Float, MAX_RUN_LENGTH> cached_table;
		for (size_t l1 = 0; l1 < MAX_RUN_LENGTH; ++l1)
		{
			cached_table[l1] = (len_to_len_transition_prob(transmitted_run, l1) / 
						ndel_prob);
		}
		for (size_t received_run = 0; received_run < MAX_RUN_LENGTH; ++received_run)
		{
			for (size_t tot_length = transmitted_run; tot_length < MAX_RUN_LENGTH; ++tot_length)
			{
				RUN_TO_RUN_PROB_CONV[transmitted_run][received_run][tot_length] = 0.0;
				for (size_t l1 = 1; l1 < received_run + 1; ++l1)
				{
					RUN_TO_RUN_PROB_CONV[transmitted_run][received_run][tot_length] += 
						cached_table[l1] * len_to_len_transition_prob(tot_length, received_run - l1);
				}
			}
		}
	}
}


std::vector<CodeWord> get_all_codewords(size_t r, size_t l){
	std::vector<CodeWord> res;
	if (r == 1)
	{
		for (size_t i = 1; i <= l; ++i)
		{
			res.push_back(CodeWord(std::vector<Run>({Run(0, i)})));
			res.push_back(CodeWord(std::vector<Run>({Run(1, i)})));
		}
		res.push_back(CodeWord());
		return res;
	}
	std::vector<CodeWord> r_minus_1 = get_all_codewords(r-1, l);
	for (auto iter = r_minus_1.begin(); iter != r_minus_1.end(); ++iter)
	{
		res.push_back(*iter);
		if (iter -> size() != r-1)
		{
			continue;
		}
		bool b = iter -> rbegin() -> value;
		for (size_t i = 1; i <= l - iter->total_length; ++i)
		{
			res.push_back((*iter) + Run(1 - b, i));
		}
	}
	return res;
}