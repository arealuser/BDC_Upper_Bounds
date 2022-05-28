#pragma once
#include "utils.h"
#include <array>
#include <vector>
#include <list>
#include <numeric>
#include <cassert>


const size_t MAX_RUN_LENGTH = 128;

// A 2D array mapping the length of a transmitted run and a received run to the probability that they will be happen in the channel.
extern std::array<std::array<Float, MAX_RUN_LENGTH>, MAX_RUN_LENGTH> RUN_TO_RUN_PROB;
// A transposed copy of this array for reducing cache misses.
extern std::array<std::array<Float, MAX_RUN_LENGTH>, MAX_RUN_LENGTH> RUN_TO_RUN_PROB_TRANSPOSE;

extern std::array<std::array<std::array<Float, MAX_RUN_LENGTH>, MAX_RUN_LENGTH>, MAX_RUN_LENGTH> RUN_TO_RUN_PROB_CONV;

inline Float len_to_len_transition_prob(const size_t& transmitted, const size_t& received, bool transposed=false){
	if (transposed)
	{
		return RUN_TO_RUN_PROB_TRANSPOSE[received][transmitted];
	} else{
		return RUN_TO_RUN_PROB[transmitted][received];
	}
}


/*
Compute the probability distribution of pairs of transmitted runs and received runs 
	(i.e. compute Pr(Bin(n, 1-d) = k) for all n,k <= MAX_RUN_LENGTH) and storest the results in RUN_TO_RUN_PROB
*/
void initialize_channel(Float deletion_probability);

struct Run{
	/*
	Our code will be made up of "runs" of 0s or 1s. Each run will have a positive length and a boolean value
	*/
	bool value;
	size_t length;
	inline Run(bool val=0, size_t len=0): value(val), length(len) {}
	inline Run(const Run& other): value(other.value), length(other.length) {}
	inline Run(Run&& other): value(other.value), length(other.length) {}
	inline bool operator== (const Run& other) const {return other.value == value and other.length == length;}
};

inline Float run_to_run_transition_prob(const Run& transmitted, const Run& received){
	if (transmitted.value != received.value)
	{
		return 0.0;
	}
	return len_to_len_transition_prob(transmitted.length, received.length);
}




// typedef std::vector<Run> CodeWord;

class CodeWord: public std::vector<Run>
{
public:
	size_t total_length;
	// template <typename ...Args>
	// inline CodeWord(Args... args): std::vector<Run>(args...), 
	// 							   total_length(accumulate(begin(), end(), 0, [](size_t s, const Run& r){return s + r.length;})) {}
	inline CodeWord(const std::vector<Run>& runs): std::vector<Run>(runs.begin(), runs.end()),
									total_length(accumulate(begin(), end(), 0, [](size_t s, const Run& r){return s + r.length;})) {}
	inline CodeWord(): std::vector<Run>(0), total_length(0) {}
	inline CodeWord(Run r): std::vector<Run>({r}), total_length(r.length) {}
	inline CodeWord(const CodeWord& other): std::vector<Run>(other.begin(), other.end()), total_length(other.total_length) {}
	inline CodeWord(CodeWord&& other): std::vector<Run>(std::move((std::vector<Run>) other)), total_length(other.total_length) {}
	// CodeWord& operator+= (const CodeWord& other)
	inline CodeWord& operator+= (const Run& r){
		push_back(r);
		total_length += r.length;
		return *this;
	}
	inline CodeWord operator+ (const Run& r) const{
		CodeWord other = CodeWord(*this);
		return other += r;
	}
};

/*
Returns the probability that the entire codeword was deleted.
*/
inline Float prob_all_were_deleted(const CodeWord& transmitted){
	Float res = 1.0;
	for (auto iter = transmitted.begin(); iter != transmitted.end(); ++iter)
	{
		res *= len_to_len_transition_prob(iter -> length, 0);
	}
	return res;
}

/*
Returns the probability that two codewords of the same number of runs to be the input/output of the channel.
This case is easy to compute, because the result is simply prod(P(t_i -> r_i))
*/
inline Float get_transition_prob_same_lengths(const CodeWord& transmitted, const CodeWord& received, bool transposed=false){
	// If both codewords are empty then the transition occurs w.p. 1.
	if (transmitted.size() == 0)
	{
		return 1.0;
	}

	// If the codewords don't match, then one cannot have produced the other.
	if (transmitted.begin() -> value != received.begin() -> value) 
	{
		return 0.0;
	}

	// If the codewords have the same number of runs, then they must also be alligned.
	Float res = 1.0;
	auto iter1 = transmitted.begin();
	auto iter2 = received.begin();
	for (; iter1 != transmitted.end(); ++iter1, ++iter2)
	{
		res *= len_to_len_transition_prob(iter1 -> length, iter2 -> length, transposed);
	}
	return res;
}

template <typename _InputIter>
inline void print_codeword(_InputIter begin, _InputIter end){
	printf("[");
	bool first = true;
	for(auto iter = begin; iter != end; ++iter){
		if (not first)
		{
			printf(", ");
		}
		printf("(%d ^ %lu)", iter -> value, iter -> length);
		first = false;
	}
	printf("]\n");
}

/*
Given possible transmitted and received codewords, computes the probability that the channel will transform one into the other.
*/
inline Float get_transition_prob(const CodeWord& transmitted, const CodeWord& received, bool verbose=false, bool transposed=false){

	if (verbose)
	{
		printf("Transmitted codeword =");
		print_codeword(transmitted.begin(), transmitted.end());
		printf("Received codeword =");
		print_codeword(received.begin(), received.end());
	}

	// If we received the same number of runs that we have transmitted, then we can just compare the runs one after
	// 	another, since they must be aligned to appear at all.
	if (transmitted.size() == received.size())
	{
		return get_transition_prob_same_lengths(transmitted, received, transposed);
	}
	// We cannot receive more runs than we transmit.
	if (transmitted.size() < received.size())
	{
		return 0.0;
	}
	if (received.size() == 0)
	{
		return prob_all_were_deleted(transmitted);
	}
	// Other special cases we can easily handle very efficiently:
	if (received.size() == transmitted.size() - 1)
	{
		if (received[0].value == transmitted[0].value)
		{
			Float product = len_to_len_transition_prob(transmitted[transmitted.size()-1].length, 0, true);
			for (size_t i = 0; i < received.size(); ++i)
			{
				product *= run_to_run_transition_prob(transmitted[i], received[i]);
			}
			return product;
		} else{
			Float product = len_to_len_transition_prob(transmitted[0].length, 0);
			for (size_t i = 0; i < received.size(); ++i)
			{
				product *= run_to_run_transition_prob(transmitted[i+1], received[i]);
			}
			return product;
		}
	}
	if ((received.size() == transmitted.size() - 2) and (received[0].value == transmitted[1].value))
	{
		Float product = len_to_len_transition_prob(transmitted[transmitted.size()-1].length, 0, true) * 
						len_to_len_transition_prob(transmitted[0].length, 0, true);
		if (verbose)
		{
			printf("%.3f%%\n", product);
		}
		for (size_t i = 0; i < received.size(); ++i)
		{
			product *= run_to_run_transition_prob(transmitted[i+1], received[i]);
		}
		if (verbose)
		{
			printf("%.3f%%\n", product);
		}
		return product;
	}

	
	// The number of runs that must have been deleted to get the number of received runs.
	size_t num_dels = transmitted.size() - received.size();

	// The probability that before we started recieving the ith run of the received codeword, j of the runs of the transmitted
	//  codeword were deleted.
	std::array<std::array<Float, 15>, 15> dynamic_programming_state;
	for (auto iter1 = dynamic_programming_state.begin(); iter1 != dynamic_programming_state.end(); ++iter1)
	{
		for (auto iter2 = iter1->begin(); iter2 != iter1->end(); ++iter2)
		{
			*iter2 = 0.0;
		}
	}

	dynamic_programming_state[0][0] = 1;
	for (size_t del_index = 1; del_index < num_dels + 1; ++del_index)
	{
		dynamic_programming_state[0][del_index] = dynamic_programming_state[0][del_index - 1] * 
				len_to_len_transition_prob(transmitted[del_index - 1].length, 0, true);
	}
	for (size_t del_index = 0; del_index < num_dels + 1; ++del_index)
	{
		dynamic_programming_state[0][del_index] *= (1 - len_to_len_transition_prob(transmitted[del_index].length, 0, true));
	}

	/*
	Use dynamic programming to compute dynamic_programming_state[i=rec_index][j=del_index]:
	*/
	for (size_t rec_index = 0; rec_index < received.size(); ++rec_index)
	{
		for (size_t del_index = 0; del_index <= num_dels; ++del_index)
		{
			// The relevant index in the transmitted codeword:
			size_t trans_index = rec_index + del_index;
			// If the received and transmitted codewords have different values, then they cannot be the correct alignment.
			if (received[rec_index].value != transmitted[trans_index].value)
			{
				continue;
			}

			size_t tot_length = 0;
			Float probability = dynamic_programming_state[rec_index][del_index];

			// The probability that the next block was not deleted:
			Float next_block_wasnt_deleted = 1.0;
			if (rec_index < received.size() - 1)
			{
				next_block_wasnt_deleted = 1 - len_to_len_transition_prob(transmitted[trans_index+1].length, 0, true);
			}

			if (verbose)
			{
				printf("Updating dynamic_programming_state[%lu][%lu] = %.3f%% + %.3f%% * %.3f%% * %.3f%% / %.3f\n", 
								rec_index+1, del_index, 100*dynamic_programming_state[rec_index+1][del_index],
								100*probability, 100*next_block_wasnt_deleted, 
								100*run_to_run_transition_prob(transmitted[trans_index], received[rec_index]),
								(1 - len_to_len_transition_prob(transmitted[trans_index].length, 0)));
				printf("%lu, %lu\n", rec_index, del_index);
			}
			// For del index not to change, we need for the next transmitted block not to be deleted:
			dynamic_programming_state[rec_index+1][del_index] += 
					probability * 
					next_block_wasnt_deleted * 
					run_to_run_transition_prob(transmitted[trans_index], received[rec_index]) / 
					(1 - len_to_len_transition_prob(transmitted[trans_index].length, 0, true));
			for (size_t next_del_index = del_index+1; next_del_index < num_dels + 1; ++next_del_index)
			{
				// The index in the transmitted message relevant to the next possible deletion index:
				size_t next_trans_index = rec_index + next_del_index;

				// If the next transmitted run doesn't have the same value as the current received one, then it must be deleted:
				if (received[rec_index].value != transmitted[next_trans_index].value)
				{
					probability *= len_to_len_transition_prob(transmitted[next_trans_index].length, 0);
				} else{
					// We increase the number of bits transmitted within this received block:
					tot_length += transmitted[next_trans_index].length;
					
					// The probability that the next block was not deleted:
					Float next_block_wasnt_deleted = 1.0;
					if (rec_index < received.size() - 1)
					{
						next_block_wasnt_deleted = 1 - len_to_len_transition_prob(transmitted[next_trans_index+1].length, 0);
					}

					dynamic_programming_state[rec_index+1][next_del_index] += 
							RUN_TO_RUN_PROB_CONV[transmitted[trans_index].length][received[rec_index].length][tot_length] *
							probability * next_block_wasnt_deleted;

					// We enumerate over the number of bits from the received run that came from the first block:
					// for (size_t l1 = 1; l1 < received[rec_index].length + 1; ++l1)
					// {
					// 	if (verbose)
					// 	{
					// 		printf("Updating dynamic_programming_state[%lu][%lu] = %.3f%% + %.3f%% * %.3f%% * %.3f%% * %.3f%%\n", 
					// 						rec_index+1, next_del_index, 100*dynamic_programming_state[rec_index+1][next_del_index],
					// 						100*(len_to_len_transition_prob(transmitted[trans_index].length, l1) / 
					// 							(1 - len_to_len_transition_prob(transmitted[trans_index].length, 0))), 
					// 						100*len_to_len_transition_prob(tot_length, received[rec_index].length - l1), 
					// 						100*probability, 100*next_block_wasnt_deleted);
					// 		printf("%lu, %lu\n", rec_index, del_index);
					// 	}

					// 	dynamic_programming_state[rec_index+1][next_del_index] +=
					// 		(len_to_len_transition_prob(transmitted[trans_index].length, l1) / 
					// 			(1 - len_to_len_transition_prob(transmitted[trans_index].length, 0))) *
					// 		len_to_len_transition_prob(tot_length, received[rec_index].length - l1) * probability *
					// 		next_block_wasnt_deleted;

					// 	if (verbose)
					// 	{
					// 		printf("to %.3f%%\n", 100*dynamic_programming_state[rec_index+1][next_del_index]);
					// 		printf("%lu, %lu\n", rec_index, del_index);
					// 	}
					// }
				}
			}
		}
	}

	if (verbose)
	{
		printf("dynamic_programming_state = [\n");
		for (size_t i = 0; i < received.size() + 1; ++i)
		{
			printf("\t[");
			for (size_t j = 0; j < num_dels + 1; ++j)
			{
				printf("%.3f%%\t", 100*dynamic_programming_state[i][j]);
			}
			printf("]\n");
		}
		printf("]\n");
	}

	Float result = 0.0;
	for (size_t del_index = num_dels-1; del_index < num_dels+1; ++del_index)
	{
		Float base_probability = dynamic_programming_state[received.size()][del_index];
		Float deletion_probability = 1.0;
		size_t trans_index = del_index + received.size();
		for (size_t trans_index2 = trans_index; trans_index2 < transmitted.size(); ++trans_index2)
		{
			deletion_probability *= len_to_len_transition_prob(transmitted[trans_index2].length, 0);
		}
		if (verbose)
		{
			printf("The %luth contribution to res is %.3f%% * %.3f%% = %.3f%%\n", del_index+1, base_probability*100,
				deletion_probability*100, base_probability*deletion_probability*100);
		}
			
		result += base_probability * deletion_probability;
	}
	if (verbose)
	{
		printf("result = %.3f%%\n", 100*result);
	}
	return result;
}



/*
Returns a list of all the codewords with at most r runs and a total length of at most l
*/
std::vector<CodeWord> get_all_codewords(size_t r, size_t l);



