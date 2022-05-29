#include "channel.h"
#include "bit_channel.h"
#include <algorithm>
#include <ctime>
#include <cassert>
#include <cmath>

int main(int argc, char const *argv[])
{
	// CodeWord transmitted = {Run(0, 2), Run(1, 1), Run(0, 2), Run(1, 1), Run(0, 2), Run(1, 1), Run(0, 2)};
	// CodeWord received = {Run(0,4)};

	// if (argc != 2)
	// {
	// 	printf("Usage: %s output_file\n", argv[0]);
	// 	exit(1);
	// }

	// FILE* fout = fopen(argv[1], "w");

	CodeWord transmitted = std::vector<Run>({Run(0, 9), Run(1, 1), Run(0, 1)});
	CodeWord received = std::vector<Run>({Run(0, 5)});

	Float deletion_probability = 0.9;

	initialize_channel(deletion_probability);
	printf("The probablity of going from the transmitted to the received codeword is %.4f%%.\n", 
		100*get_transition_prob(transmitted, received));

	printf("The bit_channel probablity of going from the transmitted to the received codeword is %.4f%%.\n", 
		100*get_bit_transition_prob(convert_to_bit_word(transmitted), convert_to_bit_word(received), deletion_probability));
	// return 0;


	auto all_codewords = get_all_codewords(2, 4);
	printf("Got a total of %lu codewords.\n", all_codewords.size());

	auto transmitted_codewords = get_all_codewords(3, 40);
	std::vector<BitCodeWord> bit_transmitted_codewords(transmitted_codewords.size());
	auto t_iter = transmitted_codewords.begin();
	std::generate(bit_transmitted_codewords.begin(), bit_transmitted_codewords.end(), [&](){return convert_to_bit_word(*(t_iter++));});
	for (size_t i = 0; i < transmitted_codewords.size(); ++i)
	{
		assert(convert_to_run_word(bit_transmitted_codewords[i]) == transmitted_codewords[i]);
	}

	auto received_codewords = get_all_codewords(3, 20);
	std::vector<BitCodeWord> bit_received_codewords(received_codewords.size());
	auto r_iter = received_codewords.begin();
	std::generate(bit_received_codewords.begin(), bit_received_codewords.end(), [&](){return convert_to_bit_word(*(r_iter++));});
	for (size_t i = 0; i < received_codewords.size(); ++i)
	{
		assert(convert_to_run_word(bit_received_codewords[i]) == received_codewords[i]);
	}

	printf("%lu possible transmitted codewords\n", transmitted_codewords.size());
	printf("%lu possible received codewords\n", received_codewords.size());
	Float s = 0;
	// size_t total_length = 0;
	auto t0 = clock();
	size_t counter = 0;

	auto bit_trans_iter = bit_transmitted_codewords.begin();
	for(auto trans_iter = transmitted_codewords.begin(); trans_iter != transmitted_codewords.end(); ++trans_iter){
		std::vector<Float> probs(received_codewords.size()), probs2(received_codewords.size()), cumsum(received_codewords.size());
		auto rec_iter = received_codewords.begin();
		std::generate(probs.begin(), probs.end(), 
			[&](){
				return get_transition_prob(*trans_iter, *(rec_iter++));
		});

		auto brec_iter = bit_received_codewords.begin();
		std::generate(probs2.begin(), probs2.end(), 
			[&](){
				return get_bit_transition_prob(*bit_trans_iter, *(brec_iter++), deletion_probability);
		});


		for (size_t ir = 0; ir < received_codewords.size(); ++ir)
		{
			if (std::abs(probs[ir] - probs2[ir]) > 1E-6)
			{
				printf("Transmitted Codeword:\n");
				print_codeword(trans_iter->begin(), trans_iter->end());
				printf("Received Codeword[%lu]:\n", ir);
				print_codeword(received_codewords[ir].begin(), received_codewords[ir].end());
				printf("prob1 = %f, prob2 = %f, diff = %f\n", probs[ir], probs2[ir], probs[ir] - probs2[ir]);
				get_bit_transition_prob(*bit_trans_iter, bit_received_codewords[ir], deletion_probability, true);
				assert(false);
			}
		}

		// assert(probs == probs2);

		// auto brec_iter = bit_received_codewords.begin();
		// 
		// std::sort(probs.begin(), probs.end());
		// std::partial_sum(probs.rbegin(), probs.rend(), cumsum.begin());
		// if(not ((*cumsum.rbegin() >= 1.0-1E-6) and (*cumsum.rbegin() <= 1.0+1E-6))){
		// 	print_codeword(trans_iter->begin(), trans_iter->end());
		// 	get_transition_prob(*trans_iter, CodeWord(std::vector<Run>({Run(0,1)})), true);
		// 	get_transition_prob(*trans_iter, CodeWord(std::vector<Run>({Run(1,1)})), true);
		// 	print_codeword(trans_iter->begin(), trans_iter->end());
		// 	printf("%f\n", *cumsum.rbegin() - 1.0);
		// 	assert((*cumsum.rbegin() >= 1.0-1E-6) and (*cumsum.rbegin() <= 1.0+1E-6));
		// }

		// for (size_t num_important = 0; num_important < cumsum.size(); ++num_important)
		// {
		// 	if (cumsum[num_important] > 1.0 - 1E-6)
		// 	{
		// 		fprintf(fout, "%lu\n", num_important);
		// 		break;
		// 	}
		// }

		
		s += *probs.rbegin() + *probs2.rbegin();
		if (++counter % (1 + (transmitted_codewords.size()) / 100) == 0)
		{
			printf("TQDM: %lu / %lu iterations (%.1f%%) in %.1f seconds\r", counter, transmitted_codewords.size(), 
																			(100.0 * counter) / transmitted_codewords.size(),
																			((float) (clock() - t0)) / CLOCKS_PER_SEC); fflush(stdout);
		}
		++bit_trans_iter;
	}
	printf("\n");
	printf("%.1f\n", s);
	printf("This took %.1f seconds.\n", ((float) (clock() - t0)) / CLOCKS_PER_SEC); fflush(stdout);
	printf("This means we have 2 ** %.1f iterations per second\n", 
		log((((Float) transmitted_codewords.size()) * (received_codewords.size())) / 
			(((Float) (clock() - t0)) / CLOCKS_PER_SEC)) / log(2));
	return 0;
}