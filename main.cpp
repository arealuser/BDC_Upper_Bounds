#include "channel.h"
#include "bit_channel.h"
#include "baa.h"
#include <algorithm>
#include <ctime>
#include <cassert>
#include <cmath>

int main(int argc, char const *argv[])
{

	Float deletion_probability = 0.9;
	initialize_channel(deletion_probability);

	auto transmitted_codewords = get_all_codewords(3, 50);
	auto received_codewords = get_all_codewords(3, 20);

	printf("%lu (= 2^%.1f) possible transmitted codewords\n", transmitted_codewords.size(), log(transmitted_codewords.size()) / log(2));
	printf("%lu (= 2^%.1f) possible received codewords\n", received_codewords.size(), log(received_codewords.size()) / log(2));
	printf("In total P_jk has 2^%.1f entries\n", (log(transmitted_codewords.size()) + log(received_codewords.size())) / log(2));
	Float s = 0;
	size_t total_length = 0;
	auto t0 = clock();

	std::vector<Float> Q;
	Q.resize(transmitted_codewords.size());
	std::for_each(Q.begin(), Q.end(), [transmitted_codewords](Float& Q){Q = 1.0 / transmitted_codewords.size();});

	for (int i = 0; i < 5; ++i)
	{
		printf("Running the %dth step of the BAA algorithm (%.2f seconds)...\n", i+1, ((float) (clock() - t0)) / CLOCKS_PER_SEC);
		Q = do_full_baa_step(transmitted_codewords, received_codewords, Q);
	}
	return 0;

	auto t1 = clock();
	size_t counter = 0;
	for(auto trans_iter = transmitted_codewords.begin(); trans_iter != transmitted_codewords.end(); ++trans_iter){
		std::vector<Float> probs(received_codewords.size()), probs2(received_codewords.size()), cumsum(received_codewords.size());
		auto rec_iter = received_codewords.begin();
		std::generate(probs.begin(), probs.end(), 
			[&](){
				return get_transition_prob(*trans_iter, *(rec_iter++));
		});
		
		s += *probs.rbegin() + *probs2.rbegin();
		if (++counter % (1 + (transmitted_codewords.size()) / 100) == 0)
		{
			printf("TQDM: %lu / %lu iterations (%.1f%%) in %.1f seconds\r", counter, transmitted_codewords.size(), 
																			(100.0 * counter) / transmitted_codewords.size(),
																			((float) (clock() - t1)) / CLOCKS_PER_SEC); fflush(stdout);
		}
	}
	printf("\n");
	printf("%.1f\n", s);
	printf("This took %.1f seconds.\n", ((float) (clock() - t1)) / CLOCKS_PER_SEC); fflush(stdout);
	printf("This means we have 2 ** %.1f iterations per second\n", 
		log((((Float) transmitted_codewords.size()) * (received_codewords.size())) / 
			(((Float) (clock() - t1)) / CLOCKS_PER_SEC)) / log(2));
	return 0;
}