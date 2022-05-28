#include "channel.h"
#include "bit_channel.h"
#include "parallelized_baa.h"
#include "baa.h"
#include <algorithm>
#include <ctime>
#include <cassert>
#include <cmath>
#include <boost/range/combine.hpp>


int main(int argc, char const *argv[])
{
	Float deletion_probability = 0.5;
	initialize_channel(deletion_probability);

	auto transmitted_codewords = get_all_codewords(10, 10);
	auto received_codewords = get_all_codewords(10, 10);

	printf("%lu (= 2^%.1f) possible transmitted codewords\n", transmitted_codewords.size(), log(transmitted_codewords.size()) / log(2));
	printf("%lu (= 2^%.1f) possible received codewords\n", received_codewords.size(), log(received_codewords.size()) / log(2));
	printf("In total P_jk has 2^%.1f entries\n", (log(transmitted_codewords.size()) + log(received_codewords.size())) / log(2));
	size_t total_length = 0;
	auto t0 = clock();

	std::vector<Float> Q;
	Q.resize(transmitted_codewords.size());
	std::for_each(Q.begin(), Q.end(), [transmitted_codewords](Float& Q){Q = 1.0 / transmitted_codewords.size();});

	std::vector<std::vector<Float> > Qs = {Q};


	for (int i = 0; i < 2; ++i)
	{
		printf("Running the %dth step of the BAA algorithm (%.2f seconds)...\n", i+1, ((float) (clock() - t0)) / CLOCKS_PER_SEC);
		printf("total probs = %.2f%%\t", 100 * (std::accumulate(Q.begin(), Q.end(), 0.0)));
		printf("min prob = %.2f%%\t", 100 * (*std::min_element(Q.begin(), Q.end())));
		printf("max prob = %.2f%%\n", 100 * (*std::max_element(Q.begin(), Q.end())));

		printf("total probs = %.2f%%\t", 100 * (std::accumulate(Qs[0].begin(), Qs[0].end(), 0.0)));
		printf("min prob = %.2f%%\t", 100 * (*std::min_element(Qs[0].begin(), Qs[0].end())));
		printf("max prob = %.2f%%\n", 100 * (*std::max_element(Qs[0].begin(), Qs[0].end())));

		auto rates = compute_rate_parallelized(transmitted_codewords, received_codewords, Qs);
		printf("The current rate is %f\n", *rates.begin());

		auto rate = compute_rate_naive(transmitted_codewords, received_codewords, Q);
		printf("The current rate is %f\n", rate);

		Q = do_baa_step_naive(transmitted_codewords, received_codewords, Q);
		// Q = do_full_baa_step(transmitted_codewords, received_codewords, Q);
		Qs = do_full_baa_step_parallelized(transmitted_codewords, received_codewords, Qs);
	}

		

	std::vector<Float> Q2 = Qs[0];
	Float s = 0;
	for(auto q_12 : boost::combine(Q, Q2)){
		Float q1, q2;
		boost::tie(q1, q2) = q_12;
		s += abs(q1 - q2);
	}

	printf("%f\n", s);

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