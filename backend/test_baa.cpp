#include "channel.h"
#include "bit_channel.h"
#include "parallelized_baa.h"
#include "baa.h"
#include <algorithm>
#include <ctime>
#include <cassert>
#include <cmath>


int main()
{
	Float deletion_probability = 0.5;
	initialize_channel(deletion_probability);

	auto transmitted_codewords = get_all_codewords(10, 10);
	auto received_codewords = get_all_codewords(10, 10);

	printf("%lu (= 2^%.1f) possible transmitted codewords\n", transmitted_codewords.size(), log(transmitted_codewords.size()) / log(2));
	printf("%lu (= 2^%.1f) possible received codewords\n", received_codewords.size(), log(received_codewords.size()) / log(2));
	printf("In total P_jk has 2^%.1f entries\n", (log(transmitted_codewords.size()) + log(received_codewords.size())) / log(2));
	auto t0 = clock();

	std::vector<Float> Q;
	Q.resize(transmitted_codewords.size());
	std::for_each(Q.begin(), Q.end(), [transmitted_codewords](Float& Q){Q = 1.0 / transmitted_codewords.size();});

	std::vector<std::vector<Float> > Qs; Qs.resize(1); Qs[0].resize(transmitted_codewords.size());
	std::for_each(Qs[0].begin(), Qs[0].end(), [transmitted_codewords](Float& Q){Q = 1.0 / transmitted_codewords.size();});


	for (int i = 0; i < 3; ++i)
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

		rate = compute_rate_naive(transmitted_codewords, received_codewords, Qs[0]);
		printf("The current rate is %f\n", rate);

		Q = do_baa_step_naive(transmitted_codewords, received_codewords, Q);
		// Q = do_full_baa_step(transmitted_codewords, received_codewords, Q);
		Qs = do_full_baa_step_parallelized(transmitted_codewords, received_codewords, Qs);
	}

		

	std::vector<Float> Q2 = Qs[0];
	Float s = 0;
	for(size_t i = 0; i < Q.size(); ++i){
		s += std::abs(Q[i] - Q2[i]);
	}

	printf("TVD between implementations: %f%%\n", s * 100);
	assert(s < 1E-6);
	return 0;
}