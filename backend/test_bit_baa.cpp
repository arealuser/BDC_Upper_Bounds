#include "channel.h"
#include "bit_channel.h"
#include "parallelized_baa.h"
#include "bit_baa.h"
#include <algorithm>
#include <ctime>
#include <cassert>
#include <cmath>
#include <boost/range/combine.hpp>


int main()
{
	Float deletion_probability = 0.5;
	initialize_channel(deletion_probability);
	constexpr size_t in_len = 10;
	constexpr size_t out_len = 5;
	initialize_bit_channel(deletion_probability, in_len, out_len, false);

	auto transmitted_codewords = get_all_bit_codewords(in_len);
	auto received_codewords = get_all_bit_codewords(out_len);

	printf("%lu (= 2^%.1f) possible transmitted codewords\n", transmitted_codewords.size(), log(transmitted_codewords.size()) / log(2));
	printf("%lu (= 2^%.1f) possible received codewords\n", received_codewords.size(), log(received_codewords.size()) / log(2));
	printf("In total P_jk has 2^%.1f entries\n", (log(transmitted_codewords.size()) + log(received_codewords.size())) / log(2));
	auto t0 = clock();

	std::vector<Float> Q;
	Q.resize(transmitted_codewords.size());
	std::for_each(Q.begin(), Q.end(), [transmitted_codewords](Float& Q){Q = 1.0 / transmitted_codewords.size();});


	for (int i = 0; i < 101; ++i)
	{
		if (i % 20 == 0)
		{
			printf("Running the %dth step of the BAA algorithm (%.2f seconds)...\n", i+1, ((float) (clock() - t0)) / CLOCKS_PER_SEC);
			printf("total probs = %.2f%%\t", 100 * (std::accumulate(Q.begin(), Q.end(), 0.0)));
			printf("min prob = %.2f%%\t", 100 * (*std::min_element(Q.begin(), Q.end())));
			printf("max prob = %.2f%%\n", 100 * (*std::max_element(Q.begin(), Q.end())));

			auto rate = compute_rate(transmitted_codewords, received_codewords, Q);
			printf("The current rate is %f\n", rate);
			auto log_dens = compute_all_log_Wjk_den(transmitted_codewords, received_codewords, Q);
			Float rate2 = compute_bit_rate_efficient(transmitted_codewords, received_codewords, log_dens, Q);
			printf("Efficient rate computation: %f\n", rate2);
			assert(std::abs(rate-rate2) < 1E-6);
		}
		Q = do_full_baa_step(transmitted_codewords, received_codewords, Q);
	}

		

	// std::vector<Float> Q2 = Qs[0];
	// Float s = 0;
	// for(auto q_12 : boost::combine(Q, Q2)){
	// 	Float q1, q2;
	// 	boost::tie(q1, q2) = q_12;
	// 	s += std::abs(q1 - q2);
	// }

	// printf("TVD between implementations: %f%%\n", s * 100);

	// s = 0;
	// for(size_t i; i < transmitted_codewords.size(); ++i){
	// 	if (Q[i] >= Qs[0][i])
	// 	{
	// 		s += Q[i] - Qs[0][i];
	// 	}else{
	// 		s -= Q[i] - Qs[0][i];
	// 	}
	// }

	// printf("TVD between implementations: %f%%\n", s * 100);
	return 0;
}