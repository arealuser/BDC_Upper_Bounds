#include "channel.h"
#include "bit_channel.h"
#include "bit_baa_fast.h"
#include "parallelized_baa.h"
#include "bit_baa.h"
#include <algorithm>
#include <ctime>
#include <cassert>
#include <cmath>


int main()
{
	Float deletion_probability = 0.5;
	initialize_channel(deletion_probability);
	constexpr size_t in_len = 15;
	constexpr size_t out_len = 5;
	initialize_bit_channel(deletion_probability, in_len, out_len, false);

	auto transmitted_codewords = get_all_bit_codewords(in_len);
	std::sort(transmitted_codewords.begin(), transmitted_codewords.end(), 
		[](const EfficientBitCodeWord& a, const EfficientBitCodeWord& b) {return a < b;});
	auto received_codewords = get_all_bit_codewords(out_len);
	std::sort(received_codewords.begin(), received_codewords.end(), 
		[](const EfficientBitCodeWord& a, const EfficientBitCodeWord& b) {return a < b;});

	std::vector<EfficientBitCodeWord> transmitted_codewords_efficient(transmitted_codewords.begin(), transmitted_codewords.end());
	std::sort(transmitted_codewords_efficient.begin(), transmitted_codewords_efficient.end());
	transmitted_codewords_efficient = get_transmitted_codewords_symmetries(transmitted_codewords_efficient);
	std::vector<EfficientBitCodeWord> received_codewords_efficient(received_codewords.begin(), received_codewords.end());
	std::sort(received_codewords_efficient.begin(), received_codewords_efficient.end());

	printf("%lu (= 2^%.1f) possible transmitted codewords\n", transmitted_codewords.size(), log(transmitted_codewords.size()) / log(2));
	printf("%lu (= 2^%.1f) possible received codewords\n", received_codewords.size(), log(received_codewords.size()) / log(2));
	printf("In total P_jk has 2^%.1f entries\n", (log(transmitted_codewords.size()) + log(received_codewords.size())) / log(2));
	auto t0 = clock();

	std::vector<Float> Q;
	Q.resize(transmitted_codewords_efficient.size());
	std::for_each(Q.begin(), Q.end(), [transmitted_codewords_efficient](Float& Q){Q = 1.0 / transmitted_codewords_efficient.size();});

	std::vector<Float> Q2;
	Q2.resize(transmitted_codewords.size());
	std::for_each(Q2.begin(), Q2.end(), [transmitted_codewords](Float& Q2){Q2 = 1.0 / transmitted_codewords.size();});


	for (int i = 0; i < 151; ++i)
	{
		// Float tvd1 = 0.0;
		// Float tvd2 = 0.0;
		// for (size_t j = 0; j < Q.size(); ++j)
		// {
		// 	tvd1 += std::abs(Q[j] - (2 * Q2[j*2]));
		// 	tvd2 += std::abs(Q2[2*j + 1] - Q2[j*2]);
		// }
		// printf("%f, %f\n", tvd1, tvd2);
		if ((i % 30) == 0)
		{
			printf("Running the %dth step of the BAA algorithm (%.2f seconds)...\n", i+1, ((float) (clock() - t0)) / CLOCKS_PER_SEC);
			printf("total probs = %.2f%%\t", 100 * (std::accumulate(Q.begin(), Q.end(), 0.0)));
			printf("min prob = %.2f%%\t", 100 * (*std::min_element(Q.begin(), Q.end())));
			printf("max prob = %.2f%%\n", 100 * (*std::max_element(Q.begin(), Q.end())));

			// auto rate = compute_rate(transmitted_codewords_efficient, received_codewords_efficient, Q);
			// printf("The current rate is %f\n", rate);
			auto log_dens = compute_all_log_Wjk_den(transmitted_codewords_efficient, received_codewords_efficient, Q);
			std::vector<Float> dens; dens.resize(log_dens.size());
			for (size_t j = 0; j < dens.size(); ++j)
			{
				dens[j] = exp(log_dens[j]);
			}
			printf("sum(dens)=%f\n", std::accumulate(dens.begin(), dens.end(), 0.0));
			printf("sum(Q)=%f\n", std::accumulate(Q.begin(), Q.end(), 0.0));
			Float rate2 = compute_bit_rate_efficient(transmitted_codewords_efficient, received_codewords_efficient, log_dens, Q);
			printf("Efficient rate computation: %f\n", rate2 / log(2));
			// assert(std::abs(rate-rate2) < 1E-6);
		}
		Q = do_full_baa_step(transmitted_codewords_efficient, received_codewords_efficient, Q);
		// Q2 = do_full_baa_step(transmitted_codewords, received_codewords, Q2);
	}
	return 0;
}