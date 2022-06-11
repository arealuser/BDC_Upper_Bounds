#include "bit_baa.h"
#include <algorithm>
#include <boost/range/combine.hpp>
#include <cmath>


struct Sum
{
    void operator()(Float n) { sum += n; }
    Float sum{0};
};

std::vector<Float> do_full_baa_step(const std::vector<BitCodeWord>& transmitted, const std::vector<BitCodeWord>& received, 
	const std::vector<Float>& Q_i){
	std::vector<Float> log_W_jk = compute_all_log_Wjk_den (transmitted, received, Q_i);
	std::vector<Float> log_alphas = compute_all_log_alpha_k (transmitted, received, Q_i, log_W_jk);

	// Align alphas so that they are not all small and none of them are too huge for more accurate numerics:
	Float max_log_alpha = *std::max_element(log_alphas.begin(), log_alphas.end());
	std::for_each(log_alphas.begin(), log_alphas.end(), [max_log_alpha](Float &log_alpha){ log_alpha = exp(log_alpha - max_log_alpha);});
	std::vector<Float> alphas = std::move(log_alphas);

	// Normalize the alphas by their sum:
	Float alpha_sum = std::accumulate(alphas.begin(), alphas.end(), 0.0);
	std::for_each(alphas.begin(), alphas.end(), [alpha_sum](Float &alpha){ alpha /= alpha_sum;});

	return alphas;
}

std::vector<Float> compute_all_log_Wjk_den (const std::vector<BitCodeWord>& transmitted, const std::vector<BitCodeWord>& received, 
	const std::vector<Float>& Q_i){
	// Iteratively call compute_log_Wjk_den for each possible received codeword.
	std::vector<Float> log_Wjk_den;
	log_Wjk_den.reserve(received.size());
	for(auto rec_iter = received.begin(); rec_iter != received.end(); ++rec_iter){
		log_Wjk_den.push_back(compute_log_Wjk_den(transmitted, *rec_iter, Q_i));
	}
	return log_Wjk_den;
}

std::vector<Float> compute_all_log_alpha_k (const std::vector<BitCodeWord>& transmitted, const std::vector<BitCodeWord>& received, 
	const std::vector<Float>& Q_i, const std::vector<Float>& log_W_jk_den){
	// Iteratively call compute_log_alpha_k for each possible transmitted codeword.
	std::vector<Float> log_alphas;
	log_alphas.reserve(transmitted.size());
	auto Q_iter = Q_i.begin();
	for(auto trans_iter = transmitted.begin(); trans_iter != transmitted.end(); ++trans_iter){
		log_alphas.push_back(compute_log_alpha_k(*trans_iter, received, *(Q_iter++), log_W_jk_den));
	}
	return log_alphas;
}


Float compute_log_Wjk_den (const std::vector<BitCodeWord>& transmitted, const BitCodeWord& received, const std::vector<Float>& Q_i){
	std::vector<Float> probs_col = compute_Pjk_col(transmitted, received);
	Float denominator = 0.0;
	for (auto pr_Qi : boost::combine(probs_col, Q_i))
	{
		Float P_ji, Q_i;
		boost::tie(P_ji, Q_i) = pr_Qi;
		denominator += P_ji * Q_i;
	}
	return log(denominator);
}


Float compute_log_alpha_k (const BitCodeWord& transmitted, const std::vector<BitCodeWord>& received, 
	Float Q_k, const std::vector<Float>& log_W_jk_den){
	Float log_Q_k = log(Q_k);
	std::vector<Float> probs_row = compute_Pjk_row(transmitted, received);
	Float log_alpha = 0.0;
	for(auto pr_den : boost::combine(probs_row, log_W_jk_den)){
		Float P_jk, log_den;
		boost::tie(P_jk, log_den) = pr_den;
		// Reduces the runtime of the algorithm by 30% with minimal change to the outcome.
		if (P_jk < 1E-12)
		{
			continue;
		}
		log_alpha += P_jk * (log_Q_k + log(P_jk) - log_den);
		// log_alpha += P_jk;
	}
	return log_alpha;
}



std::vector<Float> compute_Pjk_row(const BitCodeWord& transmitted, const std::vector<BitCodeWord>& received){
	std::vector<Float> res; res.reserve(received.size());
	for(auto rec_iter = received.begin(); rec_iter != received.end(); ++rec_iter){
		res.push_back(get_bit_transition_prob(transmitted, *rec_iter));
	}
	return res;
}

std::vector<Float> compute_Pjk_col(const std::vector<BitCodeWord>& transmitted, const BitCodeWord& received){
	std::vector<Float> res; res.reserve(transmitted.size());
	for(auto trans_iter = transmitted.begin(); trans_iter != transmitted.end(); ++trans_iter){
		res.push_back(get_bit_transition_prob(*trans_iter, received));
	}
	return res;
}


Float compute_bit_rate_efficient(const std::vector<BitCodeWord>& transmitted, const std::vector<BitCodeWord>& received, 
	const std::vector<Float>& log_W_jk_den, const std::vector<Float>& Q_i){
	size_t n_I = transmitted.size();
	size_t n_J = received.size();
	assert(transmitted.size() == Q_i.size());
	Float rate = 0.0;

	for(size_t i = 0; i < n_I; ++i){
		auto probs_row = compute_Pjk_row(transmitted[i], received);
		std::vector<Float> log_probs_row = probs_row;
		for_each(log_probs_row.begin(), log_probs_row.end(), [](Float& x){x = log(x);});
		Float Qk = Q_i[i];
		for(size_t j = 0; j < n_J; ++j){
			Float log_den = log_W_jk_den[j];
			Float P_jk = probs_row[j];
			Float log_P_jk = log_probs_row[j];
			if (P_jk < 1E-20)
			{
				continue;
			}
			// printf("%lu\t%lu\t%.2f%%\t%.2f%%\t%.2f\t%.2f\n", i, j, 100*Qk, 100*P_jk, log_P_jk, log_den);
			rate += Qk * P_jk * (log_P_jk - log_den);
		}
	}
	return rate;
}


Float compute_rate(const std::vector<BitCodeWord>& transmitted, const std::vector<BitCodeWord>& received, 
	const std::vector<Float>& Q_i){
	std::vector<std::vector<Float> > prob_table;

	size_t n_I = transmitted.size();	prob_table.resize(n_I);
	size_t n_J = received.size();		for_each(prob_table.begin(), prob_table.end(), [n_J](std::vector<Float>& row){
											row.resize(n_J);
	});

	for (size_t i = 0; i < n_I; ++i)
	{
		for (size_t j = 0; j < n_J; ++j)
		{
			prob_table[i][j] = get_bit_transition_prob(transmitted[i], received[j]);
		}
	}

	std::vector<Float> denominator; denominator.resize(n_J);
	for (size_t j = 0; j < n_J; ++j)
	{
		denominator[j] = 0.0;
		for (size_t i = 0; i < n_I; ++i)
		{
			denominator[j] += prob_table[i][j] * Q_i[i];
		}
	}

	Float rate = 0.0;
	for (size_t k = 0; k < n_I; ++k)
	{
		for (size_t j = 0; j < n_J; ++j)
		{
			if (prob_table[k][j] < 1E-30)
			{
				continue;
			}
			if (std::isnan(denominator[j]) or denominator[j] < 1E-50)
			{
				denominator[j] = 1E-50;
			}
			// printf("%lu\t%lu\t%.2f%%\t%.2f%%\t%.2f\t%.2f\n", k, j, 100*Q_i[k], 100*prob_table[k][j], log(prob_table[k][j]), log(denominator[j]));
			rate += Q_i[k] * prob_table[k][j] * log(prob_table[k][j] / denominator[j]);
		}
	}

	return rate;
}