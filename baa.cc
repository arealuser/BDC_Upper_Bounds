#include "baa.h"
#include <algorithm>
#include <boost/range/combine.hpp>
#include <cmath>


struct Sum
{
    void operator()(Float n) { sum += n; }
    Float sum{0};
};

std::vector<Float> do_full_baa_step(const std::vector<CodeWord>& transmitted, const std::vector<CodeWord>& received, 
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

std::vector<Float> compute_all_log_Wjk_den (const std::vector<CodeWord>& transmitted, const std::vector<CodeWord>& received, 
	const std::vector<Float>& Q_i){
	// Iteratively call compute_log_Wjk_den for each possible received codeword.
	std::vector<Float> log_Wjk_den;
	log_Wjk_den.reserve(received.size());
	for(auto rec_iter = received.begin(); rec_iter != received.end(); ++rec_iter){
		log_Wjk_den.push_back(compute_log_Wjk_den(transmitted, *rec_iter, Q_i));
	}
	return log_Wjk_den;
}

std::vector<Float> compute_all_log_alpha_k (const std::vector<CodeWord>& transmitted, const std::vector<CodeWord>& received, 
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


Float compute_log_Wjk_den (const std::vector<CodeWord>& transmitted, const CodeWord& received, const std::vector<Float>& Q_i){
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


Float compute_log_alpha_k (const CodeWord& transmitted, const std::vector<CodeWord>& received, 
	Float Q_k, const std::vector<Float>& log_W_jk_den){
	Float log_Q_k = log(Q_k);
	std::vector<Float> probs_row = compute_Pjk_row(transmitted, received);
	Float log_alpha = 0.0;
	for(auto pr_den : boost::combine(probs_row, log_W_jk_den)){
		Float P_jk, log_den;
		boost::tie(P_jk, log_den) = pr_den;
		log_alpha += P_jk * (log_Q_k + log(P_jk) - log_den);
	}
	return log_alpha;
}



std::vector<Float> compute_Pjk_row(const CodeWord& transmitted, const std::vector<CodeWord>& received){
	std::vector<Float> res; res.reserve(received.size());
	for(auto rec_iter = received.begin(); rec_iter != received.end(); ++rec_iter){
		res.push_back(get_transition_prob(transmitted, *rec_iter, false));
	}
	return res;
}

std::vector<Float> compute_Pjk_col(const std::vector<CodeWord>& transmitted, const CodeWord& received){
	std::vector<Float> res; res.reserve(transmitted.size());
	for(auto trans_iter = transmitted.begin(); trans_iter != transmitted.end(); ++trans_iter){
		res.push_back(get_transition_prob(*trans_iter, received, false, true));
	}
	return res;
}