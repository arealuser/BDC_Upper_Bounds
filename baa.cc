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


std::vector<Float> do_baa_step_naive(const std::vector<CodeWord>& transmitted, const std::vector<CodeWord>& received, 
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
			prob_table[i][j] = get_transition_prob(transmitted[i], received[j]);
		}
		if (not (abs(std::accumulate(prob_table[i].begin(), prob_table[i].end(), 0.0) - 1) < 1E-10))
		{
			printf("%f\n", std::accumulate(prob_table[i].begin(), prob_table[i].end(), 0.0));
			printf("%f\n", prob_table[i][i]);
			printf("%lu\n", i);
			auto print = [](const Run& r){printf("(%lu, %d)\t", r.length, r.value);};
			printf("Annoying codeword:\n");
			for_each(transmitted[i].begin(), transmitted[i].end(), print);
			printf("\n");
			for (size_t j = 0; j < n_J; ++j)
			{
				if (prob_table[i][j] > 0.1)
				{
					printf("Accomplice (%f, %lu):\n", prob_table[i][j], j);
					get_transition_prob(transmitted[i], received[j], true);
					for_each(transmitted[j].begin(), transmitted[j].end(), print);
					printf("\n");
				}
			}
		}
		assert(abs(std::accumulate(prob_table[i].begin(), prob_table[i].end(), 0.0) - 1) < 1E-10);
		assert(abs(prob_table[i][i] - 1) < 1E-10);
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

	std::vector<Float> log_alphas; log_alphas.resize(n_I);
	for (size_t k = 0; k < n_I; ++k)
	{
		log_alphas[k] = 0.0;
		for (size_t j = 0; j < n_J; ++j)
		{
			if (prob_table[k][j] < 1E-12)
			{
				continue;
			}
			if (std::isnan(denominator[j]) or denominator[j] < 1E-50)
			{
				denominator[j] = 1E-50;
			}
			log_alphas[k] += prob_table[k][j] * log(Q_i[k] * prob_table[k][j] / denominator[j]);
		}
	}

	Float max_log_alpha = *max_element(log_alphas.begin(), log_alphas.end());
	for_each(log_alphas.begin(), log_alphas.end(), [max_log_alpha](Float& log_alpha){
		log_alpha = exp(log_alpha - max_log_alpha);
	});

	Float sum = std::accumulate(log_alphas.begin(), log_alphas.end(), 0.0);
	for_each(log_alphas.begin(), log_alphas.end(), [sum](Float& alpha){
		alpha /= sum;
	});

	return log_alphas;
}



Float compute_rate_naive(const std::vector<CodeWord>& transmitted, const std::vector<CodeWord>& received, 
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
			prob_table[i][j] = get_transition_prob(transmitted[i], received[j]);
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
			if (prob_table[k][j] < 1E-12)
			{
				continue;
			}
			if (std::isnan(denominator[j]) or denominator[j] < 1E-50)
			{
				denominator[j] = 1E-50;
			}
			rate += Q_i[k] * prob_table[k][j] * log(prob_table[k][j] / denominator[j]);
		}
	}

	return rate;
}