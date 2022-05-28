#include "parallelized_baa.h"

std::vector<std::vector<Float> > compute_all_log_Wjk_den_parallelized (const std::vector<CodeWord>& transmitted, const std::vector<CodeWord>& received, 
	const std::vector<std::vector<Float> >& Q_is){
	size_t n_Q = Q_is.size();
	// size_t n_I = transmitted.size();
	size_t n_J = received.size();

	std::vector<std::vector<Float> > res; res.resize(n_Q);
	for(auto riter = res.begin(); riter != res.end(); ++riter){
		riter -> resize(n_J);
	}

	for (size_t i_J = 0; i_J < n_J; ++i_J)
	{
		std::vector<Float> probs_col = compute_Pjk_col(transmitted, received[i_J]);
		for (size_t i_Q = 0; i_Q < n_Q; ++i_Q)
		{
			res[i_Q][i_J] = log(std::inner_product(probs_col.begin(), probs_col.end(), Q_is[i_Q].begin(), 0.0));
		}
	}
	return res;
}


std::vector<std::vector<Float> > compute_all_log_alpha_k_parallelized (const std::vector<CodeWord>& transmitted, const std::vector<CodeWord>& received, 
	const std::vector<std::vector<Float> >& Q_is, const std::vector<std::vector<Float> >& log_W_jk_den){
	size_t n_Q = Q_is.size();
	size_t n_I = transmitted.size();
	size_t n_J = received.size();

	std::vector<std::vector<Float> > log_alphass; log_alphass.resize(n_Q);
	for(auto riter = log_alphass.begin(); riter != log_alphass.end(); ++riter){
		riter -> resize(n_I);
		for(auto riter2 = (*riter).begin(); riter2 != (*riter).end(); ++riter2){
			*riter2 = 0.0;
		}
	}

	for (size_t i_I = 0; i_I < n_I; ++i_I)
	{
		std::vector<Float> probs_row = compute_Pjk_row(transmitted[i_I], received);
		std::vector<Float> log_probs_row = probs_row;
		for_each(log_probs_row.begin(), log_probs_row.end(), [](Float& x){x = log(x);});
		for (size_t i_Q = 0; i_Q < n_Q; ++i_Q)
		{
			Float log_Qk = log(Q_is[i_Q][i_I]);
			if (std::isnan(log_Qk))
			{
				log_Qk = -700.0;
			}
			for (size_t i_J = 0; i_J < n_J; ++i_J)
			{
				Float P_jk = probs_row[i_J];
				if (P_jk < 1E-12)
				{
					continue;
				}
				Float log_den = log_W_jk_den[i_Q][i_J];
				Float log_P_jk = log_probs_row[i_J];
				log_alphass[i_Q][i_I] += P_jk * (log_P_jk + log_Qk - log_den);
			}
		}
	}

	return log_alphass;

}


std::vector<std::vector<Float> > do_full_baa_step_parallelized(const std::vector<CodeWord>& transmitted, 
	const std::vector<CodeWord>& received, const std::vector<std::vector<Float> >& Q_is){
	auto log_W_jk_den = compute_all_log_Wjk_den_parallelized(transmitted, received, Q_is);
	auto log_alphass = compute_all_log_alpha_k_parallelized(transmitted, received, Q_is, log_W_jk_den);

	size_t n_Q = Q_is.size();
	for (size_t i_Q = 0; i_Q < n_Q; ++i_Q)
	{
		// std::vector<Float>& log_alphas = log_alphass[i_Q];
		Float max_log_alpha = *std::max_element(log_alphass[i_Q].begin(), log_alphass[i_Q].end());
		std::for_each(log_alphass[i_Q].begin(), log_alphass[i_Q].end(), 
			[max_log_alpha](Float& log_alpha){log_alpha = exp(log_alpha - max_log_alpha);});
		Float alpha_total = std::accumulate(log_alphass[i_Q].begin(), log_alphass[i_Q].end(), 0.0);
		std::for_each(log_alphass[i_Q].begin(), log_alphass[i_Q].end(), 
			[alpha_total](Float& log_alpha){log_alpha /= alpha_total;});
	}
	return log_alphass;
}


std::vector<Float> compute_rate_parallelized(const std::vector<CodeWord>& transmitted, const std::vector<CodeWord>& received, 
	const std::vector<std::vector<Float> >& Q_is){
	size_t n_Q = Q_is.size();
	size_t n_I = transmitted.size();
	size_t n_J = received.size();

	std::vector<Float> res; res.resize(n_Q);
	auto log_W_jk_den = compute_all_log_Wjk_den_parallelized(transmitted, received, Q_is);

	for (size_t i_I = 0; i_I < n_I; ++i_I)
	{
		std::vector<Float> probs_row = compute_Pjk_row(transmitted[i_I], received);
		std::vector<Float> log_probs_row = probs_row;
		for_each(log_probs_row.begin(), log_probs_row.end(), [](Float& x){x = log(x);});
		for (size_t i_Q = 0; i_Q < n_Q; ++i_Q)
		{
			Float Qk = Q_is[i_Q][i_I];
			// Float log_Qk = log(Qk);
			for (size_t i_J = 0; i_J < n_J; ++i_J)
			{
				Float P_jk = probs_row[i_J];
				if (P_jk < 1E-12)
				{
					continue;
				}
				Float log_den = log_W_jk_den[i_Q][i_J];
				Float log_P_jk = log_probs_row[i_J];
				res[i_Q] += Qk * P_jk * (log_P_jk - log_den);
			}
		}
	}
	return res;
}
