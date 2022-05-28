#pragma once
#include "baa.h"

/*
Computes the denominator of multiple W_jk entries parallelized over many input distributions Q_i. 
This is a function that depends on the transition probabilities out of all of the transmitted codewords.
In general, when distributing, this should be called with all of the transmitted codewords and part of the received ones.
*/
std::vector<std::vector<Float> > compute_all_log_Wjk_den_parallelized (const std::vector<CodeWord>& transmitted, const std::vector<CodeWord>& received, 
	const std::vector<std::vector<Float> >& Q_is);

/*
Computes the values of alphas (which determine the probabilities in the next BAA step).
When distributing, this should be called with a subset of the transmitted codewords and all of the received ones.
*/
std::vector<std::vector<Float> > compute_all_log_alpha_k_parallelized (const std::vector<CodeWord>& transmitted, const std::vector<CodeWord>& received, 
	const std::vector<std::vector<Float> >& Q_is, const std::vector<std::vector<Float> >& log_W_jk_den);

/*
Performs a full BAA step on the given input and output alphabets, with the given initial distribution Q_i.
*/
std::vector<std::vector<Float> > do_full_baa_step_parallelized(const std::vector<CodeWord>& transmitted, 
	const std::vector<CodeWord>& received, const std::vector<std::vector<Float> >& Q_is);

