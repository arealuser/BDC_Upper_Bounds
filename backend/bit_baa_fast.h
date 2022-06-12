#pragma once
#include "utils.h"
#include "bit_channel.h"


/*
Performs a full BAA step on the given input and output alphabets, with the given initial distribution Q_i.
*/
std::vector<Float> do_full_baa_step(const std::vector<EfficientBitCodeWord>& transmitted, const std::vector<EfficientBitCodeWord>& received, 
	const std::vector<Float>& Q_i);


/*
Computes the amount of information from the given distribution on the given transmitted codewords.
This is a straightforward implementation that computes a full transition prob table and is useful only for debugging purposes.
*/
Float compute_rate(const std::vector<EfficientBitCodeWord>& transmitted, const std::vector<EfficientBitCodeWord>& received, 
	const std::vector<Float>& Q_i);

/*
Computes the amount of information from the given distribution on the given transmitted codewords.
For disributing purposes it is possible to run this with only some of the codewords and then to sum over the possibilities.
*/
Float compute_bit_rate_efficient_fast(const std::vector<EfficientBitCodeWord>& transmitted, const std::vector<EfficientBitCodeWord>& received, 
	const std::vector<Float>& log_W_jk_den, const std::vector<Float>& Q_i);


/*
Computes the denominator of multiple W_jk entries. 
This is a function that depends on the transition probabilities out of all of the transmitted codewords.
In general, when distributing, this should be called with all of the transmitted codewords and part of the received ones.
*/
std::vector<Float> compute_all_log_Wjk_den (const std::vector<EfficientBitCodeWord>& transmitted, const std::vector<EfficientBitCodeWord>& received, 
	const std::vector<Float>& Q_i);

/*
Computes the values of alphas (which determine the probabilities in the next BAA step).
When distributing, this should be called with a subset of the transmitted codewords and all of the received ones.
*/
std::vector<Float> compute_all_log_alpha_k (const std::vector<EfficientBitCodeWord>& transmitted, const std::vector<EfficientBitCodeWord>& received, 
	const std::vector<Float>& Q_i, const std::vector<Float>& log_W_jk_den);




std::vector<Float> compute_Pjk_row(const EfficientBitCodeWord& transmitted, const std::vector<EfficientBitCodeWord>& received);

std::vector<Float> compute_Pjk_col(const std::vector<EfficientBitCodeWord>& transmitted, const EfficientBitCodeWord& received);

Float compute_log_Wjk_den (const std::vector<EfficientBitCodeWord>& transmitted, const EfficientBitCodeWord& received, const std::vector<Float>& Q_i);
Float compute_log_alpha_k (const EfficientBitCodeWord& transmitted, const std::vector<EfficientBitCodeWord>& received, 
	Float Q_k, const std::vector<Float>& log_W_jk_den);