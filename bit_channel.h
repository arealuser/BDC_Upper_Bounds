#pragma once
#include "channel.h"


typedef std::vector<uint8_t> BitCodeWord;

/*
Uses the dynamic programming algorithm to determine the probability that the transmitted code-word will
	be transformed into the recieved one by a deletion channel with deletion probability deletion_prob.
*/
Float get_bit_transition_prob(const BitCodeWord& transmitted, const BitCodeWord& recieved, Float deletion_prob, bool verbose=false);

/*
Converts a BitCodeWord into a normal CodeWord
*/
CodeWord convert_to_run_word(const BitCodeWord& bit_code);

/*
Converts a CodeWord into a BitCodeWord.
*/
BitCodeWord convert_to_bit_word(const CodeWord& run_word);

