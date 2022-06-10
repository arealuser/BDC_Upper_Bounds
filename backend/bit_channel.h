#pragma once
#include "channel.h"


typedef std::vector<uint8_t> BitCodeWord;

/*
Uses the dynamic programming algorithm to determine the probability that the transmitted code-word will
	be transformed into the recieved one by a deletion channel with deletion probability deletion_prob.
*/
Float get_bit_transition_prob(const BitCodeWord& transmitted, const BitCodeWord& recieved, bool verbose=false);

/*
Converts a BitCodeWord into a normal CodeWord
*/
CodeWord convert_to_run_word(const BitCodeWord& bit_code);

/*
Converts a CodeWord into a BitCodeWord.
*/
BitCodeWord convert_to_bit_word(const CodeWord& run_word);


std::vector<BitCodeWord> get_all_bit_codewords(size_t len, bool up_to=false);

void initialize_bit_channel(Float deletion_prob, size_t in_len, size_t out_len, bool up_to);

uint64_t btc_to_num(const BitCodeWord& codeword);
BitCodeWord num_to_btc(uint64_t num, size_t len);


/*
Saves the given array of codewords to the given file.
Saves each codeword as 2 consecutive uint64_ts, where the first is the length of the codeword and the second is 
	its bitwise representation.
*/
void save_bit_codewords_to_file(FILE* out_file, const std::vector<BitCodeWord> codewords);

/*
Loads an array of codewords from the given file.
If a range of values is specified then only these values are loaded from the array.
Otherwise the entire array is loaded.
*/
std::vector<BitCodeWord> load_bit_codewords_from_file(FILE* in_file, size_t from = 0, size_t to = -1);