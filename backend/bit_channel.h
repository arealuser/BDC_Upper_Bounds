#pragma once
#include "channel.h"
#include "cached_transition_probs.h"


typedef std::vector<uint8_t> BitCodeWord;

struct EfficientBitCodeWord;

void initialize_bit_channel(Float deletion_prob, size_t in_len, size_t out_len, bool up_to, bool load_cache=true);

template <typename _InputIter>
uint64_t btc_to_num(_InputIter first, _InputIter last);

uint64_t btc_to_num(const BitCodeWord& codeword);
BitCodeWord num_to_btc(uint64_t num, size_t len);

uint64_t btc_to_idx(const BitCodeWord& codeword);
uint64_t btc_to_idx(const EfficientBitCodeWord& codeword);
uint64_t num_to_idx(uint64_t num, size_t len);

/*
Converts a BitCodeWord into a normal CodeWord
*/
CodeWord convert_to_run_word(const BitCodeWord& bit_code);

/*
Converts a CodeWord into a BitCodeWord.
*/
BitCodeWord convert_to_bit_word(const CodeWord& run_word);


struct EfficientBitCodeWord
{
	uint64_t num;
	size_t len;
	inline EfficientBitCodeWord(uint64_t n, size_t l) : num(n), len(l) {}
	inline EfficientBitCodeWord(const BitCodeWord& word) : num(btc_to_num(word)), len(word.size()) {}
	inline EfficientBitCodeWord(const CodeWord& word) : num(btc_to_num(convert_to_bit_word(word))), len(word.total_length) {}

	/*
	Uses a smart ordering to allow a simple sort of the array of efficient bit codwords to set equivalent codewords next
		to one another for easier implementation of symmetry speed-ups.
	*/
	friend bool operator< (const EfficientBitCodeWord& a, const EfficientBitCodeWord& b);

	friend EfficientBitCodeWord operator~ (const EfficientBitCodeWord& a);
};

/*
Uses the dynamic programming algorithm to determine the probability that the transmitted code-word will
	be transformed into the recieved one by a deletion channel with deletion probability deletion_prob.
*/
Float get_bit_transition_prob(const BitCodeWord& transmitted, const BitCodeWord& recieved, bool verbose=false,
	bool use_cache=false);

Float get_bit_transition_prob_fast(const EfficientBitCodeWord& transmitted, const EfficientBitCodeWord& recieved, bool verbose=false);


/*
Uses the dynamic programming algorithm to determine the number of ways the transmitted codeword could have been transformed
	into the received one. By multiplying the result of this function with the correct normalization factor (which depends
	only on the lengths of the transmitted and received codewords), on can compute the transition probability.
This function is also used to produce a cached version of the results which allows us to compute the transition probability
	very quickly.
*/
size_t get_num_transition_possibilities(const BitCodeWord& transmitted, const BitCodeWord& recieved, bool verbose=false);


/*
Uses the cached values of the transition counts (as produced by get_num_transition_possibilities)
	on halves of transmitted codewords in order to quickly compute the transition count on whole codewords.
*/

size_t get_num_transition_possibilities_using_cache(const BitCodeWord& transmitted, const BitCodeWord& recieved, bool verbose=false);
size_t get_num_transition_possibilities_using_cache_fast(const EfficientBitCodeWord& transmitted, const EfficientBitCodeWord& recieved, bool verbose=false);




std::vector<BitCodeWord> get_all_bit_codewords(size_t len, bool up_to=false);




/*
Saves the given array of codewords to the given file.
Saves each codeword as 2 consecutive uint64_ts, where the first is the length of the codeword and the second is 
	its bitwise representation.
*/
void save_bit_codewords_to_file(FILE* out_file, const std::vector<EfficientBitCodeWord> codewords);

/*
Loads an array of codewords from the given file.
If a range of values is specified then only these values are loaded from the array.
Otherwise the entire array is loaded.
*/
std::vector<BitCodeWord> load_bit_codewords_from_file(FILE* in_file, size_t from = 0, size_t to = -1);
std::vector<EfficientBitCodeWord> load_bit_codewords_from_file_fast(FILE* in_file, size_t from = 0, size_t to = -1);

#include "bit_channel.inl"
