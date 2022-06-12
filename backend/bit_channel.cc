#include "bit_channel.h"
#include "transition_counts/cached_transition_count.h"
#include <cmath>


std::vector<std::vector<Float> > _normalization_factors;
Float _deletion_prob = 0.0;

void initialize_bit_channel(Float deletion_prob, size_t in_len, size_t out_len, bool up_to){
	_deletion_prob = deletion_prob;
	_normalization_factors.resize(in_len+1);
	for (size_t i_len = 0; i_len < in_len+1; ++i_len)
	{
		_normalization_factors[i_len].resize(out_len+1);
		std::vector<Float> prob_out_len_k;
		prob_out_len_k.resize(out_len+1);
		prob_out_len_k[0] = pow(deletion_prob, i_len);
		std::vector<size_t> num_out_len_k;
		num_out_len_k.resize(out_len+1);
		num_out_len_k[0] = 1;
		for (size_t k = 1; k < out_len+1; ++k)
		{
			num_out_len_k[k] = ((i_len - k + 1) * num_out_len_k[k-1]) / (k);
			if (up_to)
			{
				prob_out_len_k[k] = pow(deletion_prob, i_len-k) * pow(1 - deletion_prob, k);
			}
		}
		if (not up_to)
		{
			_normalization_factors[i_len][out_len] = 1. / num_out_len_k[out_len];
			continue;
		}
		Float total = std::inner_product(num_out_len_k.begin(), num_out_len_k.end(), prob_out_len_k.begin(), 0.0);
		for (size_t k = 0; k < out_len+1; ++k)
		{
			_normalization_factors[i_len][k] = prob_out_len_k[k] / total;
		}
	}
		
}

size_t get_num_transition_possibilities(const BitCodeWord& transmitted, const BitCodeWord& recieved, bool verbose){
	size_t st = transmitted.size() + 1;
	size_t sr = recieved.size() + 1;

	std::vector<size_t> dynamic_programming_state;
	dynamic_programming_state.resize(st * sr);

	for (size_t it = 0; it < st; ++it)
	{
		dynamic_programming_state[it] = 1;
	}
	for (size_t ir = 1; ir < sr; ++ir)
	{
		for (size_t it = 0; it < st; ++it)
		{
			if (it == 0)
			{
				continue;
			}
			if (it + (sr - ir) < sr)
			{
				continue;
			}

			dynamic_programming_state[(ir*st)+it] = dynamic_programming_state[(ir*st)+it-1];
			if (transmitted[it - 1] == recieved[ir - 1])
			{
				dynamic_programming_state[(ir * st) + it] += dynamic_programming_state[((ir-1) * st) + (it-1)];
			}
		}
	}
	return dynamic_programming_state[(sr * st) - 1];
}


Float get_bit_transition_prob(const BitCodeWord& transmitted, const BitCodeWord& recieved, bool verbose, bool use_cache){
	size_t st = transmitted.size() + 1;
	size_t sr = recieved.size() + 1;

	Float base_prob = _normalization_factors[st - 1][sr - 1];
	
	size_t count;
	if (verbose)
	{
		printf("use_cache=%d\n", use_cache);
	}
	if (use_cache)
	{
		count = get_num_transition_possibilities_using_cache(transmitted, recieved, verbose);
	} else {
		count = get_num_transition_possibilities(transmitted, recieved, verbose);
	}
	if (verbose)
	{
		printf("%lu, %lu\n", transmitted.size(), recieved.size());
		printf("transmitted = [");
		for (size_t i = 0; i < transmitted.size(); ++i)
		{
			printf("%d, ", transmitted[i]);
		}
		printf("]\n");

		printf("recieved = [");
		for (size_t i = 0; i < recieved.size(); ++i)
		{
			printf("%d, ", recieved[i]);
		}
		printf("]\n");
		printf("base_prob = %f, count = %lu, final_prob = %f\n", 
			base_prob, count, count * base_prob);
	}
	return base_prob * count;
}


Float get_bit_transition_prob_fast(const EfficientBitCodeWord& transmitted, const EfficientBitCodeWord& recieved, bool verbose){
	size_t st = transmitted.len + 1;
	size_t sr = recieved.len + 1;

	Float base_prob = _normalization_factors[st - 1][sr - 1];
	
	size_t count;
	count = get_num_transition_possibilities_using_cache_fast(transmitted, recieved, verbose);
	if (verbose)
	{
		printf("%lu, %lu\n", transmitted.len, recieved.len);
		printf("transmitted = %016lx\n", num_to_idx(transmitted.num, transmitted.len));

		printf("recieved = %016lx\n", num_to_idx(recieved.num, recieved.len));
		printf("base_prob = %f, count = %lu, final_prob = %f\n", 
			base_prob, count, count * base_prob);
	}
	return base_prob * count;
}


CodeWord convert_to_run_word(const BitCodeWord& bit_code){
	std::vector<Run> res;
	if (bit_code.size() == 0)
	{
		return res;
	}

	res.push_back(Run(*bit_code.begin(), 0));
	for(auto iter = bit_code.begin(); iter != bit_code.end(); ++iter){
		if (*iter == res.rbegin() -> value)
		{
			res.rbegin() -> length++;
		} else{
			res.push_back(Run(*iter, 1));
		}
	}
	return res;
}


BitCodeWord convert_to_bit_word(const CodeWord& run_word){
	BitCodeWord res;
	for(auto iter = run_word.begin(); iter != run_word.end(); ++iter){
		for (size_t i = 0; i < iter -> length; ++i)
		{
			res.push_back(iter -> value);
		}
	}
	return res;
}


std::vector<BitCodeWord> get_all_bit_codewords(size_t len, bool up_to){
	if (len == 0)
	{
		BitCodeWord empty;
		std::vector<BitCodeWord> res = {empty};
		return res;
	}
	else
	{
		std::vector<BitCodeWord> recursion = get_all_bit_codewords(len - 1);
		std::vector<BitCodeWord> res;
		for(auto iter = recursion.begin(); iter != recursion.end(); ++iter){
			if ((iter -> size()) != (len-1))
			{
				continue;
			}
			BitCodeWord word1 = *iter;
			BitCodeWord word2 = *iter;
			word1.push_back(1);
			word2.push_back(0);
			res.push_back(word1);
			res.push_back(word2);
			if (up_to)
			{
				res.push_back(*iter);
			}
		}
		return res;
	}
}



uint64_t btc_to_num(const BitCodeWord& codeword){
	uint64_t res = 0;
	for (size_t i = 0; i < codeword.size(); ++i)
	{
		res <<= 1;
		res ^= codeword[i];
	}
	return res;
}
BitCodeWord num_to_btc(uint64_t num, size_t len){
	BitCodeWord res;
	res.reserve(len);
	for (size_t i = 0; i < len; ++i)
	{
		res.push_back((num >> i) & 0x01);
	}
	return res;
}



void save_bit_codewords_to_file(FILE* out_file, const std::vector<BitCodeWord> codewords){
	for (const auto& btc : codewords)
	{
		uint64_t len = btc.size();
		uint64_t num = btc_to_num(btc);
		fwrite(&len, sizeof(len), 1, out_file);
		fwrite(&num, sizeof(num), 1, out_file);
	}
}

std::vector<BitCodeWord> load_bit_codewords_from_file(FILE* in_file, size_t from, size_t to){
	constexpr size_t buff_size = 128;
	uint64_t buffer[2*buff_size];
	std::vector<BitCodeWord> res;
	size_t num_read, total_read = 0;
	fseek(in_file, 2*sizeof(uint64_t)*from, SEEK_CUR);
	while(num_read = fread(buffer, 2*sizeof(uint64_t), std::min(buff_size, (to - from - total_read)), in_file)){
		total_read += num_read;
		for (size_t i = 0; i < num_read; ++i)
		{
			res.push_back(num_to_btc(buffer[(2*i)+1], buffer[(2*i)]));
		}
	}
	return res;
}



std::vector<EfficientBitCodeWord> load_bit_codewords_from_file_fast(FILE* in_file, size_t from, size_t to){
	constexpr size_t buff_size = 128;
	uint64_t buffer[2*buff_size];
	std::vector<EfficientBitCodeWord> res;
	size_t num_read, total_read = 0;
	fseek(in_file, 2*sizeof(uint64_t)*from, SEEK_CUR);
	while(num_read = fread(buffer, 2*sizeof(uint64_t), std::min(buff_size, (to - from - total_read)), in_file)){
		total_read += num_read;
		for (size_t i = 0; i < num_read; ++i)
		{
			res.push_back(EfficientBitCodeWord(buffer[(2*i)+1], buffer[(2*i)]));
		}
	}
	return res;
}


uint64_t btc_to_idx(const BitCodeWord& codeword){
	return num_to_idx(btc_to_num(codeword), codeword.size());
}
uint64_t num_to_idx(uint64_t num, size_t len){
	return num ^ (1ULL << len);
}




size_t get_num_transition_possibilities_using_cache(const BitCodeWord& transmitted, const BitCodeWord& recieved, bool verbose){
	size_t n = transmitted.size();
	size_t k = recieved.size();
	size_t n1 = (n+1) / 2;
	size_t n2 = n - n1;
	size_t total = 0;
	size_t trans_num1 = btc_to_num(transmitted.begin(), transmitted.begin() + n1);
	size_t trans_num2 = btc_to_num(transmitted.begin() + n1, transmitted.end());
	for (size_t i = 0; i <= k; ++i)
	{
		size_t k1 = i;
		size_t k2 = k - i;
		if ((k1 > n1) or (k2 > n2))
		{
			if (verbose)
			{
				printf("Skipping n1=%lu\tk1=%lu\tn2=%lu\tk2=%lu\n", n1, k1, n2, k2);
			}
			continue;
		}
		size_t rec_num1 = btc_to_num(recieved.begin(), recieved.begin() + k1);
		size_t rec_num2 = btc_to_num(recieved.begin() + k1, recieved.end());
		if (verbose)
		{
			printf("n1=%lu\t k1=%lu\t n2=%lu\t k2=%lu\t trans_num1=%lu\t trans_num2=%lu\t rec_num1=%lu\t rec_num2=%lu\n", n1, k1, n2, k2,
				trans_num1, trans_num2, rec_num1, rec_num2);
			fflush(stdout);
		}
		size_t count1 = get_transition_count_cache(n1, trans_num1, k1, rec_num1);
		size_t count2 = get_transition_count_cache(n2, trans_num2, k2, rec_num2);
		total += count1*count2;
		if (verbose)
		{
			printf("n=%lu\tk=%lu\tcount1=%lu\tcount2=%lutotal=%lu\n", n, k, count1, count2, total);
		}
	}

	return total;
}

size_t get_num_transition_possibilities_using_cache_fast(const EfficientBitCodeWord& transmitted, const EfficientBitCodeWord& recieved, bool verbose){
	size_t n = transmitted.len;
	size_t k = recieved.len;
	size_t n1 = (n+1) / 2;
	size_t n2 = n - n1;
	size_t total = 0;
	size_t trans_num1 = transmitted.num >> n2;
	size_t trans_num2 = transmitted.num  & ((1 << n2) - 1);
	for (size_t i = 0; i <= k; ++i)
	{
		size_t k1 = i;
		size_t k2 = k - i;
		if ((k1 > n1) or (k2 > n2))
		{
			if (verbose)
			{
				printf("Skipping n1=%lu\tk1=%lu\tn2=%lu\tk2=%lu\n", n1, k1, n2, k2);
			}
			continue;
		}
		size_t rec_num1 = recieved.num >> k2;
		size_t rec_num2 = recieved.num & ((1 << k2) - 1);
		if (verbose)
		{
			printf("n1=%lu\t k1=%lu\t n2=%lu\t k2=%lu\t trans_num1=%lu\t trans_num2=%lu\t rec_num1=%lu\t rec_num2=%lu\n", n1, k1, n2, k2,
				trans_num1, trans_num2, rec_num1, rec_num2);
			fflush(stdout);
		}
		size_t count1 = get_transition_count_cache(n1, trans_num1, k1, rec_num1);
		size_t count2 = get_transition_count_cache(n2, trans_num2, k2, rec_num2);
		total += count1*count2;
		if (verbose)
		{
			printf("count1=%lu\tcount2=%lu\ttotal=%lu\n", count1, count2, total);
		}
	}

	return total;
}



template <typename _InputIter>
uint64_t btc_to_num(_InputIter first, _InputIter last){
	if (first == last)
	{
		return 0;
	}
	size_t res = 0;
	for(; first != last; ++first){
		res <<= 1;
		res ^= *first;
	}
	return res;
}