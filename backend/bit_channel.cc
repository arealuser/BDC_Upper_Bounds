#include "bit_channel.h"
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






void save_bit_codewords_to_file(FILE* out_file, const std::vector<EfficientBitCodeWord> codewords){
	for (const auto& btc : codewords)
	{
		uint64_t len = btc.len;
		uint64_t num = btc.num;
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

