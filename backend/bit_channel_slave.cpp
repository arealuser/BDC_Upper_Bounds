// #include "bit_baa.h"
#include "bit_baa_fast.h"
#include <cstring>

const char* GENERATE_CODEWORDS = "gen_codewords";
const char* COMPUTE_DENOMS = "denominators";
const char* COMPUTE_ALPHAS = "alphas";
const char* COMPUTE_RATE = "rate";

void generate_codewords(bool up_to, size_t max_len, const char* output_file_name){
	auto codewords = get_all_bit_codewords(max_len, up_to);
	FILE* output_file = try_to_open_file(output_file_name, "wb");
	save_bit_codewords_to_file(output_file, codewords);
	fclose(output_file);
}

void compute_denominators(const char* transmitted_codewords_filename, const char* received_codewords_filename, 
	size_t start, size_t end, const char* Q_array_filename, Float deletion_probability, const char* output_file_name, 
	size_t input_len, size_t output_len, bool up_to){

	FILE* transmitted_codewords_file = try_to_open_file(transmitted_codewords_filename, "rb");
	FILE* received_codewords_file = try_to_open_file(received_codewords_filename, "rb");
	FILE* Q_array_file = try_to_open_file(Q_array_filename, "rb");
	FILE* output_file = try_to_open_file(output_file_name, "wb");

	auto transmitted_codewords = load_bit_codewords_from_file_fast(transmitted_codewords_file);
	auto received_codewords = load_bit_codewords_from_file_fast(received_codewords_file, start, end);
	auto Q = load_1d_array_from_file(Q_array_file);

	initialize_bit_channel(deletion_probability, input_len, output_len, up_to);

	auto denominators = compute_all_log_Wjk_den (transmitted_codewords, received_codewords, Q);
	write_1d_array_to_file(output_file, denominators);
	fclose(output_file); fclose(transmitted_codewords_file); fclose(received_codewords_file); fclose(Q_array_file);
}


void compute_alphas(const char* transmitted_codewords_filename, const char* received_codewords_filename,
	size_t start, size_t end, const char* Q_array_filename, const char* denominators_filename, 
	size_t input_len, size_t output_len, bool up_to, Float deletion_probability, const char* output_file_name){

	FILE* transmitted_codewords_file = try_to_open_file(transmitted_codewords_filename, "rb");
	FILE* received_codewords_file = try_to_open_file(received_codewords_filename, "rb");
	FILE* Q_array_file = try_to_open_file(Q_array_filename, "rb");
	FILE* denominators_file = try_to_open_file(denominators_filename, "rb");
	FILE* output_file = try_to_open_file(output_file_name, "wb");

	assert(start <= end);
	auto transmitted_codewords = load_bit_codewords_from_file_fast(transmitted_codewords_file, start, end);
	auto received_codewords = load_bit_codewords_from_file_fast(received_codewords_file);
	std::vector<Float> Q;
	Q.resize(end - start);
	{
		auto Q_all = load_1d_array_from_file(Q_array_file);
		assert(Q_all.size() >= end);
		for (size_t i = 0; i < (end-start); ++i)
		{
			Q[i] = Q_all[start+i];
		}
	}

	auto denominators = load_1d_array_from_file(denominators_file);

	initialize_bit_channel(deletion_probability, input_len, output_len, up_to);

	auto alphas = compute_all_log_alpha_k (transmitted_codewords, received_codewords, Q, denominators);

	write_1d_array_to_file(output_file, alphas);

	fclose(output_file); fclose(transmitted_codewords_file); fclose(received_codewords_file); 
	fclose(Q_array_file); fclose(denominators_file);
}

void compute_rate(const char* transmitted_codewords_filename, const char* received_codewords_filename,
	size_t start, size_t end, const char* Q_array_filename, const char* denominators_filename, 
	size_t input_len, size_t output_len, bool up_to, Float deletion_probability, const char* output_file_name){

	FILE* transmitted_codewords_file = try_to_open_file(transmitted_codewords_filename, "rb");
	FILE* received_codewords_file = try_to_open_file(received_codewords_filename, "rb");
	FILE* Q_array_file = try_to_open_file(Q_array_filename, "rb");
	FILE* denominators_file = try_to_open_file(denominators_filename, "rb");
	FILE* output_file = try_to_open_file(output_file_name, "wb");

	assert(start <= end);
	auto transmitted_codewords = load_bit_codewords_from_file_fast(transmitted_codewords_file, start, end);
	auto received_codewords = load_bit_codewords_from_file_fast(received_codewords_file);
	std::vector<Float> Q;
	Q.resize(end - start);
	{
		auto Q_all = load_1d_array_from_file(Q_array_file);
		assert(Q_all.size() >= end);
		for (size_t i = 0; i < (end-start); ++i)
		{
			Q[i] = Q_all[start+i];
		}
	}

	auto denominators = load_1d_array_from_file(denominators_file);

	initialize_bit_channel(deletion_probability, input_len, output_len, up_to);

	Float rate = compute_bit_rate_efficient_fast(transmitted_codewords, received_codewords, denominators, Q);
	std::vector<Float> rate_as_array = {rate};

	write_1d_array_to_file(output_file, rate_as_array);

	fclose(output_file); fclose(transmitted_codewords_file); fclose(received_codewords_file); 
	fclose(Q_array_file); fclose(denominators_file);
}


int main(int argc, char const *argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s command [options]\n", argv[0]);
		exit(1);
	}
	if (!strcmp(argv[1], GENERATE_CODEWORDS))
	{
		// Generate a file with all of the codewords. This is done so that different iterations of the computation
		//    can have the same codewords and won't need to produce them all.
		if (argc != 5)
		{
			fprintf(stderr, "Usage %s %s up_to max_len output_file_name\n", argv[0], argv[1]);
			exit(1);
		}
		bool up_to = atoi(argv[2]);
		size_t max_len = atol(argv[3]);
		const char* output_file_name = argv[4];
		generate_codewords(up_to, max_len, output_file_name);
		return 0;
	} else if(!strcmp(argv[1], COMPUTE_DENOMS)){
		// Compute the denominators of the W_jk parameters. These are used in the logarithm portion of the formula for the alphas
		//    and require a summation over all possible transmitted codewords.
		if (argc != 12)
		{
			fprintf(stderr, 
				"Usage %s %s transmitted_codewords_file received_codewords_file start end Q_array_file deletion_probability output_file input_len output_len up_to\n", 
				argv[0], argv[1]);
			exit(1);
		}
		const char* transmitted_codewords_filename = argv[2];
		const char* received_codewords_filename = argv[3];
		size_t start = atol(argv[4]);
		size_t end = atol(argv[5]);
		const char* Q_array_filename = argv[6];
		Float deletion_probability = atof(argv[7]);
		const char* output_file_name = argv[8];
		size_t input_len = atol(argv[9]);
		size_t output_len = atol(argv[10]);
		bool up_to = atoi(argv[11]);

		compute_denominators(transmitted_codewords_filename, received_codewords_filename, start, end, Q_array_filename, 
			deletion_probability, output_file_name, input_len, output_len, up_to);
	} else if(!strcmp(argv[1], COMPUTE_ALPHAS)){
		if (argc != 13)
		{
			fprintf(stderr, 
				"Usage %s %s transmitted_codewords_file received_codewords_file start end Q_array_file deletion_probability log_dens_file output_file input_len output_len up_to\n", 
				argv[0], argv[1]);
			exit(1);
		}

		const char* transmitted_codewords_filename = argv[2];
		const char* received_codewords_filename = argv[3];
		size_t start = atol(argv[4]);
		size_t end = atol(argv[5]);
		const char* Q_array_filename = argv[6];
		Float deletion_probability = atof(argv[7]);
		const char* log_dens_filename = argv[8];
		const char* output_file_name = argv[9];
		size_t input_len = atol(argv[10]);
		size_t output_len = atol(argv[11]);
		bool up_to = atoi(argv[12]);

		compute_alphas(transmitted_codewords_filename, received_codewords_filename,
			start, end, Q_array_filename, log_dens_filename, input_len, output_len, up_to,
			deletion_probability, output_file_name);

	} else if(!strcmp(argv[1], COMPUTE_RATE)){

		if (argc != 13)
		{
			fprintf(stderr, 
				"Usage %s %s transmitted_codewords_file received_codewords_file start end Q_array_file deletion_probability log_dens_file output_file input_len output_len up_to\n", 
				argv[0], argv[1]);
			exit(1);
		}

		const char* transmitted_codewords_filename = argv[2];
		const char* received_codewords_filename = argv[3];
		size_t start = atol(argv[4]);
		size_t end = atol(argv[5]);
		const char* Q_array_filename = argv[6];
		Float deletion_probability = atof(argv[7]);
		const char* log_dens_filename = argv[8];
		const char* output_file_name = argv[9];
		size_t input_len = atol(argv[10]);
		size_t output_len = atol(argv[11]);
		bool up_to = atoi(argv[12]);

		compute_rate(transmitted_codewords_filename, received_codewords_filename,
			start, end, Q_array_filename, log_dens_filename, input_len, output_len, up_to,
			deletion_probability, output_file_name);

	} else{
		fprintf(stderr, "Unknown command %s.\n", argv[1]);
		fprintf(stderr, "Try running with %s %s %s or %s instead.\n", 
			GENERATE_CODEWORDS, COMPUTE_DENOMS, COMPUTE_ALPHAS, COMPUTE_RATE);
		exit(3);
	}
	return 0;
}