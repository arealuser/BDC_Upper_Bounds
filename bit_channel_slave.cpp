#include "bit_baa.h"
#include <cstring>

const char* GENERATE_CODEWORDS = "gen_codewords";
const char* COMPUTE_DENOMS = "denominators";


FILE* try_to_open_file(const char* filename, const char* mode){
	FILE* fptr = fopen(filename, mode);
	if (fptr == NULL)
	{
		fprintf(stderr, "Error: failed to open file %s.\n", filename);
		exit(2);
	}
	return fptr;
}

void generate_codewords(bool up_to, size_t max_len, const char* output_file_name){
	auto codewords = get_all_bit_codewords(max_len, up_to);
	FILE* output_file = try_to_open_file(output_file_name, "wb");
	save_bit_codewords_to_file(output_file, codewords);
}

void compute_denominators(const char* transmitted_codewords_filename, const char* received_codewords_filename, 
	size_t start, size_t end, const char* Q_array_filename, Float deletion_probability, const char* output_file_name, 
	size_t input_len, size_t output_len, bool up_to){

	FILE* transmitted_codewords_file = try_to_open_file(transmitted_codewords_filename, "rb");
	FILE* received_codewords_file = try_to_open_file(received_codewords_filename, "rb");
	FILE* Q_array_file = try_to_open_file(Q_array_filename, "rb");
	FILE* output_file = try_to_open_file(output_file_name, "wb");

	auto transmitted_codewords = load_bit_codewords_from_file(transmitted_codewords_file);
	auto received_codewords = load_bit_codewords_from_file(received_codewords_file, start, end);
	auto Q = load_1d_array_from_file(Q_array_file);

	initialize_bit_channel(deletion_probability, input_len, output_len, up_to);

	auto denominators = compute_all_log_Wjk_den (transmitted_codewords, received_codewords, Q);
	write_1d_array_to_file(output_file, denominators);

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
	}
	return 0;
}