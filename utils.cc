#include "utils.h"
#include <cassert>

std::vector<std::vector<Float>> load_array_from_file(FILE* in_file){
	uint32_t shape[2];
	assert(fread(shape, sizeof(uint32_t), 2, in_file) == 2);
	printf("shape = %u, %u\n", shape[0], shape[1]);

	std::vector<std::vector<Float>> res; res.resize(shape[0]);
	for(auto riter = res.begin(); riter != res.end(); ++riter){
		riter -> resize(shape[1]);
		size_t len = fread(riter -> data(), sizeof(Float), shape[1], in_file);
		assert(len == shape[1]);
	}
	return res;
}

void write_array_to_file(FILE* out_file, const std::vector<std::vector<Float>>& array){
	uint32_t shape[2];
	shape[0] = (uint32_t) array.size();
	shape[1] = (uint32_t) array[0].size();
	fwrite(shape, sizeof(*shape), 2, out_file);

	for(auto iter = array.begin(); iter != array.end(); ++iter){
		assert(fwrite(iter -> data(), sizeof(Float), shape[1], out_file) == shape[1]);
	}
}
