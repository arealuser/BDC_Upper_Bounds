#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include "utils.h"

inline std::string get_cache_filename(size_t n, size_t k){
	const std::string BASE_PATH = __SOURCE_PATH__;
	const std::string PATH_SEPARATOR = "/";
	const std::string CACHE_FOLDER_NAME = "transition_counts";
	const std::string FORMAT = "cache_";
	const std::string FORMAT_SEPARATOR = "_";

	return BASE_PATH + PATH_SEPARATOR + CACHE_FOLDER_NAME + PATH_SEPARATOR + FORMAT + std::to_string(n) + FORMAT_SEPARATOR + std::to_string(k);
}


inline void _save_data_to_cache_file(FILE* cache_file, const std::vector<Int>& data){
	fwrite(data.data(), sizeof(Int), data.size(), cache_file);
}

inline void save_data_to_cache_file(size_t n, size_t k, const std::vector<Int>& data){
	std::string filename = get_cache_filename(n, k);
	FILE* cache_file = try_to_open_file(filename.data(), "wb");
	_save_data_to_cache_file(cache_file, data);
	fclose(cache_file);
}


inline std::vector<Int> _load_data_from_cache_file(FILE* cache_file){
	std::vector<Int> v;
	Int buf[1024];
	while (size_t len = fread(buf, sizeof(Int), 1024, cache_file))
		v.insert(v.end(), buf, buf + len);
	return v;
}

inline std::vector<Int> load_data_from_cache_file(size_t n, size_t k){
	std::string filename = get_cache_filename(n, k);
	FILE* cache_file = try_to_open_file(filename.data(), "rb");
	auto res = _load_data_from_cache_file(cache_file);
	fclose(cache_file);
	return res;
}
