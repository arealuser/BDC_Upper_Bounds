#pragma once

typedef double Float;
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdio>

std::vector<std::vector<Float>> load_array_from_file(FILE* in_file);
void write_array_to_file(FILE* out_file, const std::vector<std::vector<Float>>& array);
