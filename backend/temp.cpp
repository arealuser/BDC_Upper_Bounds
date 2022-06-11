// #include <vector>
// #include <cstdio>
// #include <numeric>
// #include <ctime>
// // #include "utils.h"
// // #include <iostream>
// #include <algorithm>

// class foo: public std::vector<int>{
// public:
// 	int sum;

// 	template <typename ...Args>
// 	foo(Args... args): std::vector<int>(args...), sum(std::accumulate(begin(), end(), 0)) {}
// };


struct Sum
{
    void operator()(int n) { sum += n; }
    int sum{0};
};
 
int main(int argc, char const *argv[])
{

    // std::array<std::array<int, 2>, 2> a = {1, 2, 3, 4};
    // printf("%d\n", a[1][1]);

    // if (argc > 1)
    // {
    //     FILE* in_file = fopen(argv[1], "r");
    //     auto arr = load_array_from_file(in_file);
    //     printf("array shape = (%lu, %lu)\n", arr.size(), arr[0].size());
    //     printf("arr =\n");
    //     auto print = [](const Float& n) { std::cout << " " << n; };
    //     auto print_vec = [print](const std::vector<Float>& v) { std::for_each(v.begin(), v.end(), print); std::cout << std::endl;};
    //     std::for_each(arr.begin(), arr.end(), print_vec);
    //     if (argc > 2)
    //     {
    //         FILE* out_file = fopen(argv[2], "w");
    //         write_array_to_file(out_file, arr);
    //     }
    // }
    

    // auto t1 = clock();
    // double value = 1.00000001;
    // volatile double res = 1.0;
    // for (size_t i = 0; i < 1000000000ULL; ++i)
    // {
    //     res *= value;
    // }
    // printf("%f\n", res);
    // printf("This took %.1f seconds.\n", ((float) (clock() - t1)) / CLOCKS_PER_SEC); fflush(stdout);
    // 

    // for (auto num : nums)
    // {
    // 	printf("%d\n", num);
    // }
 
    // std::vector<int> nums{3, 4, 2, 8, 15, 267};
    // std::vector<int> nums2 = nums;

    // auto print = [](const int& n) { std::cout << " " << n; };
 
    // std::cout << "before:";
    // std::for_each(nums.cbegin(), nums.cend(), print);
    // std::cout << '\n';
 
    // std::for_each(nums.begin(), nums.end(), [](int &n){ n++; });
 
    // // calls Sum::operator() for each number
    // Sum s = std::for_each(nums.begin(), nums.end(), Sum());
 
    // std::cout << "after: ";
    // std::for_each(nums.cbegin(), nums.cend(), print);
    // std::cout << '\n';
    // std::cout << "sum: " << s.sum << '\n';
}

// int main(int argc, char const *argv[])
// {

// 	// foo f = (std::vector<int>({1, 2, 3, 4}));
// 	// foo f;
// 	// f.push_back(1);
// 	// f.push_back(2);
// 	// f.push_back()
// 	// printf("%d\n", f.sum);
// 	return 0;
// }
