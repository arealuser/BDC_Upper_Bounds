#include <vector>
#include <cstdio>
#include <numeric>
#include <ctime>
#include "utils.h"
#include <iostream>
#include <algorithm>

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
 
int main()
{

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
 
    std::vector<int> nums{3, 4, 2, 8, 15, 267};
    std::vector<int> nums2 = nums;

    auto print = [](const int& n) { std::cout << " " << n; };
 
    std::cout << "before:";
    std::for_each(nums.cbegin(), nums.cend(), print);
    std::cout << '\n';
 
    std::for_each(nums.begin(), nums.end(), [](int &n){ n++; });
 
    // calls Sum::operator() for each number
    Sum s = std::for_each(nums.begin(), nums.end(), Sum());
 
    std::cout << "after: ";
    std::for_each(nums.cbegin(), nums.cend(), print);
    std::cout << '\n';
    std::cout << "sum: " << s.sum << '\n';
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
