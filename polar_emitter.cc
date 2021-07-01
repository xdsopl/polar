/*
Code emitter for successive cancellation decoding of polar codes

Copyright 2021 Ahmet Inan <xdsopl@gmail.com>
*/

#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <functional>
#include "polar_freezer.hh"

class PolarEmitter
{
	static int frozen_count(const uint8_t *frozen, int level)
	{
		int count = 0;
		for (int i = 0; i < (1<<level); ++i)
			count += frozen[i];
		return count;
	}
	static void compile(const uint8_t *frozen, int level)
	{
		assert(level > 0);
		int count = frozen_count(frozen, level);
		if (count == 1<<level) {
			std::cout << "rate0<" << level << ">(sft, hrd, msg);" << std::endl;
		} else if (count == 0) {
			std::cout << "rate1<" << level << ">(sft, hrd, msg);" << std::endl;
			std::cout << "msg += 1 << " << level << ";" << std::endl;
		} else if (count == (1<<level)-1 && !frozen[(1<<level)-1]) {
			std::cout << "rep<" << level << ">(sft, hrd, msg);" << std::endl;
			std::cout << "msg += 1;" << std::endl;
		} else if (count == 1 && frozen[0]) {
			std::cout << "spc<" << level << ">(sft, hrd, msg);" << std::endl;
			std::cout << "msg += (1 << " << level << ") - 1;" << std::endl;
		} else {
			std::cout << "left<" << level << ">(sft, hrd, msg);" << std::endl;
			compile(frozen, level-1);
			std::cout << "right<" << level << ">(sft, hrd, msg);" << std::endl;
			std::cout << "hrd += 1 << " << level-1 << ";" << std::endl;
			compile(frozen+(1<<(level-1)), level-1);
			std::cout << "hrd -= 1 << " << level-1 << ";" << std::endl;
			std::cout << "comb<" << level << ">(sft, hrd, msg);" << std::endl;
		}
	}
public:
	void operator()(const uint8_t *frozen, int level)
	{
		compile(frozen, level);
	}
};

int main()
{
	const int M = 14;
	const int N = 1 << M;
	auto frozen = new uint8_t[N];
	double erasure_probability = 0.5;
	int K = (1 - erasure_probability) * N;
	double design_SNR = 10 * std::log10(-std::log(erasure_probability));
	std::cerr << "design SNR: " << design_SNR << std::endl;
	auto freeze = new PolarCodeConst0<M>;
	double better_SNR = design_SNR + 1.59175;
	std::cerr << "better SNR: " << better_SNR << std::endl;
	double probability = std::exp(-pow(10.0, better_SNR / 10));
	(*freeze)(frozen, M, K, probability);
	PolarEmitter emit;
	emit(frozen, M);
	return 0;
}

