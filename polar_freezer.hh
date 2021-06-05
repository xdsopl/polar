/*
Bit freezers for polar codes

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#pragma once

#include <algorithm>

class PolarFreezer
{
	static const int U = 2;
	static void freeze(uint8_t *bits, long double pe, long double th, int i, int h)
	{
		if (h) {
			freeze(bits, pe * (2-pe), th, i, h/2);
			freeze(bits, pe * pe, th, i+h, h/2);
		} else {
			bits[i>>U] |= (pe>th) << (i&((1<<U)-1));
		}
	}
public:
	int operator()(uint8_t *frozen_bits, int level, long double erasure_probability = 0.5L, long double freezing_threshold = 0.5L)
	{
		int length = 1 << level;
		for (int i = 0; i < 1<<(level-U); ++i)
			frozen_bits[i] = 0;
		freeze(frozen_bits, erasure_probability, freezing_threshold, 0, length / 2);
		int K = 0;
		for (int i = 0; i < length; ++i)
			K += !(frozen_bits[i>>U]&(1<<(i&((1<<U)-1))));
		return K;
	}
};

template <int MAX_M>
class PolarCodeConst0
{
	static const int U = 2;
	void compute(long double pe, int i, int h)
	{
		if (h) {
			compute(pe * (2-pe), i, h/2);
			compute(pe * pe, i+h, h/2);
		} else {
			prob[i] = pe;
		}
	}
	long double prob[1<<MAX_M];
	int index[1<<MAX_M];
public:
	void operator()(uint8_t *frozen_bits, int level, int K, long double erasure_probability = std::exp(-1.L))
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		compute(erasure_probability, 0, length / 2);
		for (int i = 0; i < length; ++i)
			index[i] = i;
		std::nth_element(index, index+K, index+length, [this](int a, int b){ return prob[a] < prob[b]; });
		for (int i = 0; i < 1<<(level-U); ++i)
			frozen_bits[i] = 0;
		for (int i = K; i < length; ++i)
			frozen_bits[index[i]>>U] |= 1 << (index[i]&((1<<U)-1));
	}
};

