/*
Test bench for successive cancellation decoding of polar codes on the binary erasure channel

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#include <random>
#include <chrono>
#include <cassert>
#include <iostream>
#include <functional>

template <int M>
class PolarTransform
{
	static const int N = 1 << M;
public:
	void operator()(int8_t *output, const int8_t *input)
	{
		for (int i = 0; i < N; i += 2) {
			int in0 = input[i];
			int in1 = input[i+1];
			output[i] = in0 * in1;
			output[i+1] = in1;
		}
		for (int h = 2; h < N; h *= 2)
			for (int i = 0; i < N; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					output[j] *= output[j+h];
	}
};

template <int M>
class PolarEncoder
{
	static const int N = 1 << M;
	static const int U = 2;
public:
	void operator()(int8_t *codeword, const int8_t *message, const uint8_t *frozen)
	{
		for (int i = 0; i < N; i += 2) {
			int msg0 = frozen[i>>U]&(1<<(i&((1<<U)-1))) ? 1 : *message++;
			int msg1 = frozen[(i+1)>>U]&(1<<((i+1)&((1<<U)-1))) ? 1 : *message++;
			codeword[i] = msg0 * msg1;
			codeword[i+1] = msg1;
		}
		for (int h = 2; h < N; h *= 2)
			for (int i = 0; i < N; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					codeword[j] *= codeword[j+h];
	}
};

template <int M>
class PolarSysEnc
{
	static const int N = 1 << M;
	static const int U = 2;
public:
	void operator()(int8_t *codeword, const int8_t *message, const uint8_t *frozen)
	{
		for (int i = 0; i < N; i += 2) {
			int msg0 = frozen[i>>U]&(1<<(i&((1<<U)-1))) ? 1 : *message++;
			int msg1 = frozen[(i+1)>>U]&(1<<((i+1)&((1<<U)-1))) ? 1 : *message++;
			codeword[i] = msg0 * msg1;
			codeword[i+1] = msg1;
		}
		for (int h = 2; h < N; h *= 2)
			for (int i = 0; i < N; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					codeword[j] *= codeword[j+h];
		for (int i = 0; i < N; i += 2) {
			int msg0 = frozen[i>>U]&(1<<(i&((1<<U)-1))) ? 1 : codeword[i];
			int msg1 = frozen[(i+1)>>U]&(1<<((i+1)&((1<<U)-1))) ? 1 : codeword[i+1];
			codeword[i] = msg0 * msg1;
			codeword[i+1] = msg1;
		}
		for (int h = 2; h < N; h *= 2)
			for (int i = 0; i < N; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					codeword[j] *= codeword[j+h];
	}
};

class PolarFreezer
{
	static const int U = 2;
	static void freeze(uint8_t *bits, double pe, double th, int i, int h)
	{
		if (h) {
			freeze(bits, pe * (2-pe), th, i, h/2);
			freeze(bits, pe * pe, th, i+h, h/2);
		} else {
			bits[i>>U] |= (pe>th) << (i&((1<<U)-1));
		}
	}
public:
	int operator()(uint8_t *frozen_bits, int level, double erasure_probability = 0.5, double freezing_threshold = 0.5)
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

template<typename TYPE>
int popcnt(TYPE x)
{
	int cnt = 0;
	while (x) {
		++cnt;
		x &= x-1;
	}
	return cnt;
}

class PolarCompiler
{
	static const int U = 2;
	static uint8_t leaf(int frozen)
	{
		return frozen;
	}
	static uint8_t node(int func, int level)
	{
		return (func << 5) | level;
	}
	static uint8_t left(int level)
	{
		return node(1, level);
	}
	static uint8_t right(int level)
	{
		return node(2, level);
	}
	static uint8_t combine(int level)
	{
		return node(3, level);
	}
	static uint8_t rate0(int level)
	{
		return node(4, level);
	}
	static uint8_t rate1(int level)
	{
		return node(5, level);
	}
	static uint8_t rep(int level)
	{
		return node(6, level);
	}
	static int frozen_count(const uint8_t *frozen, int level)
	{
		int count = 0;
		for (int i = 0; i < 1<<(level-U); ++i)
			count += popcnt(frozen[i]);
		return count;
	}
	static void compile(uint8_t **program, const uint8_t *frozen, int level)
	{
		if (level > U) {
			int count = frozen_count(frozen, level);
			if (count == 1<<level) {
				*(*program)++ = rate0(level);
			} else if (count == 0) {
				*(*program)++ = rate1(level);
			} else if (count == (1<<level)-1 && frozen[(1<<(level-U))-1] == (1<<((1<<U)-1))-1) {
				*(*program)++ = rep(level);
			} else {
				*(*program)++ = left(level);
				compile(program, frozen, level-1);
				*(*program)++ = right(level);
				compile(program, frozen+(1<<(level-1-U)), level-1);
				*(*program)++ = combine(level);
			}
		} else {
			*(*program)++ = leaf(*frozen);
		}
	}
public:
	int operator()(uint8_t *program, const uint8_t *frozen, int level)
	{
		uint8_t *first = program;
		*program++ = level;
		compile(&program, frozen, level);
		*program++ = 255;
		return program - first;
	}
};

template <int MAX_M>
class PolarDecoder
{
	static const int U = 2;
	static int8_t signum(int8_t v)
	{
		return (v > 0) - (v < 0);
	}
	static int8_t qabs(int8_t a)
	{
		return std::abs(std::max<int8_t>(a, -127));
	}
	static int8_t qadd(int8_t a, int8_t b)
	{
		return std::min<int16_t>(std::max<int16_t>(int16_t(a) + int16_t(b), -128), 127);
	}
	static int8_t qmul(int8_t a, int8_t b)
	{
		return std::min<int16_t>(std::max<int16_t>(int16_t(a) * int16_t(b), -128), 127);
	}
	static int8_t prod(int8_t a, int8_t b)
	{
		return signum(a) * signum(b) * std::min(qabs(a), qabs(b));
	}
	static int8_t madd(int8_t a, int8_t b, int8_t c)
	{
		return qadd(qmul(a, b), c);
	}
	void leaf0(int8_t **msg, int index)
	{
#if 0
		int8_t soft20 = soft[(1<<U)+0];
		int8_t soft21 = soft[(1<<U)+1];
		int8_t soft22 = soft[(1<<U)+2];
		int8_t soft23 = soft[(1<<U)+3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t soft00 = prod(soft10, soft11);
		int8_t hard0 = signum(soft00);
		*(*msg)++ = hard0;
		int8_t soft01 = madd(hard0, soft10, soft11);
		int8_t hard1 = signum(soft01);
		*(*msg)++ = hard1;
		hard0 *= hard1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t soft02 = prod(soft12, soft13);
		int8_t hard2 = signum(soft02);
		*(*msg)++ = hard2;
		int8_t soft03 = madd(hard2, soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = hard3;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
#else
		int8_t hard0 = signum(soft[(1<<U)+0]);
		int8_t hard1 = signum(soft[(1<<U)+1]);
		int8_t hard2 = signum(soft[(1<<U)+2]);
		int8_t hard3 = signum(soft[(1<<U)+3]);
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
		hard0 *= hard1;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		*(*msg)++ = hard0;
		*(*msg)++ = hard1;
		*(*msg)++ = hard2;
		*(*msg)++ = hard3;
#endif
	}
	void leaf1(int8_t **msg, int index)
	{
		int8_t soft20 = soft[(1<<U)+0];
		int8_t soft21 = soft[(1<<U)+1];
		int8_t soft22 = soft[(1<<U)+2];
		int8_t soft23 = soft[(1<<U)+3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t hard0 = 1;
		int8_t soft01 = qadd(soft10, soft11);
		int8_t hard1 = signum(soft01);
		*(*msg)++ = hard1;
		hard0 *= hard1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t soft02 = prod(soft12, soft13);
		int8_t hard2 = signum(soft02);
		*(*msg)++ = hard2;
		int8_t soft03 = madd(hard2, soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = hard3;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf2(int8_t **msg, int index)
	{
		int8_t soft20 = soft[(1<<U)+0];
		int8_t soft21 = soft[(1<<U)+1];
		int8_t soft22 = soft[(1<<U)+2];
		int8_t soft23 = soft[(1<<U)+3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t soft00 = prod(soft10, soft11);
		int8_t hard0 = signum(soft00);
		*(*msg)++ = hard0;
		int8_t hard1 = 1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t soft02 = prod(soft12, soft13);
		int8_t hard2 = signum(soft02);
		*(*msg)++ = hard2;
		int8_t soft03 = madd(hard2, soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = hard3;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf3(int8_t **msg, int index)
	{
		int8_t soft20 = soft[(1<<U)+0];
		int8_t soft21 = soft[(1<<U)+1];
		int8_t soft22 = soft[(1<<U)+2];
		int8_t soft23 = soft[(1<<U)+3];
		int8_t hard0 = 1;
		int8_t hard1 = 1;
		int8_t soft12 = qadd(soft20, soft22);
		int8_t soft13 = qadd(soft21, soft23);
		int8_t soft02 = prod(soft12, soft13);
		int8_t hard2 = signum(soft02);
		*(*msg)++ = hard2;
		int8_t soft03 = madd(hard2, soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = hard3;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf4(int8_t **msg, int index)
	{
		int8_t soft20 = soft[(1<<U)+0];
		int8_t soft21 = soft[(1<<U)+1];
		int8_t soft22 = soft[(1<<U)+2];
		int8_t soft23 = soft[(1<<U)+3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t soft00 = prod(soft10, soft11);
		int8_t hard0 = signum(soft00);
		*(*msg)++ = hard0;
		int8_t soft01 = madd(hard0, soft10, soft11);
		int8_t hard1 = signum(soft01);
		*(*msg)++ = hard1;
		hard0 *= hard1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t hard2 = 1;
		int8_t soft03 = qadd(soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = hard3;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf5(int8_t **msg, int index)
	{
		int8_t soft20 = soft[(1<<U)+0];
		int8_t soft21 = soft[(1<<U)+1];
		int8_t soft22 = soft[(1<<U)+2];
		int8_t soft23 = soft[(1<<U)+3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t hard0 = 1;
		int8_t soft01 = qadd(soft10, soft11);
		int8_t hard1 = signum(soft01);
		*(*msg)++ = hard1;
		hard0 *= hard1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t hard2 = 1;
		int8_t soft03 = qadd(soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = hard3;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf6(int8_t **msg, int index)
	{
		int8_t soft20 = soft[(1<<U)+0];
		int8_t soft21 = soft[(1<<U)+1];
		int8_t soft22 = soft[(1<<U)+2];
		int8_t soft23 = soft[(1<<U)+3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t soft00 = prod(soft10, soft11);
		int8_t hard0 = signum(soft00);
		*(*msg)++ = hard0;
		int8_t hard1 = 1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t hard2 = 1;
		int8_t soft03 = qadd(soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = hard3;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf7(int8_t **msg, int index)
	{
#if 0
		int8_t soft20 = soft[(1<<U)+0];
		int8_t soft21 = soft[(1<<U)+1];
		int8_t soft22 = soft[(1<<U)+2];
		int8_t soft23 = soft[(1<<U)+3];
		int8_t hard0 = 1;
		int8_t hard1 = 1;
		int8_t soft12 = qadd(soft20, soft22);
		int8_t soft13 = qadd(soft21, soft23);
		int8_t hard2 = 1;
		int8_t soft03 = qadd(soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = hard3;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
#else
		int8_t soft20 = soft[(1<<U)+0];
		int8_t soft21 = soft[(1<<U)+1];
		int8_t soft22 = soft[(1<<U)+2];
		int8_t soft23 = soft[(1<<U)+3];
		int8_t soft12 = qadd(soft20, soft22);
		int8_t soft13 = qadd(soft21, soft23);
		int8_t soft03 = qadd(soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = hard3;
		hard[index+0] = hard3;
		hard[index+1] = hard3;
		hard[index+2] = hard3;
		hard[index+3] = hard3;
#endif
	}
	void leaf8(int8_t **msg, int index)
	{
		int8_t soft20 = soft[(1<<U)+0];
		int8_t soft21 = soft[(1<<U)+1];
		int8_t soft22 = soft[(1<<U)+2];
		int8_t soft23 = soft[(1<<U)+3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t soft00 = prod(soft10, soft11);
		int8_t hard0 = signum(soft00);
		*(*msg)++ = hard0;
		int8_t soft01 = madd(hard0, soft10, soft11);
		int8_t hard1 = signum(soft01);
		*(*msg)++ = hard1;
		hard0 *= hard1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t soft02 = prod(soft12, soft13);
		int8_t hard2 = signum(soft02);
		*(*msg)++ = hard2;
		int8_t hard3 = 1;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf9(int8_t **msg, int index)
	{
		int8_t soft20 = soft[(1<<U)+0];
		int8_t soft21 = soft[(1<<U)+1];
		int8_t soft22 = soft[(1<<U)+2];
		int8_t soft23 = soft[(1<<U)+3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t hard0 = 1;
		int8_t soft01 = qadd(soft10, soft11);
		int8_t hard1 = signum(soft01);
		*(*msg)++ = hard1;
		hard0 *= hard1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t soft02 = prod(soft12, soft13);
		int8_t hard2 = signum(soft02);
		*(*msg)++ = hard2;
		int8_t hard3 = 1;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf10(int8_t **msg, int index)
	{
		int8_t soft20 = soft[(1<<U)+0];
		int8_t soft21 = soft[(1<<U)+1];
		int8_t soft22 = soft[(1<<U)+2];
		int8_t soft23 = soft[(1<<U)+3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t soft00 = prod(soft10, soft11);
		int8_t hard0 = signum(soft00);
		*(*msg)++ = hard0;
		int8_t hard1 = 1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t soft02 = prod(soft12, soft13);
		int8_t hard2 = signum(soft02);
		*(*msg)++ = hard2;
		int8_t hard3 = 1;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf11(int8_t **msg, int index)
	{
		int8_t soft20 = soft[(1<<U)+0];
		int8_t soft21 = soft[(1<<U)+1];
		int8_t soft22 = soft[(1<<U)+2];
		int8_t soft23 = soft[(1<<U)+3];
		int8_t hard0 = 1;
		int8_t hard1 = 1;
		int8_t soft12 = qadd(soft20, soft22);
		int8_t soft13 = qadd(soft21, soft23);
		int8_t soft02 = prod(soft12, soft13);
		int8_t hard2 = signum(soft02);
		*(*msg)++ = hard2;
		int8_t hard3 = 1;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf12(int8_t **msg, int index)
	{
		int8_t soft20 = soft[(1<<U)+0];
		int8_t soft21 = soft[(1<<U)+1];
		int8_t soft22 = soft[(1<<U)+2];
		int8_t soft23 = soft[(1<<U)+3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t soft00 = prod(soft10, soft11);
		int8_t hard0 = signum(soft00);
		*(*msg)++ = hard0;
		int8_t soft01 = madd(hard0, soft10, soft11);
		int8_t hard1 = signum(soft01);
		*(*msg)++ = hard1;
		hard0 *= hard1;
		int8_t hard2 = 1;
		int8_t hard3 = 1;
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf13(int8_t **msg, int index)
	{
		int8_t soft20 = soft[(1<<U)+0];
		int8_t soft21 = soft[(1<<U)+1];
		int8_t soft22 = soft[(1<<U)+2];
		int8_t soft23 = soft[(1<<U)+3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t hard0 = 1;
		int8_t soft01 = qadd(soft10, soft11);
		int8_t hard1 = signum(soft01);
		*(*msg)++ = hard1;
		hard0 *= hard1;
		int8_t hard2 = 1;
		int8_t hard3 = 1;
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf14(int8_t **msg, int index)
	{
		int8_t soft20 = soft[(1<<U)+0];
		int8_t soft21 = soft[(1<<U)+1];
		int8_t soft22 = soft[(1<<U)+2];
		int8_t soft23 = soft[(1<<U)+3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t soft00 = prod(soft10, soft11);
		int8_t hard0 = signum(soft00);
		*(*msg)++ = hard0;
		int8_t hard1 = 1;
		int8_t hard2 = 1;
		int8_t hard3 = 1;
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf15(int8_t **, int index)
	{
		int8_t hard0 = 1;
		int8_t hard1 = 1;
		int8_t hard2 = 1;
		int8_t hard3 = 1;
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	template <int level>
	void left(int8_t **, int)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			soft[i+length/2] = prod(soft[i+length], soft[i+length/2+length]);
	}
	template <int level>
	void right(int8_t **, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			soft[i+length/2] = madd(hard[index+i], soft[i+length], soft[i+length/2+length]);
	}
	template <int level>
	void combine(int8_t **, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			hard[index+i] *= hard[index+i+length/2];
	}
	template <int level>
	void rate0(int8_t **, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length; ++i)
			hard[index+i] = 1;
	}
	template <int level>
	void rate1(int8_t **msg, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length; ++i)
			hard[index+i] = signum(soft[i+length]);
		for (int i = 0; i < length; i += 2) {
			soft[i] = hard[index+i] * hard[index+i+1];
			soft[i+1] = hard[index+i+1];
		}
		for (int h = 2; h < length; h *= 2)
			for (int i = 0; i < length; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					soft[j] *= soft[j+h];
		for (int i = 0; i < length; ++i)
			*(*msg)++ = soft[i];
	}
	template <int level>
	void rep(int8_t **msg, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int h = length; h; h /= 2)
			for (int i = 0; i < h/2; ++i)
				soft[i+h/2] = qadd(soft[i+h], soft[i+h/2+h]);
		int hardi = signum(soft[1]);
		*(*msg)++ = hardi;
		for (int i = 0; i < length; ++i)
			hard[index+i] = hardi;
	}
	int8_t soft[1U<<(MAX_M+1)];
	int8_t hard[1U<<MAX_M];
public:
	void operator()(int8_t *message, const int8_t *codeword, const uint8_t *program)
	{
		int level = *program++;
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length; ++i)
			soft[i+length] = codeword[i];
		int idx = 0;
		int8_t **msg = &message;
		while (*program != 255) {
			switch (*program++) {
			case 0: leaf0(msg, idx); break;
			case 1: leaf1(msg, idx); break;
			case 2: leaf2(msg, idx); break;
			case 3: leaf3(msg, idx); break;
			case 4: leaf4(msg, idx); break;
			case 5: leaf5(msg, idx); break;
			case 6: leaf6(msg, idx); break;
			case 7: leaf7(msg, idx); break;
			case 8: leaf8(msg, idx); break;
			case 9: leaf9(msg, idx); break;
			case 10: leaf10(msg, idx); break;
			case 11: leaf11(msg, idx); break;
			case 12: leaf12(msg, idx); break;
			case 13: leaf13(msg, idx); break;
			case 14: leaf14(msg, idx); break;
			case 15: leaf15(msg, idx); break;
			case (1<<5)+3: left<3>(msg, idx); break;
			case (1<<5)+4: left<4>(msg, idx); break;
			case (1<<5)+5: left<5>(msg, idx); break;
			case (1<<5)+6: left<6>(msg, idx); break;
			case (1<<5)+7: left<7>(msg, idx); break;
			case (1<<5)+8: left<8>(msg, idx); break;
			case (1<<5)+9: left<9>(msg, idx); break;
			case (1<<5)+10: left<10>(msg, idx); break;
			case (1<<5)+11: left<11>(msg, idx); break;
			case (1<<5)+12: left<12>(msg, idx); break;
			case (1<<5)+13: left<13>(msg, idx); break;
			case (1<<5)+14: left<14>(msg, idx); break;
			case (1<<5)+15: left<15>(msg, idx); break;
			case (1<<5)+16: left<16>(msg, idx); break;
			case (1<<5)+17: left<17>(msg, idx); break;
			case (1<<5)+18: left<18>(msg, idx); break;
			case (1<<5)+19: left<19>(msg, idx); break;
			case (1<<5)+20: left<20>(msg, idx); break;
			case (1<<5)+21: left<21>(msg, idx); break;
			case (1<<5)+22: left<22>(msg, idx); break;
			case (1<<5)+23: left<23>(msg, idx); break;
			case (1<<5)+24: left<24>(msg, idx); break;
			case (1<<5)+25: left<25>(msg, idx); break;
			case (1<<5)+26: left<26>(msg, idx); break;
			case (1<<5)+27: left<27>(msg, idx); break;
			case (1<<5)+28: left<28>(msg, idx); break;
			case (1<<5)+29: left<29>(msg, idx); break;
			case (1<<5)+30: left<30>(msg, idx); break;
			case (2<<5)+3: right<3>(msg, idx); idx += 1<<(3-1); break;
			case (2<<5)+4: right<4>(msg, idx); idx += 1<<(4-1); break;
			case (2<<5)+5: right<5>(msg, idx); idx += 1<<(5-1); break;
			case (2<<5)+6: right<6>(msg, idx); idx += 1<<(6-1); break;
			case (2<<5)+7: right<7>(msg, idx); idx += 1<<(7-1); break;
			case (2<<5)+8: right<8>(msg, idx); idx += 1<<(8-1); break;
			case (2<<5)+9: right<9>(msg, idx); idx += 1<<(9-1); break;
			case (2<<5)+10: right<10>(msg, idx); idx += 1<<(10-1); break;
			case (2<<5)+11: right<11>(msg, idx); idx += 1<<(11-1); break;
			case (2<<5)+12: right<12>(msg, idx); idx += 1<<(12-1); break;
			case (2<<5)+13: right<13>(msg, idx); idx += 1<<(13-1); break;
			case (2<<5)+14: right<14>(msg, idx); idx += 1<<(14-1); break;
			case (2<<5)+15: right<15>(msg, idx); idx += 1<<(15-1); break;
			case (2<<5)+16: right<16>(msg, idx); idx += 1<<(16-1); break;
			case (2<<5)+17: right<17>(msg, idx); idx += 1<<(17-1); break;
			case (2<<5)+18: right<18>(msg, idx); idx += 1<<(18-1); break;
			case (2<<5)+19: right<19>(msg, idx); idx += 1<<(19-1); break;
			case (2<<5)+20: right<20>(msg, idx); idx += 1<<(20-1); break;
			case (2<<5)+21: right<21>(msg, idx); idx += 1<<(21-1); break;
			case (2<<5)+22: right<22>(msg, idx); idx += 1<<(22-1); break;
			case (2<<5)+23: right<23>(msg, idx); idx += 1<<(23-1); break;
			case (2<<5)+24: right<24>(msg, idx); idx += 1<<(24-1); break;
			case (2<<5)+25: right<25>(msg, idx); idx += 1<<(25-1); break;
			case (2<<5)+26: right<26>(msg, idx); idx += 1<<(26-1); break;
			case (2<<5)+27: right<27>(msg, idx); idx += 1<<(27-1); break;
			case (2<<5)+28: right<28>(msg, idx); idx += 1<<(28-1); break;
			case (2<<5)+29: right<29>(msg, idx); idx += 1<<(29-1); break;
			case (2<<5)+30: right<30>(msg, idx); idx += 1<<(30-1); break;
			case (3<<5)+3: idx -= 1<<(3-1); combine<3>(msg, idx); break;
			case (3<<5)+4: idx -= 1<<(4-1); combine<4>(msg, idx); break;
			case (3<<5)+5: idx -= 1<<(5-1); combine<5>(msg, idx); break;
			case (3<<5)+6: idx -= 1<<(6-1); combine<6>(msg, idx); break;
			case (3<<5)+7: idx -= 1<<(7-1); combine<7>(msg, idx); break;
			case (3<<5)+8: idx -= 1<<(8-1); combine<8>(msg, idx); break;
			case (3<<5)+9: idx -= 1<<(9-1); combine<9>(msg, idx); break;
			case (3<<5)+10: idx -= 1<<(10-1); combine<10>(msg, idx); break;
			case (3<<5)+11: idx -= 1<<(11-1); combine<11>(msg, idx); break;
			case (3<<5)+12: idx -= 1<<(12-1); combine<12>(msg, idx); break;
			case (3<<5)+13: idx -= 1<<(13-1); combine<13>(msg, idx); break;
			case (3<<5)+14: idx -= 1<<(14-1); combine<14>(msg, idx); break;
			case (3<<5)+15: idx -= 1<<(15-1); combine<15>(msg, idx); break;
			case (3<<5)+16: idx -= 1<<(16-1); combine<16>(msg, idx); break;
			case (3<<5)+17: idx -= 1<<(17-1); combine<17>(msg, idx); break;
			case (3<<5)+18: idx -= 1<<(18-1); combine<18>(msg, idx); break;
			case (3<<5)+19: idx -= 1<<(19-1); combine<19>(msg, idx); break;
			case (3<<5)+20: idx -= 1<<(20-1); combine<20>(msg, idx); break;
			case (3<<5)+21: idx -= 1<<(21-1); combine<21>(msg, idx); break;
			case (3<<5)+22: idx -= 1<<(22-1); combine<22>(msg, idx); break;
			case (3<<5)+23: idx -= 1<<(23-1); combine<23>(msg, idx); break;
			case (3<<5)+24: idx -= 1<<(24-1); combine<24>(msg, idx); break;
			case (3<<5)+25: idx -= 1<<(25-1); combine<25>(msg, idx); break;
			case (3<<5)+26: idx -= 1<<(26-1); combine<26>(msg, idx); break;
			case (3<<5)+27: idx -= 1<<(27-1); combine<27>(msg, idx); break;
			case (3<<5)+28: idx -= 1<<(28-1); combine<28>(msg, idx); break;
			case (3<<5)+29: idx -= 1<<(29-1); combine<29>(msg, idx); break;
			case (3<<5)+30: idx -= 1<<(30-1); combine<30>(msg, idx); break;
			case (4<<5)+3: rate0<3>(msg, idx); break;
			case (4<<5)+4: rate0<4>(msg, idx); break;
			case (4<<5)+5: rate0<5>(msg, idx); break;
			case (4<<5)+6: rate0<6>(msg, idx); break;
			case (4<<5)+7: rate0<7>(msg, idx); break;
			case (4<<5)+8: rate0<8>(msg, idx); break;
			case (4<<5)+9: rate0<9>(msg, idx); break;
			case (4<<5)+10: rate0<10>(msg, idx); break;
			case (4<<5)+11: rate0<11>(msg, idx); break;
			case (4<<5)+12: rate0<12>(msg, idx); break;
			case (4<<5)+13: rate0<13>(msg, idx); break;
			case (4<<5)+14: rate0<14>(msg, idx); break;
			case (4<<5)+15: rate0<15>(msg, idx); break;
			case (4<<5)+16: rate0<16>(msg, idx); break;
			case (4<<5)+17: rate0<17>(msg, idx); break;
			case (4<<5)+18: rate0<18>(msg, idx); break;
			case (4<<5)+19: rate0<19>(msg, idx); break;
			case (4<<5)+20: rate0<20>(msg, idx); break;
			case (4<<5)+21: rate0<21>(msg, idx); break;
			case (4<<5)+22: rate0<22>(msg, idx); break;
			case (4<<5)+23: rate0<23>(msg, idx); break;
			case (4<<5)+24: rate0<24>(msg, idx); break;
			case (4<<5)+25: rate0<25>(msg, idx); break;
			case (4<<5)+26: rate0<26>(msg, idx); break;
			case (4<<5)+27: rate0<27>(msg, idx); break;
			case (4<<5)+28: rate0<28>(msg, idx); break;
			case (4<<5)+29: rate0<29>(msg, idx); break;
			case (4<<5)+30: rate0<30>(msg, idx); break;
			case (5<<5)+3: rate1<3>(msg, idx); break;
			case (5<<5)+4: rate1<4>(msg, idx); break;
			case (5<<5)+5: rate1<5>(msg, idx); break;
			case (5<<5)+6: rate1<6>(msg, idx); break;
			case (5<<5)+7: rate1<7>(msg, idx); break;
			case (5<<5)+8: rate1<8>(msg, idx); break;
			case (5<<5)+9: rate1<9>(msg, idx); break;
			case (5<<5)+10: rate1<10>(msg, idx); break;
			case (5<<5)+11: rate1<11>(msg, idx); break;
			case (5<<5)+12: rate1<12>(msg, idx); break;
			case (5<<5)+13: rate1<13>(msg, idx); break;
			case (5<<5)+14: rate1<14>(msg, idx); break;
			case (5<<5)+15: rate1<15>(msg, idx); break;
			case (5<<5)+16: rate1<16>(msg, idx); break;
			case (5<<5)+17: rate1<17>(msg, idx); break;
			case (5<<5)+18: rate1<18>(msg, idx); break;
			case (5<<5)+19: rate1<19>(msg, idx); break;
			case (5<<5)+20: rate1<20>(msg, idx); break;
			case (5<<5)+21: rate1<21>(msg, idx); break;
			case (5<<5)+22: rate1<22>(msg, idx); break;
			case (5<<5)+23: rate1<23>(msg, idx); break;
			case (5<<5)+24: rate1<24>(msg, idx); break;
			case (5<<5)+25: rate1<25>(msg, idx); break;
			case (5<<5)+26: rate1<26>(msg, idx); break;
			case (5<<5)+27: rate1<27>(msg, idx); break;
			case (5<<5)+28: rate1<28>(msg, idx); break;
			case (5<<5)+29: rate1<29>(msg, idx); break;
			case (5<<5)+30: rate1<30>(msg, idx); break;
			case (6<<5)+3: rep<3>(msg, idx); break;
			case (6<<5)+4: rep<4>(msg, idx); break;
			case (6<<5)+5: rep<5>(msg, idx); break;
			case (6<<5)+6: rep<6>(msg, idx); break;
			case (6<<5)+7: rep<7>(msg, idx); break;
			case (6<<5)+8: rep<8>(msg, idx); break;
			case (6<<5)+9: rep<9>(msg, idx); break;
			case (6<<5)+10: rep<10>(msg, idx); break;
			case (6<<5)+11: rep<11>(msg, idx); break;
			case (6<<5)+12: rep<12>(msg, idx); break;
			case (6<<5)+13: rep<13>(msg, idx); break;
			case (6<<5)+14: rep<14>(msg, idx); break;
			case (6<<5)+15: rep<15>(msg, idx); break;
			case (6<<5)+16: rep<16>(msg, idx); break;
			case (6<<5)+17: rep<17>(msg, idx); break;
			case (6<<5)+18: rep<18>(msg, idx); break;
			case (6<<5)+19: rep<19>(msg, idx); break;
			case (6<<5)+20: rep<20>(msg, idx); break;
			case (6<<5)+21: rep<21>(msg, idx); break;
			case (6<<5)+22: rep<22>(msg, idx); break;
			case (6<<5)+23: rep<23>(msg, idx); break;
			case (6<<5)+24: rep<24>(msg, idx); break;
			case (6<<5)+25: rep<25>(msg, idx); break;
			case (6<<5)+26: rep<26>(msg, idx); break;
			case (6<<5)+27: rep<27>(msg, idx); break;
			case (6<<5)+28: rep<28>(msg, idx); break;
			case (6<<5)+29: rep<29>(msg, idx); break;
			case (6<<5)+30: rep<30>(msg, idx); break;
			default:
				assert(false);
			}
		}
	}
};

int main()
{
	const int M = 20;
	const int N = 1 << M;
	const int U = 2; // unrolled at level 2
	const bool systematic = true;
	std::random_device rd;
	typedef std::default_random_engine generator;
	typedef std::uniform_int_distribution<int> distribution;
	auto data = std::bind(distribution(0, 1), generator(rd()));
	auto frozen = new uint8_t[N>>U];
	auto codeword = new int8_t[N];
	PolarFreezer freeze;
	double erasure_probability = 0.5;
	double freezing_threshold = 0.5;
	int K = freeze(frozen, M, erasure_probability, freezing_threshold);
	std::cerr << "Polar(" << N << ", " << K << ")" << std::endl;
	auto message = new int8_t[K];
	auto decoded = new int8_t[K];
	PolarEncoder<M> encode;
	auto program = new uint8_t[N];
	PolarCompiler compile;
	int length = compile(program, frozen, M);
	std::cerr << "program length = " << length << std::endl;
	std::cerr << "sizeof(PolarDecoder<M>) = " << sizeof(PolarDecoder<M>) << std::endl;
	auto decode = new PolarDecoder<M>;

#if 0
	auto epos = std::bind(distribution(0, N-1), generator(rd()));
	std::cerr << "errors erasures msec" << std::endl;
	for (int loop = 0; loop < 100; ++loop) {
		for (int i = 0; i < K; ++i)
			message[i] = 1 - 2 * data();
		if (systematic) {
			PolarSysEnc<M> sysenc;
			sysenc(codeword, message, frozen);
			for (int i = 0, j = 0; i < N; ++i)
				if (!(frozen[i>>U]&(1<<(i&((1<<U)-1)))))
					assert(codeword[i] == message[j++]);
		} else {
			encode(codeword, message, frozen);
		}
		for (int i = 0; i < erasure_probability * N; ++i)
			codeword[epos()] = 0;
		auto start = std::chrono::system_clock::now();
		(*decode)(decoded, codeword, program);
		auto end = std::chrono::system_clock::now();
		auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		if (systematic) {
			encode(codeword, decoded, frozen);
			for (int i = 0, j = 0; i < N; ++i)
				if (!(frozen[i>>U]&(1<<(i&((1<<U)-1)))))
					decoded[j++] = codeword[i];
		}
		int erasures = 0;
		for (int i = 0; i < K; ++i)
			erasures += !decoded[i];
		int errors = 0;
		for (int i = 0; i < K; ++i)
			errors += decoded[i] * message[i] < 0;
		std::cout << errors << " " << erasures << " " << msec.count() << std::endl;
	}
#else
	auto orig = new int8_t[N];
	auto noisy = new int8_t[N];
	auto symb = new double[N];
	if (freezing_threshold != 0.5)
		std::cerr << "freezing_threshold not 0.5. SNR design calculations might be wrong." << std::endl;
	double design_SNR = 10 * std::log10(-std::log(erasure_probability));
	std::cerr << "designed for: " << design_SNR << " SNR" << std::endl;
	double low_SNR = std::floor(design_SNR-3);
	double high_SNR = std::ceil(design_SNR+2);
	double min_SNR = high_SNR, max_mbs = 0;
	for (double SNR = low_SNR; SNR <= high_SNR; SNR += 0.1) {
		//double mean_signal = 0;
		double sigma_signal = 1;
		double mean_noise = 0;
		double sigma_noise = std::sqrt(sigma_signal * sigma_signal / (2 * std::pow(10, SNR / 10)));

		typedef std::normal_distribution<double> normal;
		auto awgn = std::bind(normal(mean_noise, sigma_noise), generator(rd()));

		int64_t awgn_errors = 0;
		int64_t quantization_erasures = 0;
		int64_t uncorrected_errors = 0;
		int64_t ambiguity_erasures = 0;
		int64_t loops = 10;
		double avg_mbs = 0;
		for (int l = 0; l < loops; ++l) {
			for (int i = 0; i < K; ++i)
				message[i] = 1 - 2 * data();

			if (systematic) {
				PolarSysEnc<M> sysenc;
				sysenc(codeword, message, frozen);
				for (int i = 0, j = 0; i < N; ++i)
					if (!(frozen[i>>U]&(1<<(i&((1<<U)-1)))))
						assert(codeword[i] == message[j++]);
			} else {
				encode(codeword, message, frozen);
			}

			for (int i = 0; i < N; ++i)
				orig[i] = codeword[i];

			for (int i = 0; i < N; ++i)
				symb[i] = codeword[i];

			for (int i = 0; i < N; ++i)
				symb[i] += awgn();

			// $LLR=log(\frac{p(x=+1|y)}{p(x=-1|y)})$
			// $p(x|\mu,\sigma)=\frac{1}{\sqrt{2\pi}\sigma}}e^{-\frac{(x-\mu)^2}{2\sigma^2}}$
			double DIST = 2; // BPSK
			double fact = DIST / (sigma_noise * sigma_noise);
			for (int i = 0; i < N; ++i)
				codeword[i] = std::min<double>(std::max<double>(std::nearbyint(fact * symb[i]), -128), 127);

			for (int i = 0; i < N; ++i)
				noisy[i] = codeword[i];

			auto start = std::chrono::system_clock::now();
			(*decode)(decoded, codeword, program);
			auto end = std::chrono::system_clock::now();
			auto usec = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
			double mbs = (double)K / usec.count();
			avg_mbs += mbs;

			if (systematic) {
				encode(codeword, decoded, frozen);
				for (int i = 0, j = 0; i < N; ++i)
					if (!(frozen[i>>U]&(1<<(i&((1<<U)-1)))))
						decoded[j++] = codeword[i];
			}

			for (int i = 0; i < N; ++i)
				awgn_errors += noisy[i] * orig[i] < 0;
			for (int i = 0; i < N; ++i)
				quantization_erasures += !noisy[i];
			for (int i = 0; i < K; ++i)
				uncorrected_errors += decoded[i] * message[i] < 0;
			for (int i = 0; i < K; ++i)
				ambiguity_erasures += !decoded[i];
		}

		avg_mbs /= loops;
		max_mbs = std::max(max_mbs, avg_mbs);
		double bit_error_rate = (double)(uncorrected_errors + ambiguity_erasures) / (double)(K * loops);
		if (!uncorrected_errors && !ambiguity_erasures)
			min_SNR = std::min(min_SNR, SNR);

		if (0) {
			std::cerr << SNR << " Es/N0 => AWGN with standard deviation of " << sigma_noise << " and mean " << mean_noise << std::endl;
			std::cerr << awgn_errors << " errors caused by AWGN." << std::endl;
			std::cerr << quantization_erasures << " erasures caused by quantization." << std::endl;
			std::cerr << uncorrected_errors << " errors uncorrected." << std::endl;
			std::cerr << ambiguity_erasures << " ambiguity erasures." << std::endl;
			std::cerr << bit_error_rate << " bit error rate." << std::endl;
			std::cerr << avg_mbs << " megabit per second." << std::endl;
		} else {
			std::cout << SNR << " " << bit_error_rate << " " << avg_mbs << std::endl;
		}
	}
	std::cerr << "QEF at: " << min_SNR << " SNR, speed: " << max_mbs << " Mb/s." << std::endl;
#endif
	return 0;
}
