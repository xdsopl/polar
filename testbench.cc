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
class PolarFreezer
{
	static const int N = 1 << M;
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
	int operator()(uint8_t *frozen_bits, double erasure_probability = 0.5, double freezing_threshold = 0.5)
	{
		for (int i = 0; i < N>>U; ++i)
			frozen_bits[i] = 0;
		freeze(frozen_bits, erasure_probability, freezing_threshold, 0, N / 2);
		int K = 0;
		for (int i = 0; i < N; ++i)
			K += !(frozen_bits[i>>U]&(1<<(i&((1<<U)-1))));
		return K;
	}
};

template <int M>
class PolarCompiler
{
	static const int N = 1 << M;
	static const int U = 2;
	static uint32_t leaf(int func, int index)
	{
		return (func << 24) | index;
	}
	static uint32_t node(int func, int level, int index)
	{
		return (func << 29) | (level << 24) | index;
	}
	static bool all_frozen(const uint8_t *frozen, int index, int level)
	{
		if (level > U) {
			if (!all_frozen(frozen, index, level-1))
				return false;
			if (!all_frozen(frozen+(1<<(level-1-U)), index+(1<<(level-1)), level-1))
				return false;
			return true;
		}
		return *frozen == (1<<(1<<U))-1;
	}
	static bool none_frozen(const uint8_t *frozen, int index, int level)
	{
		if (level > U) {
			if (!none_frozen(frozen, index, level-1))
				return false;
			if (!none_frozen(frozen+(1<<(level-1-U)), index+(1<<(level-1)), level-1))
				return false;
			return true;
		}
		return *frozen == 0;
	}
	static void compile(uint32_t **program, const uint8_t *frozen, int index, int level)
	{
		if (level > U) {
			if (all_frozen(frozen, index, level)) {
				*(*program)++ = node(4, level, index);
			} else if (none_frozen(frozen, index, level)) {
				*(*program)++ = node(5, level, index);
			} else {
				*(*program)++ = node(1, level, index);
				compile(program, frozen, index, level-1);
				*(*program)++ = node(2, level, index);
				compile(program, frozen+(1<<(level-1-U)), index+(1<<(level-1)), level-1);
				*(*program)++ = node(3, level, index);
			}
		} else {
			*(*program)++ = leaf(*frozen, index);
		}
	}
public:
	void operator()(uint32_t *program, const uint8_t *frozen)
	{
		uint32_t *first = program;
		compile(&program, frozen, 0, M);
		*program++ = 0xffffffff;
		std::cerr << "program length = " << program - first << std::endl;
	}
};

template <int M>
class PolarDecoder
{
	static const int N = 1 << M;
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
	static void leaf0(int8_t **msg, int8_t *hard, const int8_t *soft)
	{
#if 0
		int8_t soft20 = soft[0];
		int8_t soft21 = soft[1];
		int8_t soft22 = soft[2];
		int8_t soft23 = soft[3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t soft00 = prod(soft10, soft11);
		int8_t hard0 = signum(soft00);
		*(*msg)++ = soft00;
		int8_t soft01 = madd(hard0, soft10, soft11);
		int8_t hard1 = signum(soft01);
		*(*msg)++ = soft01;
		hard0 *= hard1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t soft02 = prod(soft12, soft13);
		int8_t hard2 = signum(soft02);
		*(*msg)++ = soft02;
		int8_t soft03 = madd(hard2, soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = soft03;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[0] = hard0;
		hard[1] = hard1;
		hard[2] = hard2;
		hard[3] = hard3;
#else
		hard[0] = signum(soft[0]);
		hard[1] = signum(soft[1]);
		hard[2] = signum(soft[2]);
		hard[3] = signum(soft[3]);
		*(*msg)++ = hard[0] * hard[1] * hard[2] * hard[3];
		*(*msg)++ = hard[1] * hard[3];
		*(*msg)++ = hard[2] * hard[3];
		*(*msg)++ = hard[3];
#endif
	}
	static void leaf1(int8_t **msg, int8_t *hard, const int8_t *soft)
	{
		int8_t soft20 = soft[0];
		int8_t soft21 = soft[1];
		int8_t soft22 = soft[2];
		int8_t soft23 = soft[3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t hard0 = 1;
		int8_t soft01 = qadd(soft10, soft11);
		int8_t hard1 = signum(soft01);
		*(*msg)++ = soft01;
		hard0 *= hard1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t soft02 = prod(soft12, soft13);
		int8_t hard2 = signum(soft02);
		*(*msg)++ = soft02;
		int8_t soft03 = madd(hard2, soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = soft03;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[0] = hard0;
		hard[1] = hard1;
		hard[2] = hard2;
		hard[3] = hard3;
	}
	static void leaf2(int8_t **msg, int8_t *hard, const int8_t *soft)
	{
		int8_t soft20 = soft[0];
		int8_t soft21 = soft[1];
		int8_t soft22 = soft[2];
		int8_t soft23 = soft[3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t soft00 = prod(soft10, soft11);
		int8_t hard0 = signum(soft00);
		*(*msg)++ = soft00;
		int8_t hard1 = 1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t soft02 = prod(soft12, soft13);
		int8_t hard2 = signum(soft02);
		*(*msg)++ = soft02;
		int8_t soft03 = madd(hard2, soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = soft03;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[0] = hard0;
		hard[1] = hard1;
		hard[2] = hard2;
		hard[3] = hard3;
	}
	static void leaf3(int8_t **msg, int8_t *hard, const int8_t *soft)
	{
		int8_t soft20 = soft[0];
		int8_t soft21 = soft[1];
		int8_t soft22 = soft[2];
		int8_t soft23 = soft[3];
		int8_t hard0 = 1;
		int8_t hard1 = 1;
		int8_t soft12 = qadd(soft20, soft22);
		int8_t soft13 = qadd(soft21, soft23);
		int8_t soft02 = prod(soft12, soft13);
		int8_t hard2 = signum(soft02);
		*(*msg)++ = soft02;
		int8_t soft03 = madd(hard2, soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = soft03;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[0] = hard0;
		hard[1] = hard1;
		hard[2] = hard2;
		hard[3] = hard3;
	}
	static void leaf4(int8_t **msg, int8_t *hard, const int8_t *soft)
	{
		int8_t soft20 = soft[0];
		int8_t soft21 = soft[1];
		int8_t soft22 = soft[2];
		int8_t soft23 = soft[3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t soft00 = prod(soft10, soft11);
		int8_t hard0 = signum(soft00);
		*(*msg)++ = soft00;
		int8_t soft01 = madd(hard0, soft10, soft11);
		int8_t hard1 = signum(soft01);
		*(*msg)++ = soft01;
		hard0 *= hard1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t hard2 = 1;
		int8_t soft03 = qadd(soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = soft03;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[0] = hard0;
		hard[1] = hard1;
		hard[2] = hard2;
		hard[3] = hard3;
	}
	static void leaf5(int8_t **msg, int8_t *hard, const int8_t *soft)
	{
		int8_t soft20 = soft[0];
		int8_t soft21 = soft[1];
		int8_t soft22 = soft[2];
		int8_t soft23 = soft[3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t hard0 = 1;
		int8_t soft01 = qadd(soft10, soft11);
		int8_t hard1 = signum(soft01);
		*(*msg)++ = soft01;
		hard0 *= hard1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t hard2 = 1;
		int8_t soft03 = qadd(soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = soft03;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[0] = hard0;
		hard[1] = hard1;
		hard[2] = hard2;
		hard[3] = hard3;
	}
	static void leaf6(int8_t **msg, int8_t *hard, const int8_t *soft)
	{
		int8_t soft20 = soft[0];
		int8_t soft21 = soft[1];
		int8_t soft22 = soft[2];
		int8_t soft23 = soft[3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t soft00 = prod(soft10, soft11);
		int8_t hard0 = signum(soft00);
		*(*msg)++ = soft00;
		int8_t hard1 = 1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t hard2 = 1;
		int8_t soft03 = qadd(soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = soft03;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[0] = hard0;
		hard[1] = hard1;
		hard[2] = hard2;
		hard[3] = hard3;
	}
	static void leaf7(int8_t **msg, int8_t *hard, const int8_t *soft)
	{
		int8_t soft20 = soft[0];
		int8_t soft21 = soft[1];
		int8_t soft22 = soft[2];
		int8_t soft23 = soft[3];
		int8_t hard0 = 1;
		int8_t hard1 = 1;
		int8_t soft12 = qadd(soft20, soft22);
		int8_t soft13 = qadd(soft21, soft23);
		int8_t hard2 = 1;
		int8_t soft03 = qadd(soft12, soft13);
		int8_t hard3 = signum(soft03);
		*(*msg)++ = soft03;
		hard2 *= hard3;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[0] = hard0;
		hard[1] = hard1;
		hard[2] = hard2;
		hard[3] = hard3;
	}
	static void leaf8(int8_t **msg, int8_t *hard, const int8_t *soft)
	{
		int8_t soft20 = soft[0];
		int8_t soft21 = soft[1];
		int8_t soft22 = soft[2];
		int8_t soft23 = soft[3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t soft00 = prod(soft10, soft11);
		int8_t hard0 = signum(soft00);
		*(*msg)++ = soft00;
		int8_t soft01 = madd(hard0, soft10, soft11);
		int8_t hard1 = signum(soft01);
		*(*msg)++ = soft01;
		hard0 *= hard1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t soft02 = prod(soft12, soft13);
		int8_t hard2 = signum(soft02);
		*(*msg)++ = soft02;
		int8_t hard3 = 1;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[0] = hard0;
		hard[1] = hard1;
		hard[2] = hard2;
		hard[3] = hard3;
	}
	static void leaf9(int8_t **msg, int8_t *hard, const int8_t *soft)
	{
		int8_t soft20 = soft[0];
		int8_t soft21 = soft[1];
		int8_t soft22 = soft[2];
		int8_t soft23 = soft[3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t hard0 = 1;
		int8_t soft01 = qadd(soft10, soft11);
		int8_t hard1 = signum(soft01);
		*(*msg)++ = soft01;
		hard0 *= hard1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t soft02 = prod(soft12, soft13);
		int8_t hard2 = signum(soft02);
		*(*msg)++ = soft02;
		int8_t hard3 = 1;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[0] = hard0;
		hard[1] = hard1;
		hard[2] = hard2;
		hard[3] = hard3;
	}
	static void leaf10(int8_t **msg, int8_t *hard, const int8_t *soft)
	{
		int8_t soft20 = soft[0];
		int8_t soft21 = soft[1];
		int8_t soft22 = soft[2];
		int8_t soft23 = soft[3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t soft00 = prod(soft10, soft11);
		int8_t hard0 = signum(soft00);
		*(*msg)++ = soft00;
		int8_t hard1 = 1;
		int8_t soft12 = madd(hard0, soft20, soft22);
		int8_t soft13 = madd(hard1, soft21, soft23);
		int8_t soft02 = prod(soft12, soft13);
		int8_t hard2 = signum(soft02);
		*(*msg)++ = soft02;
		int8_t hard3 = 1;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[0] = hard0;
		hard[1] = hard1;
		hard[2] = hard2;
		hard[3] = hard3;
	}
	static void leaf11(int8_t **msg, int8_t *hard, const int8_t *soft)
	{
		int8_t soft20 = soft[0];
		int8_t soft21 = soft[1];
		int8_t soft22 = soft[2];
		int8_t soft23 = soft[3];
		int8_t hard0 = 1;
		int8_t hard1 = 1;
		int8_t soft12 = qadd(soft20, soft22);
		int8_t soft13 = qadd(soft21, soft23);
		int8_t soft02 = prod(soft12, soft13);
		int8_t hard2 = signum(soft02);
		*(*msg)++ = soft02;
		int8_t hard3 = 1;
		hard0 *= hard2;
		hard1 *= hard3;
		hard[0] = hard0;
		hard[1] = hard1;
		hard[2] = hard2;
		hard[3] = hard3;
	}
	static void leaf12(int8_t **msg, int8_t *hard, const int8_t *soft)
	{
		int8_t soft20 = soft[0];
		int8_t soft21 = soft[1];
		int8_t soft22 = soft[2];
		int8_t soft23 = soft[3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t soft00 = prod(soft10, soft11);
		int8_t hard0 = signum(soft00);
		*(*msg)++ = soft00;
		int8_t soft01 = madd(hard0, soft10, soft11);
		int8_t hard1 = signum(soft01);
		*(*msg)++ = soft01;
		hard0 *= hard1;
		int8_t hard2 = 1;
		int8_t hard3 = 1;
		hard[0] = hard0;
		hard[1] = hard1;
		hard[2] = hard2;
		hard[3] = hard3;
	}
	static void leaf13(int8_t **msg, int8_t *hard, const int8_t *soft)
	{
		int8_t soft20 = soft[0];
		int8_t soft21 = soft[1];
		int8_t soft22 = soft[2];
		int8_t soft23 = soft[3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t hard0 = 1;
		int8_t soft01 = qadd(soft10, soft11);
		int8_t hard1 = signum(soft01);
		*(*msg)++ = soft01;
		hard0 *= hard1;
		int8_t hard2 = 1;
		int8_t hard3 = 1;
		hard[0] = hard0;
		hard[1] = hard1;
		hard[2] = hard2;
		hard[3] = hard3;
	}
	static void leaf14(int8_t **msg, int8_t *hard, const int8_t *soft)
	{
		int8_t soft20 = soft[0];
		int8_t soft21 = soft[1];
		int8_t soft22 = soft[2];
		int8_t soft23 = soft[3];
		int8_t soft10 = prod(soft20, soft22);
		int8_t soft11 = prod(soft21, soft23);
		int8_t soft00 = prod(soft10, soft11);
		int8_t hard0 = signum(soft00);
		*(*msg)++ = soft00;
		int8_t hard1 = 1;
		int8_t hard2 = 1;
		int8_t hard3 = 1;
		hard[0] = hard0;
		hard[1] = hard1;
		hard[2] = hard2;
		hard[3] = hard3;
	}
	static void leaf15(int8_t **, int8_t *hard, const int8_t *)
	{
		int8_t hard0 = 1;
		int8_t hard1 = 1;
		int8_t hard2 = 1;
		int8_t hard3 = 1;
		hard[0] = hard0;
		hard[1] = hard1;
		hard[2] = hard2;
		hard[3] = hard3;
	}
	template <int level>
	void node1(int8_t **, int index)
	{
		assert(level <= M);
		int length = 1 << level;
		for (int i = index; i < index+length/2; ++i)
			soft[level-1-U][i] = prod(soft[level-U][i], soft[level-U][i+length/2]);
	}
	template <int level>
	void node2(int8_t **, int index)
	{
		assert(level <= M);
		int length = 1 << level;
		for (int i = index; i < index+length/2; ++i)
			soft[level-1-U][i+length/2] = madd(hard[i], soft[level-U][i], soft[level-U][i+length/2]);
	}
	template <int level>
	void node3(int8_t **, int index)
	{
		assert(level <= M);
		int length = 1 << level;
		for (int i = index; i < index+length/2; ++i)
			hard[i] *= hard[i+length/2];
	}
	template <int level>
	void rate0(int8_t **, int index)
	{
		assert(level <= M);
		int length = 1 << level;
		for (int i = index; i < index+length; ++i)
			hard[i] = 1;
	}
	template <int level>
	void rate1(int8_t **msg, int index)
	{
		assert(level <= M);
		int length = 1 << level;
		for (int i = index; i < index+length; ++i)
			hard[i] = signum(soft[level-U][i]);
		for (int i = 0; i < length; i += 2) {
			(*msg)[i] = hard[index+i] * hard[index+i+1];
			(*msg)[i+1] = hard[index+i+1];
		}
		for (int h = 2; h < length; h *= 2)
			for (int i = 0; i < length; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					(*msg)[j] *= (*msg)[j+h];
		*msg += length;
	}
	int8_t soft[M+1-U][N];
	int8_t hard[N];
	void decode(int8_t **msg, int func, int index)
	{
		switch (func) {
		case 0: leaf0(msg, hard+index, soft[0]+index); break;
		case 1: leaf1(msg, hard+index, soft[0]+index); break;
		case 2: leaf2(msg, hard+index, soft[0]+index); break;
		case 3: leaf3(msg, hard+index, soft[0]+index); break;
		case 4: leaf4(msg, hard+index, soft[0]+index); break;
		case 5: leaf5(msg, hard+index, soft[0]+index); break;
		case 6: leaf6(msg, hard+index, soft[0]+index); break;
		case 7: leaf7(msg, hard+index, soft[0]+index); break;
		case 8: leaf8(msg, hard+index, soft[0]+index); break;
		case 9: leaf9(msg, hard+index, soft[0]+index); break;
		case 10: leaf10(msg, hard+index, soft[0]+index); break;
		case 11: leaf11(msg, hard+index, soft[0]+index); break;
		case 12: leaf12(msg, hard+index, soft[0]+index); break;
		case 13: leaf13(msg, hard+index, soft[0]+index); break;
		case 14: leaf14(msg, hard+index, soft[0]+index); break;
		case 15: leaf15(msg, hard+index, soft[0]+index); break;
		case 35: node1<3>(msg, index); break;
		case 36: node1<4>(msg, index); break;
		case 37: node1<5>(msg, index); break;
		case 38: node1<6>(msg, index); break;
		case 39: node1<7>(msg, index); break;
		case 40: node1<8>(msg, index); break;
		case 41: node1<9>(msg, index); break;
		case 42: node1<10>(msg, index); break;
		case 43: node1<11>(msg, index); break;
		case 44: node1<12>(msg, index); break;
		case 45: node1<13>(msg, index); break;
		case 46: node1<14>(msg, index); break;
		case 47: node1<15>(msg, index); break;
		case 48: node1<16>(msg, index); break;
		case 49: node1<17>(msg, index); break;
		case 50: node1<18>(msg, index); break;
		case 51: node1<19>(msg, index); break;
		case 52: node1<20>(msg, index); break;
		case 53: node1<21>(msg, index); break;
		case 54: node1<22>(msg, index); break;
		case 55: node1<23>(msg, index); break;
		case 56: node1<24>(msg, index); break;
		case 57: node1<25>(msg, index); break;
		case 58: node1<26>(msg, index); break;
		case 59: node1<27>(msg, index); break;
		case 60: node1<28>(msg, index); break;
		case 61: node1<29>(msg, index); break;
		case 62: node1<30>(msg, index); break;
		case 63: node1<31>(msg, index); break;
		case 67: node2<3>(msg, index); break;
		case 68: node2<4>(msg, index); break;
		case 69: node2<5>(msg, index); break;
		case 70: node2<6>(msg, index); break;
		case 71: node2<7>(msg, index); break;
		case 72: node2<8>(msg, index); break;
		case 73: node2<9>(msg, index); break;
		case 74: node2<10>(msg, index); break;
		case 75: node2<11>(msg, index); break;
		case 76: node2<12>(msg, index); break;
		case 77: node2<13>(msg, index); break;
		case 78: node2<14>(msg, index); break;
		case 79: node2<15>(msg, index); break;
		case 80: node2<16>(msg, index); break;
		case 81: node2<17>(msg, index); break;
		case 82: node2<18>(msg, index); break;
		case 83: node2<19>(msg, index); break;
		case 84: node2<20>(msg, index); break;
		case 85: node2<21>(msg, index); break;
		case 86: node2<22>(msg, index); break;
		case 87: node2<23>(msg, index); break;
		case 88: node2<24>(msg, index); break;
		case 89: node2<25>(msg, index); break;
		case 90: node2<26>(msg, index); break;
		case 91: node2<27>(msg, index); break;
		case 92: node2<28>(msg, index); break;
		case 93: node2<29>(msg, index); break;
		case 94: node2<30>(msg, index); break;
		case 95: node2<31>(msg, index); break;
		case 99: node3<3>(msg, index); break;
		case 100: node3<4>(msg, index); break;
		case 101: node3<5>(msg, index); break;
		case 102: node3<6>(msg, index); break;
		case 103: node3<7>(msg, index); break;
		case 104: node3<8>(msg, index); break;
		case 105: node3<9>(msg, index); break;
		case 106: node3<10>(msg, index); break;
		case 107: node3<11>(msg, index); break;
		case 108: node3<12>(msg, index); break;
		case 109: node3<13>(msg, index); break;
		case 110: node3<14>(msg, index); break;
		case 111: node3<15>(msg, index); break;
		case 112: node3<16>(msg, index); break;
		case 113: node3<17>(msg, index); break;
		case 114: node3<18>(msg, index); break;
		case 115: node3<19>(msg, index); break;
		case 116: node3<20>(msg, index); break;
		case 117: node3<21>(msg, index); break;
		case 118: node3<22>(msg, index); break;
		case 119: node3<23>(msg, index); break;
		case 120: node3<24>(msg, index); break;
		case 121: node3<25>(msg, index); break;
		case 122: node3<26>(msg, index); break;
		case 123: node3<27>(msg, index); break;
		case 124: node3<28>(msg, index); break;
		case 125: node3<29>(msg, index); break;
		case 126: node3<30>(msg, index); break;
		case 127: node3<31>(msg, index); break;
		case 131: rate0<3>(msg, index); break;
		case 132: rate0<4>(msg, index); break;
		case 133: rate0<5>(msg, index); break;
		case 134: rate0<6>(msg, index); break;
		case 135: rate0<7>(msg, index); break;
		case 136: rate0<8>(msg, index); break;
		case 137: rate0<9>(msg, index); break;
		case 138: rate0<10>(msg, index); break;
		case 139: rate0<11>(msg, index); break;
		case 140: rate0<12>(msg, index); break;
		case 141: rate0<13>(msg, index); break;
		case 142: rate0<14>(msg, index); break;
		case 143: rate0<15>(msg, index); break;
		case 144: rate0<16>(msg, index); break;
		case 145: rate0<17>(msg, index); break;
		case 146: rate0<18>(msg, index); break;
		case 147: rate0<19>(msg, index); break;
		case 148: rate0<20>(msg, index); break;
		case 149: rate0<21>(msg, index); break;
		case 150: rate0<22>(msg, index); break;
		case 151: rate0<23>(msg, index); break;
		case 152: rate0<24>(msg, index); break;
		case 153: rate0<25>(msg, index); break;
		case 154: rate0<26>(msg, index); break;
		case 155: rate0<27>(msg, index); break;
		case 156: rate0<28>(msg, index); break;
		case 157: rate0<29>(msg, index); break;
		case 158: rate0<30>(msg, index); break;
		case 159: rate0<31>(msg, index); break;
		case 163: rate1<3>(msg, index); break;
		case 164: rate1<4>(msg, index); break;
		case 165: rate1<5>(msg, index); break;
		case 166: rate1<6>(msg, index); break;
		case 167: rate1<7>(msg, index); break;
		case 168: rate1<8>(msg, index); break;
		case 169: rate1<9>(msg, index); break;
		case 170: rate1<10>(msg, index); break;
		case 171: rate1<11>(msg, index); break;
		case 172: rate1<12>(msg, index); break;
		case 173: rate1<13>(msg, index); break;
		case 174: rate1<14>(msg, index); break;
		case 175: rate1<15>(msg, index); break;
		case 176: rate1<16>(msg, index); break;
		case 177: rate1<17>(msg, index); break;
		case 178: rate1<18>(msg, index); break;
		case 179: rate1<19>(msg, index); break;
		case 180: rate1<20>(msg, index); break;
		case 181: rate1<21>(msg, index); break;
		case 182: rate1<22>(msg, index); break;
		case 183: rate1<23>(msg, index); break;
		case 184: rate1<24>(msg, index); break;
		case 185: rate1<25>(msg, index); break;
		case 186: rate1<26>(msg, index); break;
		case 187: rate1<27>(msg, index); break;
		case 188: rate1<28>(msg, index); break;
		case 189: rate1<29>(msg, index); break;
		case 190: rate1<30>(msg, index); break;
		case 191: rate1<31>(msg, index); break;
		default:
			assert(false);
		}
	}
public:
	void operator()(int8_t *message, const int8_t *codeword, const uint32_t *program)
	{
		for (int i = 0; i < N; ++i)
			soft[M-U][i] = codeword[i];
		while (*program != 0xffffffff) {
			int func = *program >> 24;
			int index = *program & 0x00ffffff;
			decode(&message, func, index);
			++program;
		};
	}
};

int main()
{
	const int M = 20;
	const int N = 1 << M;
	const int U = 2; // unrolled at level 2
	std::random_device rd;
	typedef std::default_random_engine generator;
	typedef std::uniform_int_distribution<int> distribution;
	auto data = std::bind(distribution(0, 1), generator(rd()));
	auto epos = std::bind(distribution(0, N-1), generator(rd()));
	auto frozen = new uint8_t[N>>U];
	auto codeword = new int8_t[N];
	PolarFreezer<M> freeze;
	double erasure_probability = 0.5;
	double freezing_threshold = 0.5;
	int K = freeze(frozen, erasure_probability, freezing_threshold);
	std::cerr << "Polar(" << N << ", " << K << ")" << std::endl;
	auto message = new int8_t[K];
	auto decoded = new int8_t[K];
	PolarEncoder<M> encode;
	auto program = new uint32_t[N];
	PolarCompiler<M> compile;
	compile(program, frozen);
	std::cerr << "sizeof(PolarDecoder<M>) = " << sizeof(PolarDecoder<M>) << std::endl;
	auto decode = new PolarDecoder<M>;

	std::cerr << "errors erasures msec" << std::endl;
	for (int loop = 0; loop < 100; ++loop) {
		for (int i = 0; i < K; ++i)
			message[i] = 1 - 2 * data();
		encode(codeword, message, frozen);
		for (int i = 0; i < erasure_probability * N; ++i)
			codeword[epos()] = 0;
		auto start = std::chrono::system_clock::now();
		(*decode)(decoded, codeword, program);
		auto end = std::chrono::system_clock::now();
		auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		int erasures = 0;
		for (int i = 0; i < K; ++i)
			erasures += !decoded[i];
		int errors = 0;
		auto signum = [](int v){ return (v > 0) - (v < 0); };
		for (int i = 0; i < K; ++i)
			errors += signum(decoded[i]) * message[i] < 0;
		std::cout << errors << " " << erasures << " " << msec.count() << std::endl;
	}

	return 0;
}
