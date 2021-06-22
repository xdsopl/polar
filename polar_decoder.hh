/*
Successive cancellation list decoding of polar codes

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#pragma once

template <typename TYPE, int M>
struct PolarTree
{
	typedef PolarHelper<TYPE> PH;
	typedef typename PH::PATH PATH;
	static const int N = 1 << M;
	static void decode(PATH *metric, TYPE *message, int *count, TYPE *hard, TYPE *soft, const uint8_t *frozen, int index)
	{
		for (int i = 0; i < N/2; ++i)
			soft[i+N/2] = PH::prod(soft[i+N], soft[i+N/2+N]);
		PolarTree<TYPE, M-1>::decode(metric, message, count, hard, soft, frozen, index);
		for (int i = 0; i < N/2; ++i)
			soft[i+N/2] = PH::madd(hard[index+i], soft[i+N], soft[i+N/2+N]);
		PolarTree<TYPE, M-1>::decode(metric, message, count, hard, soft, frozen+N/2, index+N/2);
		for (int i = 0; i < N/2; ++i)
			hard[index+i] = PH::qmul(hard[index+i], hard[index+i+N/2]);
	}
};

template <typename TYPE>
struct PolarTree<TYPE, 0>
{
	typedef PolarHelper<TYPE> PH;
	typedef typename PH::PATH PATH;
	static void decode(PATH *metric, TYPE *message, int *count, TYPE *hard, TYPE *soft, const uint8_t *frozen, int index)
	{
		(void)metric;
		if (*frozen) {
			hard[index] = PH::one();
		} else {
			hard[index] = PH::signum(soft[1]);
			message[(*count)++] = hard[index];
		}
	}
};

template <typename TYPE, int M>
class PolarDecoder
{
	typedef PolarHelper<TYPE> PH;
	typedef typename TYPE::value_type value_type;
	typedef typename PH::PATH PATH;
	static const int N = 1 << M;
	TYPE soft[2*N];
	TYPE hard[N];
public:
	void operator()(PATH *metric, TYPE *message, const value_type *codeword, const uint8_t *frozen)
	{
		int count = 0;
		metric[0] = 0;
		for (int k = 1; k < TYPE::SIZE; ++k)
			metric[k] = 1000;
		for (int i = 0; i < N; ++i)
			soft[N+i] = vdup<TYPE>(codeword[i]);
		PolarTree<TYPE, M>::decode(metric, message, &count, hard, soft, frozen, 0);
	}
};

