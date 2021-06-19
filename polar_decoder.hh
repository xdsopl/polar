/*
Successive cancellation list decoding of polar codes

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#pragma once

template <typename TYPE, int M>
class PolarTree
{
	typedef PolarHelper<TYPE> PH;
	typedef typename PH::PATH PATH;
	static const int N = 1 << M;
	TYPE temp[N/2];
	PolarTree<TYPE, M-1> left, right;
public:
	void operator()(PATH *metric, TYPE *message, int *count, TYPE *hard, const TYPE *soft, const uint8_t *frozen)
	{
		for (int i = 0; i < N/2; ++i)
			temp[i] = PH::prod(soft[i], soft[i+N/2]);
		left(metric, message, count, hard, temp, frozen);
		for (int i = 0; i < N/2; ++i)
			temp[i] = PH::madd(hard[i], soft[i], soft[i+N/2]);
		right(metric, message, count, hard+N/2, temp, frozen+N/2);
		for (int i = 0; i < N/2; ++i)
			hard[i] = PH::qmul(hard[i], hard[i+N/2]);
	}
};

template <typename TYPE>
class PolarTree<TYPE, 0>
{
	typedef PolarHelper<TYPE> PH;
	typedef typename PH::PATH PATH;
public:
	void operator()(PATH *metric, TYPE *message, int *count, TYPE *hard, const TYPE *soft, const uint8_t *frozen)
	{
		(void)metric;
		if (*frozen) {
			*hard = PH::one();
		} else {
			*hard = PH::signum(*soft);
			message[(*count)++] = *hard;
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
	TYPE hard[N];
	TYPE soft[N];
	PolarTree<TYPE, M> root;
public:
	void operator()(PATH *metric, TYPE *message, const value_type *codeword, const uint8_t *frozen)
	{
		int count = 0;
		metric[0] = 0;
		for (int k = 1; k < TYPE::SIZE; ++k)
			metric[k] = 1000;
		for (int i = 0; i < N; ++i)
			soft[i] = vdup<TYPE>(codeword[i]);
		root(metric, message, &count, hard, soft, frozen);
	}
};

