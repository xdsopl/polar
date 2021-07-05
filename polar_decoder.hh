/*
Successive cancellation list decoding of polar codes

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#pragma once

template <typename TYPE, int M0, int M>
struct PolarTree
{
	typedef PolarHelper<TYPE> PH;
	typedef typename PH::PATH PATH;
	typedef typename PH::MAP MAP;
	static const int N = 1 << M;
	static void decode(PATH *metric, TYPE *message, MAP *maps, int *count, TYPE *hard, TYPE *soft, const uint8_t *frozen, int index)
	{
		for (int i = 0; i < N/2; ++i)
			soft[i+N/2] = PH::prod(soft[i+N], soft[i+N/2+N]);
		PolarTree<TYPE, M0, M-1>::decode(metric, message, maps, count, hard, soft, frozen, index);
		for (int i = 0; i < N/2; ++i)
			soft[i+N/2] = PH::madd(hard[index+i], soft[i+N], soft[i+N/2+N]);
		PolarTree<TYPE, M0, M-1>::decode(metric, message, maps, count, hard, soft, frozen+N/2, index+N/2);
		for (int i = 0; i < N/2; ++i)
			hard[index+i] = PH::qmul(hard[index+i], hard[index+i+N/2]);
	}
};

template <typename TYPE, int M0>
struct PolarTree<TYPE, M0, 0>
{
	typedef PolarHelper<TYPE> PH;
	typedef typename PH::PATH PATH;
	typedef typename PH::MAP MAP;
	static void decode(PATH *metric, TYPE *message, MAP *maps, int *count, TYPE *hard, TYPE *soft, const uint8_t *frozen, int index)
	{
		TYPE hrd, sft = soft[1];
		if (*frozen) {
			for (int k = 0; k < TYPE::SIZE; ++k)
				if (sft.v[k] < 0)
					metric[k] -= sft.v[k];
			hrd = PH::one();
		} else {
			PATH fork[2*TYPE::SIZE];
			for (int k = 0; k < TYPE::SIZE; ++k)
				fork[k] = fork[k+TYPE::SIZE] = metric[k];
			for (int k = 0; k < TYPE::SIZE; ++k)
				if (sft.v[k] < 0)
					fork[k] -= sft.v[k];
				else
					fork[k+TYPE::SIZE] += sft.v[k];
			int perm[2*TYPE::SIZE];
			for (int k = 0; k < 2*TYPE::SIZE; ++k)
				perm[k] = k;
			std::nth_element(perm, perm+TYPE::SIZE, perm+2*TYPE::SIZE, [fork](int a, int b){ return fork[a] < fork[b]; });
			for (int k = 0; k < TYPE::SIZE; ++k)
				metric[k] = fork[perm[k]];
			MAP map;
			for (int k = 0; k < TYPE::SIZE; ++k)
				map.v[k] = perm[k] % TYPE::SIZE;
			for (int i = 0; i < (1<<M0); ++i)
				soft[i] = vshuf(soft[i], map);
			for (int i = 0; i < (1<<M0); ++i)
				hard[i] = vshuf(hard[i], map);
			for (int k = 0; k < TYPE::SIZE; ++k)
				hrd.v[k] = perm[k] < TYPE::SIZE ? 1 : -1;
			message[*count] = hrd;
			maps[*count] = map;
			++*count;
		}
		hard[index] = hrd;
	}
};

template <typename TYPE, int M>
class PolarDecoder
{
	typedef PolarHelper<TYPE> PH;
	typedef typename TYPE::value_type value_type;
	typedef typename PH::PATH PATH;
	typedef typename PH::MAP MAP;
	static const int N = 1 << M;
	TYPE soft[2*N];
	TYPE hard[N];
	MAP maps[N];
public:
	void operator()(PATH *metric, TYPE *message, const value_type *codeword, const uint8_t *frozen)
	{
		int count = 0;
		metric[0] = 0;
		for (int k = 1; k < TYPE::SIZE; ++k)
			metric[k] = 1000;
		for (int i = 0; i < N; ++i)
			soft[N+i] = vdup<TYPE>(codeword[i]);
		PolarTree<TYPE, M, M>::decode(metric, message, maps, &count, hard, soft, frozen, 0);
		MAP acc = maps[count-1];
		for (int i = count-2; i >= 0; --i) {
			message[i] = vshuf(message[i], acc);
			MAP map = maps[i];
			maps[i] = acc;
			acc = vshuf(map, acc);
		}
	}
};

