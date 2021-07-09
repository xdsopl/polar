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
	typedef typename PH::MAP MAP;
	static const int N = 1 << M;
	TYPE temp[N/2];
	PolarTree<TYPE, M-1> left, right;
public:
	MAP operator()(PATH *metric, TYPE *message, MAP *maps, int *count, TYPE *hard, const TYPE *soft, const uint8_t *frozen)
	{
		for (int i = 0; i < N/2; ++i)
			temp[i] = PH::prod(soft[i], soft[i+N/2]);
		MAP lmap = left(metric, message, maps, count, hard, temp, frozen);
		for (int i = 0; i < N/2; ++i)
			temp[i] = PH::madd(hard[i], vshuf(soft[i], lmap), vshuf(soft[i+N/2], lmap));
		MAP rmap = right(metric, message, maps, count, hard+N/2, temp, frozen+N/2);
		for (int i = 0; i < N/2; ++i)
			hard[i] = PH::qmul(vshuf(hard[i], rmap), hard[i+N/2]);
		return vshuf(lmap, rmap);
	}
};

template <typename TYPE>
class PolarTree<TYPE, 0>
{
	typedef PolarHelper<TYPE> PH;
	typedef typename PH::PATH PATH;
	typedef typename PH::MAP MAP;
public:
	MAP operator()(PATH *metric, TYPE *message, MAP *maps, int *count, TYPE *hard, const TYPE *soft, const uint8_t *frozen)
	{
		MAP map;
		TYPE hrd, sft = *soft;
		if (*frozen) {
			for (int k = 0; k < TYPE::SIZE; ++k)
				if (sft.v[k] < 0)
					metric[k] -= sft.v[k];
			hrd = PH::one();
			for (int k = 0; k < TYPE::SIZE; ++k)
				map.v[k] = k;
		} else {
			PATH fork[2*TYPE::SIZE];
			for (int k = 0; k < TYPE::SIZE; ++k)
				fork[k] = fork[k+TYPE::SIZE] = metric[k];
			for (int k = 0; k < TYPE::SIZE; ++k)
				if (soft->v[k] < 0)
					fork[k] -= sft.v[k];
				else
					fork[k+TYPE::SIZE] += sft.v[k];
			int perm[2*TYPE::SIZE];
			for (int k = 0; k < 2*TYPE::SIZE; ++k)
				perm[k] = k;
			std::nth_element(perm, perm+TYPE::SIZE, perm+2*TYPE::SIZE, [fork](int a, int b){ return fork[a] < fork[b]; });
			for (int k = 0; k < TYPE::SIZE; ++k)
				metric[k] = fork[perm[k]];
			for (int k = 0; k < TYPE::SIZE; ++k)
				map.v[k] = perm[k] % TYPE::SIZE;
			for (int k = 0; k < TYPE::SIZE; ++k)
				hrd.v[k] = perm[k] < TYPE::SIZE ? 1 : -1;
			message[*count] = hrd;
			maps[*count] = map;
			++*count;
		}
		*hard = hrd;
		return map;
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
	TYPE hard[N];
	TYPE soft[N];
	MAP maps[N];
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
		root(metric, message, maps, &count, hard, soft, frozen);
		MAP acc = maps[count-1];
		for (int i = count-2; i >= 0; --i) {
			message[i] = vshuf(message[i], acc);
			acc = vshuf(maps[i], acc);
		}
	}
};

