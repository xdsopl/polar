/*
Polar encoder for non-systematic and systematic codes

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#pragma once

template <typename TYPE, int M>
class PolarEncoder
{
	static const int N = 1 << M;
	static const int U = 2;
	typedef PolarHelper<TYPE> PH;
public:
	void operator()(TYPE *codeword, const TYPE *message, const uint8_t *frozen)
	{
		for (int i = 0; i < N; i += 2) {
			TYPE msg0 = frozen[i>>U]&(1<<(i&((1<<U)-1))) ? PH::one() : *message++;
			TYPE msg1 = frozen[(i+1)>>U]&(1<<((i+1)&((1<<U)-1))) ? PH::one() : *message++;
			codeword[i] = PH::qmul(msg0, msg1);
			codeword[i+1] = msg1;
		}
		for (int h = 2; h < N; h *= 2)
			for (int i = 0; i < N; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					codeword[j] = PH::qmul(codeword[j], codeword[j+h]);
	}
};

template <typename TYPE, int M>
class PolarSysEnc
{
	static const int N = 1 << M;
	static const int U = 2;
	typedef PolarHelper<TYPE> PH;
public:
	void operator()(TYPE *codeword, const TYPE *message, const uint8_t *frozen)
	{
		for (int i = 0; i < N; i += 2) {
			TYPE msg0 = frozen[i>>U]&(1<<(i&((1<<U)-1))) ? PH::one() : *message++;
			TYPE msg1 = frozen[(i+1)>>U]&(1<<((i+1)&((1<<U)-1))) ? PH::one() : *message++;
			codeword[i] = PH::qmul(msg0, msg1);
			codeword[i+1] = msg1;
		}
		for (int h = 2; h < N; h *= 2)
			for (int i = 0; i < N; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					codeword[j] = PH::qmul(codeword[j], codeword[j+h]);
		for (int i = 0; i < N; i += 2) {
			TYPE msg0 = frozen[i>>U]&(1<<(i&((1<<U)-1))) ? PH::one() : codeword[i];
			TYPE msg1 = frozen[(i+1)>>U]&(1<<((i+1)&((1<<U)-1))) ? PH::one() : codeword[i+1];
			codeword[i] = PH::qmul(msg0, msg1);
			codeword[i+1] = msg1;
		}
		for (int h = 2; h < N; h *= 2)
			for (int i = 0; i < N; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					codeword[j] = PH::qmul(codeword[j], codeword[j+h]);
	}
};

