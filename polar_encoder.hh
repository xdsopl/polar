/*
Polar encoder for non-systematic and systematic codes

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#pragma once

template <int M>
class PolarEncoder
{
	static const int N = 1 << M;
public:
	void operator()(int8_t *codeword, const int8_t *message, const uint8_t *frozen)
	{
		for (int i = 0; i < N; i += 2) {
			int8_t msg0 = frozen[i] ? 1 : *message++;
			int8_t msg1 = frozen[i+1] ? 1 : *message++;
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
public:
	void operator()(int8_t *codeword, const int8_t *message, const uint8_t *frozen)
	{
		for (int i = 0; i < N; i += 2) {
			int8_t msg0 = frozen[i] ? 1 : *message++;
			int8_t msg1 = frozen[i+1] ? 1 : *message++;
			codeword[i] = msg0 * msg1;
			codeword[i+1] = msg1;
		}
		for (int h = 2; h < N; h *= 2)
			for (int i = 0; i < N; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					codeword[j] *= codeword[j+h];
		for (int i = 0; i < N; i += 2) {
			int8_t msg0 = frozen[i] ? 1 : codeword[i];
			int8_t msg1 = frozen[i+1] ? 1 : codeword[i+1];
			codeword[i] = msg0 * msg1;
			codeword[i+1] = msg1;
		}
		for (int h = 2; h < N; h *= 2)
			for (int i = 0; i < N; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					codeword[j] *= codeword[j+h];
	}
};

