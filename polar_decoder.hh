/*
Successive cancellation decoding of polar codes

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#pragma once

template <int MAX_M>
class PolarDecoder
{
#ifdef __AVX2__
	static const int LEVEL = 5;
#else
	static const int LEVEL = 4;
#endif
	static const int WIDTH = 1 << LEVEL;
	static const int COUNT = 1 << (MAX_M - LEVEL);
	typedef SIMD<int8_t, WIDTH> VECT;

	static int8_t signum(int8_t v)
	{
		return (v > 0) - (v < 0);
	}
	static int8_t decide(int8_t v)
	{
		return (v >= 0) - (v < 0);
	}
	static int8_t qabs(int8_t a)
	{
		return std::abs(std::max<int8_t>(a, -127));
	}
	static int8_t qmin(int8_t a, int8_t b)
	{
		return std::min(a, b);
	}
	static int8_t qadd(int8_t a, int8_t b)
	{
		return std::min<int16_t>(std::max<int16_t>(int16_t(a) + int16_t(b), -128), 127);
	}
	static int8_t qmul(int8_t a, int8_t b)
	{
		// return std::min<int16_t>(std::max<int16_t>(int16_t(a) * int16_t(b), -128), 127);
		// only used for hard decision values anyway
		return a * b;
	}
	static int8_t prod(int8_t a, int8_t b)
	{
		return signum(a) * signum(b) * qmin(qabs(a), qabs(b));
	}
	static int8_t madd(int8_t a, int8_t b, int8_t c)
	{
		return std::min<int16_t>(std::max<int16_t>(int16_t(a) * int16_t(b) + int16_t(c), -128), 127);
	}
	static VECT one()
	{
		return vdup<VECT>(1);
	}
	static VECT zero()
	{
		return vzero<VECT>();
	}
	static VECT signum(VECT a)
	{
		return vsignum(a);
	}
	static VECT decide(VECT a)
	{
		return vreinterpret<VECT>(vorr(vmask(one()), vcltz(a)));
	}
	static VECT qabs(VECT a)
	{
		return vqabs(a);
	}
	static VECT qmin(VECT a, VECT b)
	{
		return vmin(a, b);
	}
	static VECT qadd(VECT a, VECT b)
	{
		return vqadd(a, b);
	}
	static VECT qmul(VECT a, VECT b)
	{
#ifdef __ARM_NEON__
		return vmul(a, b);
#else
		return vsign(a, b);
#endif
	}
	static VECT prod(VECT a, VECT b)
	{
#ifdef __ARM_NEON__
		return vmul(vmul(vsignum(a), vsignum(b)), vmin(vqabs(a), vqabs(b)));
#else
		return vsign(vmin(vqabs(a), vqabs(b)), vsign(vsignum(a), b));
#endif
	}
	static VECT madd(VECT a, VECT b, VECT c)
	{
#ifdef __ARM_NEON__
		return vqadd(vmul(a, vmax(b, vdup<VECT>(-127))), c);
#else
		return vqadd(vsign(vmax(b, vdup<VECT>(-127)), a), c);
#endif
	}
	template <int level>
	void left(int8_t *, int)
	{
		assert(level <= MAX_M);
		if (level <= LEVEL) {
			int length = 1 << level;
			for (int i = 0; i < length/2; ++i)
				reinterpret_cast<int8_t *>(soft)[i+length/2] =
					prod(reinterpret_cast<int8_t *>(soft)[i+length],
					reinterpret_cast<int8_t *>(soft)[i+length/2+length]);
		} else {
			int count = 1 << (level - LEVEL);
			for (int i = 0; i < count/2; ++i)
				soft[i+count/2] = prod(soft[i+count], soft[i+count/2+count]);
		}
	}
	template <int level>
	void right(int8_t *, int index)
	{
		assert(level <= MAX_M);
		if (level <= LEVEL) {
			int length = 1 << level;
			for (int i = 0; i < length/2; ++i)
				reinterpret_cast<int8_t *>(soft)[i+length/2] =
					madd(reinterpret_cast<int8_t *>(hard)[index+i],
					reinterpret_cast<int8_t *>(soft)[i+length],
					reinterpret_cast<int8_t *>(soft)[i+length/2+length]);
		} else {
			int count = 1 << (level - LEVEL);
			for (int i = 0; i < count/2; ++i)
				soft[i+count/2] = madd(hard[index/WIDTH+i], soft[i+count], soft[i+count/2+count]);
		}
	}
	template <int level>
	void comb(int8_t *, int index)
	{
		assert(level <= MAX_M);
		if (level <= LEVEL) {
			int length = 1 << level;
			for (int i = 0; i < length/2; ++i)
				reinterpret_cast<int8_t *>(hard)[index+i] =
					qmul(reinterpret_cast<int8_t *>(hard)[index+i],
					reinterpret_cast<int8_t *>(hard)[index+i+length/2]);
		} else {
			int count = 1 << (level - LEVEL);
			for (int i = 0; i < count/2; ++i)
				hard[index/WIDTH+i] = qmul(hard[index/WIDTH+i], hard[index/WIDTH+i+count/2]);
		}
	}
	template <int level>
	void rate0(int8_t *, int index)
	{
		assert(level <= MAX_M);
		if (level < LEVEL) {
			int length = 1 << level;
			for (int i = 0; i < length; ++i)
				reinterpret_cast<int8_t *>(hard)[index+i] = 1;
		} else {
			int count = 1 << (level - LEVEL);
			for (int i = 0; i < count; ++i)
				hard[index/WIDTH+i] = one();
		}
	}
	template <int level>
	void rate1(int8_t *mesg, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		if (level < LEVEL) {
			for (int i = 0; i < length; ++i)
				reinterpret_cast<int8_t *>(hard)[index+i] = signum(reinterpret_cast<int8_t *>(soft)[i+length]);
		} else {
			int count = 1 << (level - LEVEL);
			for (int i = 0; i < count; ++i)
				hard[index/WIDTH+i] = signum(soft[i+count]);
		}
		for (int i = 0; i < length; i += 2) {
			reinterpret_cast<int8_t *>(soft)[i] = qmul(reinterpret_cast<int8_t *>(hard)[index+i], reinterpret_cast<int8_t *>(hard)[index+i+1]);
			reinterpret_cast<int8_t *>(soft)[i+1] = reinterpret_cast<int8_t *>(hard)[index+i+1];
		}
		for (int h = 2; h < length; h *= 2)
			for (int i = 0; i < length; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					reinterpret_cast<int8_t *>(soft)[j] = qmul(reinterpret_cast<int8_t *>(soft)[j], reinterpret_cast<int8_t *>(soft)[j+h]);
		for (int i = 0; i < length; ++i)
			mesg[i] = reinterpret_cast<int8_t *>(soft)[i];
	}
	template <int level>
	void rep(int8_t *mesg, int index)
	{
		assert(level <= MAX_M);
		if (level < LEVEL) {
			int length = 1 << level;
			for (int h = length; h; h /= 2)
				for (int i = 0; i < h/2; ++i)
					reinterpret_cast<int8_t *>(soft)[i+h/2] =
						qadd(reinterpret_cast<int8_t *>(soft)[i+h],
						reinterpret_cast<int8_t *>(soft)[i+h/2+h]);
			int8_t hardi = signum(reinterpret_cast<int8_t *>(soft)[1]);
			*mesg = hardi;
			for (int i = 0; i < length; ++i)
				reinterpret_cast<int8_t *>(hard)[index+i] = hardi;
		} else {
			int count = 1 << (level - LEVEL);
			VECT vsum = soft[count];
			for (int i = 1; i < count; ++i)
				vsum = qadd(vsum, soft[i+count]);
			int8_t hsum = vsum.v[0];
			for (int i = 1; i < WIDTH; ++i)
				hsum = qadd(hsum, vsum.v[i]);
			int8_t hardi = signum(hsum);
			*mesg = hardi;
			for (int i = 0; i < count; ++i)
				hard[index/WIDTH+i] = vdup<VECT>(hardi);
		}
	}
	template <int level>
	void spc(int8_t *mesg, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		int8_t parity;
		if (level < LEVEL) {
			for (int i = 0; i < length; ++i)
				reinterpret_cast<int8_t *>(hard)[index+i] = decide(reinterpret_cast<int8_t *>(soft)[i+length]);
			parity = reinterpret_cast<int8_t *>(hard)[index];
			for (int i = 1; i < length; ++i)
				parity = qmul(parity, reinterpret_cast<int8_t *>(hard)[index+i]);
		} else {
			int count = 1 << (level - LEVEL);
			for (int i = 0; i < count; ++i)
				hard[index/WIDTH+i] = signum(soft[i+count]);
			VECT vpar = hard[index/WIDTH];
			for (int i = 1; i < count; ++i)
				vpar = qmul(vpar, hard[index/WIDTH+i]);
			parity = vpar.v[0];
			for (int i = 1; i < WIDTH; ++i)
				parity = qmul(parity, vpar.v[i]);
		}
		int worst = 0;
		for (int i = 1; i < length; ++i)
			if (qabs(reinterpret_cast<int8_t *>(soft)[worst+length]) > qabs(reinterpret_cast<int8_t *>(soft)[i+length]))
				worst = i;
		reinterpret_cast<int8_t *>(hard)[index+worst] = qmul(reinterpret_cast<int8_t *>(hard)[index+worst], parity);
		for (int i = 0; i < length; i += 2) {
			reinterpret_cast<int8_t *>(soft)[i] = qmul(reinterpret_cast<int8_t *>(hard)[index+i], reinterpret_cast<int8_t *>(hard)[index+i+1]);
			reinterpret_cast<int8_t *>(soft)[i+1] = reinterpret_cast<int8_t *>(hard)[index+i+1];
		}
		for (int h = 2; h < length; h *= 2)
			for (int i = 0; i < length; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					reinterpret_cast<int8_t *>(soft)[j] = qmul(reinterpret_cast<int8_t *>(soft)[j], reinterpret_cast<int8_t *>(soft)[j+h]);
		for (int i = 0; i < length-1; ++i)
			mesg[i] = reinterpret_cast<int8_t *>(soft)[i+1];
	}
	VECT soft[2*COUNT];
	VECT hard[COUNT];
public:
	void operator()(int8_t *message, const int8_t *codeword, const uint8_t *program)
	{
		int level = *program++;
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length; ++i)
			reinterpret_cast<int8_t *>(soft)[i+length] = codeword[i];
		int idx = 0, lvl = level;
		int8_t *msg = message;
		while (*program != 255) {
			switch (*program++) {
			case 0:	switch (lvl--) {
				case 2: left<2>(msg, idx); break;
				case 3: left<3>(msg, idx); break;
				case 4: left<4>(msg, idx); break;
				case 5: left<5>(msg, idx); break;
				case 6: left<6>(msg, idx); break;
				case 7: left<7>(msg, idx); break;
				case 8: left<8>(msg, idx); break;
				case 9: left<9>(msg, idx); break;
				case 10: left<10>(msg, idx); break;
				case 11: left<11>(msg, idx); break;
				case 12: left<12>(msg, idx); break;
				case 13: left<13>(msg, idx); break;
				case 14: left<14>(msg, idx); break;
				case 15: left<15>(msg, idx); break;
				case 16: left<16>(msg, idx); break;
				case 17: left<17>(msg, idx); break;
				case 18: left<18>(msg, idx); break;
				case 19: left<19>(msg, idx); break;
				case 20: left<20>(msg, idx); break;
				case 21: left<21>(msg, idx); break;
				case 22: left<22>(msg, idx); break;
				case 23: left<23>(msg, idx); break;
				case 24: left<24>(msg, idx); break;
				case 25: left<25>(msg, idx); break;
				case 26: left<26>(msg, idx); break;
				case 27: left<27>(msg, idx); break;
				case 28: left<28>(msg, idx); break;
				case 29: left<29>(msg, idx); break;
				case 30: left<30>(msg, idx); break;
				case 31: left<31>(msg, idx); break;
				default: assert(false);
				} break;
			case 1: switch (lvl+1) {
				case 2: right<2>(msg, idx); break;
				case 3: right<3>(msg, idx); break;
				case 4: right<4>(msg, idx); break;
				case 5: right<5>(msg, idx); break;
				case 6: right<6>(msg, idx); break;
				case 7: right<7>(msg, idx); break;
				case 8: right<8>(msg, idx); break;
				case 9: right<9>(msg, idx); break;
				case 10: right<10>(msg, idx); break;
				case 11: right<11>(msg, idx); break;
				case 12: right<12>(msg, idx); break;
				case 13: right<13>(msg, idx); break;
				case 14: right<14>(msg, idx); break;
				case 15: right<15>(msg, idx); break;
				case 16: right<16>(msg, idx); break;
				case 17: right<17>(msg, idx); break;
				case 18: right<18>(msg, idx); break;
				case 19: right<19>(msg, idx); break;
				case 20: right<20>(msg, idx); break;
				case 21: right<21>(msg, idx); break;
				case 22: right<22>(msg, idx); break;
				case 23: right<23>(msg, idx); break;
				case 24: right<24>(msg, idx); break;
				case 25: right<25>(msg, idx); break;
				case 26: right<26>(msg, idx); break;
				case 27: right<27>(msg, idx); break;
				case 28: right<28>(msg, idx); break;
				case 29: right<29>(msg, idx); break;
				case 30: right<30>(msg, idx); break;
				case 31: right<31>(msg, idx); break;
				default: assert(false);
				} idx += 1<<lvl; break;
			case 2: idx -= 1<<lvl; switch (++lvl) {
				case 2: comb<2>(msg, idx); break;
				case 3: comb<3>(msg, idx); break;
				case 4: comb<4>(msg, idx); break;
				case 5: comb<5>(msg, idx); break;
				case 6: comb<6>(msg, idx); break;
				case 7: comb<7>(msg, idx); break;
				case 8: comb<8>(msg, idx); break;
				case 9: comb<9>(msg, idx); break;
				case 10: comb<10>(msg, idx); break;
				case 11: comb<11>(msg, idx); break;
				case 12: comb<12>(msg, idx); break;
				case 13: comb<13>(msg, idx); break;
				case 14: comb<14>(msg, idx); break;
				case 15: comb<15>(msg, idx); break;
				case 16: comb<16>(msg, idx); break;
				case 17: comb<17>(msg, idx); break;
				case 18: comb<18>(msg, idx); break;
				case 19: comb<19>(msg, idx); break;
				case 20: comb<20>(msg, idx); break;
				case 21: comb<21>(msg, idx); break;
				case 22: comb<22>(msg, idx); break;
				case 23: comb<23>(msg, idx); break;
				case 24: comb<24>(msg, idx); break;
				case 25: comb<25>(msg, idx); break;
				case 26: comb<26>(msg, idx); break;
				case 27: comb<27>(msg, idx); break;
				case 28: comb<28>(msg, idx); break;
				case 29: comb<29>(msg, idx); break;
				case 30: comb<30>(msg, idx); break;
				case 31: comb<31>(msg, idx); break;
				default: assert(false);
				} break;
			case 3: switch (lvl) {
				case 1: rate0<1>(msg, idx); break;
				case 2: rate0<2>(msg, idx); break;
				case 3: rate0<3>(msg, idx); break;
				case 4: rate0<4>(msg, idx); break;
				case 5: rate0<5>(msg, idx); break;
				case 6: rate0<6>(msg, idx); break;
				case 7: rate0<7>(msg, idx); break;
				case 8: rate0<8>(msg, idx); break;
				case 9: rate0<9>(msg, idx); break;
				case 10: rate0<10>(msg, idx); break;
				case 11: rate0<11>(msg, idx); break;
				case 12: rate0<12>(msg, idx); break;
				case 13: rate0<13>(msg, idx); break;
				case 14: rate0<14>(msg, idx); break;
				case 15: rate0<15>(msg, idx); break;
				case 16: rate0<16>(msg, idx); break;
				case 17: rate0<17>(msg, idx); break;
				case 18: rate0<18>(msg, idx); break;
				case 19: rate0<19>(msg, idx); break;
				case 20: rate0<20>(msg, idx); break;
				case 21: rate0<21>(msg, idx); break;
				case 22: rate0<22>(msg, idx); break;
				case 23: rate0<23>(msg, idx); break;
				case 24: rate0<24>(msg, idx); break;
				case 25: rate0<25>(msg, idx); break;
				case 26: rate0<26>(msg, idx); break;
				case 27: rate0<27>(msg, idx); break;
				case 28: rate0<28>(msg, idx); break;
				case 29: rate0<29>(msg, idx); break;
				case 30: rate0<30>(msg, idx); break;
				default: assert(false);
				} break;
			case 4: switch (lvl) {
				case 1: rate1<1>(msg, idx); break;
				case 2: rate1<2>(msg, idx); break;
				case 3: rate1<3>(msg, idx); break;
				case 4: rate1<4>(msg, idx); break;
				case 5: rate1<5>(msg, idx); break;
				case 6: rate1<6>(msg, idx); break;
				case 7: rate1<7>(msg, idx); break;
				case 8: rate1<8>(msg, idx); break;
				case 9: rate1<9>(msg, idx); break;
				case 10: rate1<10>(msg, idx); break;
				case 11: rate1<11>(msg, idx); break;
				case 12: rate1<12>(msg, idx); break;
				case 13: rate1<13>(msg, idx); break;
				case 14: rate1<14>(msg, idx); break;
				case 15: rate1<15>(msg, idx); break;
				case 16: rate1<16>(msg, idx); break;
				case 17: rate1<17>(msg, idx); break;
				case 18: rate1<18>(msg, idx); break;
				case 19: rate1<19>(msg, idx); break;
				case 20: rate1<20>(msg, idx); break;
				case 21: rate1<21>(msg, idx); break;
				case 22: rate1<22>(msg, idx); break;
				case 23: rate1<23>(msg, idx); break;
				case 24: rate1<24>(msg, idx); break;
				case 25: rate1<25>(msg, idx); break;
				case 26: rate1<26>(msg, idx); break;
				case 27: rate1<27>(msg, idx); break;
				case 28: rate1<28>(msg, idx); break;
				case 29: rate1<29>(msg, idx); break;
				case 30: rate1<30>(msg, idx); break;
				default: assert(false);
				} msg += 1<<lvl; break;
			case 5: switch (lvl) {
				case 1: rep<1>(msg, idx); break;
				case 2: rep<2>(msg, idx); break;
				case 3: rep<3>(msg, idx); break;
				case 4: rep<4>(msg, idx); break;
				case 5: rep<5>(msg, idx); break;
				case 6: rep<6>(msg, idx); break;
				case 7: rep<7>(msg, idx); break;
				case 8: rep<8>(msg, idx); break;
				case 9: rep<9>(msg, idx); break;
				case 10: rep<10>(msg, idx); break;
				case 11: rep<11>(msg, idx); break;
				case 12: rep<12>(msg, idx); break;
				case 13: rep<13>(msg, idx); break;
				case 14: rep<14>(msg, idx); break;
				case 15: rep<15>(msg, idx); break;
				case 16: rep<16>(msg, idx); break;
				case 17: rep<17>(msg, idx); break;
				case 18: rep<18>(msg, idx); break;
				case 19: rep<19>(msg, idx); break;
				case 20: rep<20>(msg, idx); break;
				case 21: rep<21>(msg, idx); break;
				case 22: rep<22>(msg, idx); break;
				case 23: rep<23>(msg, idx); break;
				case 24: rep<24>(msg, idx); break;
				case 25: rep<25>(msg, idx); break;
				case 26: rep<26>(msg, idx); break;
				case 27: rep<27>(msg, idx); break;
				case 28: rep<28>(msg, idx); break;
				case 29: rep<29>(msg, idx); break;
				case 30: rep<30>(msg, idx); break;
				default: assert(false);
				} ++msg; break;
			case 6: switch (lvl) {
				case 1: spc<1>(msg, idx); break;
				case 2: spc<2>(msg, idx); break;
				case 3: spc<3>(msg, idx); break;
				case 4: spc<4>(msg, idx); break;
				case 5: spc<5>(msg, idx); break;
				case 6: spc<6>(msg, idx); break;
				case 7: spc<7>(msg, idx); break;
				case 8: spc<8>(msg, idx); break;
				case 9: spc<9>(msg, idx); break;
				case 10: spc<10>(msg, idx); break;
				case 11: spc<11>(msg, idx); break;
				case 12: spc<12>(msg, idx); break;
				case 13: spc<13>(msg, idx); break;
				case 14: spc<14>(msg, idx); break;
				case 15: spc<15>(msg, idx); break;
				case 16: spc<16>(msg, idx); break;
				case 17: spc<17>(msg, idx); break;
				case 18: spc<18>(msg, idx); break;
				case 19: spc<19>(msg, idx); break;
				case 20: spc<20>(msg, idx); break;
				case 21: spc<21>(msg, idx); break;
				case 22: spc<22>(msg, idx); break;
				case 23: spc<23>(msg, idx); break;
				case 24: spc<24>(msg, idx); break;
				case 25: spc<25>(msg, idx); break;
				case 26: spc<26>(msg, idx); break;
				case 27: spc<27>(msg, idx); break;
				case 28: spc<28>(msg, idx); break;
				case 29: spc<29>(msg, idx); break;
				case 30: spc<30>(msg, idx); break;
				default: assert(false);
				} msg += (1<<lvl)-1; break;
			default: assert(false);
			}
		}
		assert(lvl == level);
	}
};

