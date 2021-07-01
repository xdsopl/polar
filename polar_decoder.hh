/*
Successive cancellation decoding of polar codes

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#pragma once

template <typename TYPE, int MAX_M>
class PolarDecoder
{
	typedef PolarHelper<TYPE> PH;

	template <int level>
	static void left(TYPE *soft, TYPE *, TYPE *)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			soft[i+length/2] = PH::prod(soft[i+length], soft[i+length/2+length]);
	}
	template <int level>
	static void right(TYPE *soft, TYPE *hard, TYPE *)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			soft[i+length/2] = PH::madd(hard[i], soft[i+length], soft[i+length/2+length]);
	}
	template <int level>
	static void comb(TYPE *, TYPE *hard, TYPE *)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			hard[i] = PH::qmul(hard[i], hard[i+length/2]);
	}
	template <int level>
	static void rate0(TYPE *, TYPE *hard, TYPE *)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length; ++i)
			hard[i] = PH::one();
	}
	template <int level>
	static void rate1(TYPE *soft, TYPE *hard, TYPE *mesg)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length; ++i)
			hard[i] = PH::signum(soft[i+length]);
		for (int i = 0; i < length; i += 2) {
			soft[i] = PH::qmul(hard[i], hard[i+1]);
			soft[i+1] = hard[i+1];
		}
		for (int h = 2; h < length; h *= 2)
			for (int i = 0; i < length; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					soft[j] = PH::qmul(soft[j], soft[j+h]);
		for (int i = 0; i < length; ++i)
			mesg[i] = soft[i];
	}
	template <int level>
	static void rep(TYPE *soft, TYPE *hard, TYPE *mesg)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int h = length; h; h /= 2)
			for (int i = 0; i < h/2; ++i)
				soft[i+h/2] = PH::qadd(soft[i+h], soft[i+h/2+h]);
		TYPE hardi = PH::signum(soft[1]);
		*mesg = hardi;
		for (int i = 0; i < length; ++i)
			hard[i] = hardi;
	}
	template <int level>
	static void spc(TYPE *soft, TYPE *hard, TYPE *mesg)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length; ++i)
			hard[i] = PH::decide(soft[i+length]);
		TYPE parity = hard[0];
		for (int i = 1; i < length; ++i)
			parity = PH::qmul(parity, hard[i]);
		for (int i = 0; i < length; ++i)
			soft[i] = PH::qabs(soft[i+length]);
		TYPE weak = soft[0];
		for (int i = 1; i < length; ++i)
			weak = PH::qmin(weak, soft[i]);
		for (int i = 0; i < length; ++i)
			hard[i] = PH::flip(hard[i], parity, weak, soft[i]);
		for (int i = 0; i < length; i += 2) {
			soft[i] = PH::qmul(hard[i], hard[i+1]);
			soft[i+1] = hard[i+1];
		}
		for (int h = 2; h < length; h *= 2)
			for (int i = 0; i < length; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					soft[j] = PH::qmul(soft[j], soft[j+h]);
		for (int i = 0; i < length-1; ++i)
			mesg[i] = soft[i+1];
	}
	TYPE soft[1U<<(MAX_M+1)];
	TYPE hard[1U<<MAX_M];
public:
	void operator()(TYPE *message, const TYPE *codeword, const uint8_t *program)
	{
		TYPE *sft = soft, *hrd = hard, *msg = message;
		int level = *program++, lvl = level, length = 1 << level;
		assert(level <= MAX_M);
		for (int i = 0; i < length; ++i)
			soft[i+length] = codeword[i];
#if 1
#include "program.hh"
#else
		while (*program != 255) {
			switch (*program++) {
			case 0:	switch (lvl--) {
				case 2: left<2>(sft, hrd, msg); break;
				case 3: left<3>(sft, hrd, msg); break;
				case 4: left<4>(sft, hrd, msg); break;
				case 5: left<5>(sft, hrd, msg); break;
				case 6: left<6>(sft, hrd, msg); break;
				case 7: left<7>(sft, hrd, msg); break;
				case 8: left<8>(sft, hrd, msg); break;
				case 9: left<9>(sft, hrd, msg); break;
				case 10: left<10>(sft, hrd, msg); break;
				case 11: left<11>(sft, hrd, msg); break;
				case 12: left<12>(sft, hrd, msg); break;
				case 13: left<13>(sft, hrd, msg); break;
				case 14: left<14>(sft, hrd, msg); break;
				case 15: left<15>(sft, hrd, msg); break;
				case 16: left<16>(sft, hrd, msg); break;
				case 17: left<17>(sft, hrd, msg); break;
				case 18: left<18>(sft, hrd, msg); break;
				case 19: left<19>(sft, hrd, msg); break;
				case 20: left<20>(sft, hrd, msg); break;
				case 21: left<21>(sft, hrd, msg); break;
				case 22: left<22>(sft, hrd, msg); break;
				case 23: left<23>(sft, hrd, msg); break;
				case 24: left<24>(sft, hrd, msg); break;
				case 25: left<25>(sft, hrd, msg); break;
				case 26: left<26>(sft, hrd, msg); break;
				case 27: left<27>(sft, hrd, msg); break;
				case 28: left<28>(sft, hrd, msg); break;
				case 29: left<29>(sft, hrd, msg); break;
				case 30: left<30>(sft, hrd, msg); break;
				case 31: left<31>(sft, hrd, msg); break;
				default: assert(false);
				} break;
			case 1: switch (lvl+1) {
				case 2: right<2>(sft, hrd, msg); break;
				case 3: right<3>(sft, hrd, msg); break;
				case 4: right<4>(sft, hrd, msg); break;
				case 5: right<5>(sft, hrd, msg); break;
				case 6: right<6>(sft, hrd, msg); break;
				case 7: right<7>(sft, hrd, msg); break;
				case 8: right<8>(sft, hrd, msg); break;
				case 9: right<9>(sft, hrd, msg); break;
				case 10: right<10>(sft, hrd, msg); break;
				case 11: right<11>(sft, hrd, msg); break;
				case 12: right<12>(sft, hrd, msg); break;
				case 13: right<13>(sft, hrd, msg); break;
				case 14: right<14>(sft, hrd, msg); break;
				case 15: right<15>(sft, hrd, msg); break;
				case 16: right<16>(sft, hrd, msg); break;
				case 17: right<17>(sft, hrd, msg); break;
				case 18: right<18>(sft, hrd, msg); break;
				case 19: right<19>(sft, hrd, msg); break;
				case 20: right<20>(sft, hrd, msg); break;
				case 21: right<21>(sft, hrd, msg); break;
				case 22: right<22>(sft, hrd, msg); break;
				case 23: right<23>(sft, hrd, msg); break;
				case 24: right<24>(sft, hrd, msg); break;
				case 25: right<25>(sft, hrd, msg); break;
				case 26: right<26>(sft, hrd, msg); break;
				case 27: right<27>(sft, hrd, msg); break;
				case 28: right<28>(sft, hrd, msg); break;
				case 29: right<29>(sft, hrd, msg); break;
				case 30: right<30>(sft, hrd, msg); break;
				case 31: right<31>(sft, hrd, msg); break;
				default: assert(false);
				} hrd += 1<<lvl; break;
			case 2: hrd -= 1<<lvl; switch (++lvl) {
				case 2: comb<2>(sft, hrd, msg); break;
				case 3: comb<3>(sft, hrd, msg); break;
				case 4: comb<4>(sft, hrd, msg); break;
				case 5: comb<5>(sft, hrd, msg); break;
				case 6: comb<6>(sft, hrd, msg); break;
				case 7: comb<7>(sft, hrd, msg); break;
				case 8: comb<8>(sft, hrd, msg); break;
				case 9: comb<9>(sft, hrd, msg); break;
				case 10: comb<10>(sft, hrd, msg); break;
				case 11: comb<11>(sft, hrd, msg); break;
				case 12: comb<12>(sft, hrd, msg); break;
				case 13: comb<13>(sft, hrd, msg); break;
				case 14: comb<14>(sft, hrd, msg); break;
				case 15: comb<15>(sft, hrd, msg); break;
				case 16: comb<16>(sft, hrd, msg); break;
				case 17: comb<17>(sft, hrd, msg); break;
				case 18: comb<18>(sft, hrd, msg); break;
				case 19: comb<19>(sft, hrd, msg); break;
				case 20: comb<20>(sft, hrd, msg); break;
				case 21: comb<21>(sft, hrd, msg); break;
				case 22: comb<22>(sft, hrd, msg); break;
				case 23: comb<23>(sft, hrd, msg); break;
				case 24: comb<24>(sft, hrd, msg); break;
				case 25: comb<25>(sft, hrd, msg); break;
				case 26: comb<26>(sft, hrd, msg); break;
				case 27: comb<27>(sft, hrd, msg); break;
				case 28: comb<28>(sft, hrd, msg); break;
				case 29: comb<29>(sft, hrd, msg); break;
				case 30: comb<30>(sft, hrd, msg); break;
				case 31: comb<31>(sft, hrd, msg); break;
				default: assert(false);
				} break;
			case 3: switch (lvl) {
				case 1: rate0<1>(sft, hrd, msg); break;
				case 2: rate0<2>(sft, hrd, msg); break;
				case 3: rate0<3>(sft, hrd, msg); break;
				case 4: rate0<4>(sft, hrd, msg); break;
				case 5: rate0<5>(sft, hrd, msg); break;
				case 6: rate0<6>(sft, hrd, msg); break;
				case 7: rate0<7>(sft, hrd, msg); break;
				case 8: rate0<8>(sft, hrd, msg); break;
				case 9: rate0<9>(sft, hrd, msg); break;
				case 10: rate0<10>(sft, hrd, msg); break;
				case 11: rate0<11>(sft, hrd, msg); break;
				case 12: rate0<12>(sft, hrd, msg); break;
				case 13: rate0<13>(sft, hrd, msg); break;
				case 14: rate0<14>(sft, hrd, msg); break;
				case 15: rate0<15>(sft, hrd, msg); break;
				case 16: rate0<16>(sft, hrd, msg); break;
				case 17: rate0<17>(sft, hrd, msg); break;
				case 18: rate0<18>(sft, hrd, msg); break;
				case 19: rate0<19>(sft, hrd, msg); break;
				case 20: rate0<20>(sft, hrd, msg); break;
				case 21: rate0<21>(sft, hrd, msg); break;
				case 22: rate0<22>(sft, hrd, msg); break;
				case 23: rate0<23>(sft, hrd, msg); break;
				case 24: rate0<24>(sft, hrd, msg); break;
				case 25: rate0<25>(sft, hrd, msg); break;
				case 26: rate0<26>(sft, hrd, msg); break;
				case 27: rate0<27>(sft, hrd, msg); break;
				case 28: rate0<28>(sft, hrd, msg); break;
				case 29: rate0<29>(sft, hrd, msg); break;
				case 30: rate0<30>(sft, hrd, msg); break;
				default: assert(false);
				} break;
			case 4: switch (lvl) {
				case 1: rate1<1>(sft, hrd, msg); break;
				case 2: rate1<2>(sft, hrd, msg); break;
				case 3: rate1<3>(sft, hrd, msg); break;
				case 4: rate1<4>(sft, hrd, msg); break;
				case 5: rate1<5>(sft, hrd, msg); break;
				case 6: rate1<6>(sft, hrd, msg); break;
				case 7: rate1<7>(sft, hrd, msg); break;
				case 8: rate1<8>(sft, hrd, msg); break;
				case 9: rate1<9>(sft, hrd, msg); break;
				case 10: rate1<10>(sft, hrd, msg); break;
				case 11: rate1<11>(sft, hrd, msg); break;
				case 12: rate1<12>(sft, hrd, msg); break;
				case 13: rate1<13>(sft, hrd, msg); break;
				case 14: rate1<14>(sft, hrd, msg); break;
				case 15: rate1<15>(sft, hrd, msg); break;
				case 16: rate1<16>(sft, hrd, msg); break;
				case 17: rate1<17>(sft, hrd, msg); break;
				case 18: rate1<18>(sft, hrd, msg); break;
				case 19: rate1<19>(sft, hrd, msg); break;
				case 20: rate1<20>(sft, hrd, msg); break;
				case 21: rate1<21>(sft, hrd, msg); break;
				case 22: rate1<22>(sft, hrd, msg); break;
				case 23: rate1<23>(sft, hrd, msg); break;
				case 24: rate1<24>(sft, hrd, msg); break;
				case 25: rate1<25>(sft, hrd, msg); break;
				case 26: rate1<26>(sft, hrd, msg); break;
				case 27: rate1<27>(sft, hrd, msg); break;
				case 28: rate1<28>(sft, hrd, msg); break;
				case 29: rate1<29>(sft, hrd, msg); break;
				case 30: rate1<30>(sft, hrd, msg); break;
				default: assert(false);
				} msg += 1<<lvl; break;
			case 5: switch (lvl) {
				case 1: rep<1>(sft, hrd, msg); break;
				case 2: rep<2>(sft, hrd, msg); break;
				case 3: rep<3>(sft, hrd, msg); break;
				case 4: rep<4>(sft, hrd, msg); break;
				case 5: rep<5>(sft, hrd, msg); break;
				case 6: rep<6>(sft, hrd, msg); break;
				case 7: rep<7>(sft, hrd, msg); break;
				case 8: rep<8>(sft, hrd, msg); break;
				case 9: rep<9>(sft, hrd, msg); break;
				case 10: rep<10>(sft, hrd, msg); break;
				case 11: rep<11>(sft, hrd, msg); break;
				case 12: rep<12>(sft, hrd, msg); break;
				case 13: rep<13>(sft, hrd, msg); break;
				case 14: rep<14>(sft, hrd, msg); break;
				case 15: rep<15>(sft, hrd, msg); break;
				case 16: rep<16>(sft, hrd, msg); break;
				case 17: rep<17>(sft, hrd, msg); break;
				case 18: rep<18>(sft, hrd, msg); break;
				case 19: rep<19>(sft, hrd, msg); break;
				case 20: rep<20>(sft, hrd, msg); break;
				case 21: rep<21>(sft, hrd, msg); break;
				case 22: rep<22>(sft, hrd, msg); break;
				case 23: rep<23>(sft, hrd, msg); break;
				case 24: rep<24>(sft, hrd, msg); break;
				case 25: rep<25>(sft, hrd, msg); break;
				case 26: rep<26>(sft, hrd, msg); break;
				case 27: rep<27>(sft, hrd, msg); break;
				case 28: rep<28>(sft, hrd, msg); break;
				case 29: rep<29>(sft, hrd, msg); break;
				case 30: rep<30>(sft, hrd, msg); break;
				default: assert(false);
				} ++msg; break;
			case 6: switch (lvl) {
				case 1: spc<1>(sft, hrd, msg); break;
				case 2: spc<2>(sft, hrd, msg); break;
				case 3: spc<3>(sft, hrd, msg); break;
				case 4: spc<4>(sft, hrd, msg); break;
				case 5: spc<5>(sft, hrd, msg); break;
				case 6: spc<6>(sft, hrd, msg); break;
				case 7: spc<7>(sft, hrd, msg); break;
				case 8: spc<8>(sft, hrd, msg); break;
				case 9: spc<9>(sft, hrd, msg); break;
				case 10: spc<10>(sft, hrd, msg); break;
				case 11: spc<11>(sft, hrd, msg); break;
				case 12: spc<12>(sft, hrd, msg); break;
				case 13: spc<13>(sft, hrd, msg); break;
				case 14: spc<14>(sft, hrd, msg); break;
				case 15: spc<15>(sft, hrd, msg); break;
				case 16: spc<16>(sft, hrd, msg); break;
				case 17: spc<17>(sft, hrd, msg); break;
				case 18: spc<18>(sft, hrd, msg); break;
				case 19: spc<19>(sft, hrd, msg); break;
				case 20: spc<20>(sft, hrd, msg); break;
				case 21: spc<21>(sft, hrd, msg); break;
				case 22: spc<22>(sft, hrd, msg); break;
				case 23: spc<23>(sft, hrd, msg); break;
				case 24: spc<24>(sft, hrd, msg); break;
				case 25: spc<25>(sft, hrd, msg); break;
				case 26: spc<26>(sft, hrd, msg); break;
				case 27: spc<27>(sft, hrd, msg); break;
				case 28: spc<28>(sft, hrd, msg); break;
				case 29: spc<29>(sft, hrd, msg); break;
				case 30: spc<30>(sft, hrd, msg); break;
				default: assert(false);
				} msg += (1<<lvl)-1; break;
			default: assert(false);
			}
		}
#endif
		assert(lvl == level);
	}
};

