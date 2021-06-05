/*
Successive cancellation decoding of polar codes

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#pragma once

template <typename TYPE, int MAX_M>
class PolarDecoder
{
	static const int U = 2;
	typedef PolarHelper<TYPE> PH;

	void leaf0(TYPE **msg, int index)
	{
#if 0
		TYPE soft20 = soft[(1<<U)+0];
		TYPE soft21 = soft[(1<<U)+1];
		TYPE soft22 = soft[(1<<U)+2];
		TYPE soft23 = soft[(1<<U)+3];
		TYPE soft10 = PH::prod(soft20, soft22);
		TYPE soft11 = PH::prod(soft21, soft23);
		TYPE soft00 = PH::prod(soft10, soft11);
		TYPE hard0 = PH::signum(soft00);
		*(*msg)++ = hard0;
		TYPE soft01 = PH::madd(hard0, soft10, soft11);
		TYPE hard1 = PH::signum(soft01);
		*(*msg)++ = hard1;
		hard0 = PH::qmul(hard0, hard1);
		TYPE soft12 = PH::madd(hard0, soft20, soft22);
		TYPE soft13 = PH::madd(hard1, soft21, soft23);
		TYPE soft02 = PH::prod(soft12, soft13);
		TYPE hard2 = PH::signum(soft02);
		*(*msg)++ = hard2;
		TYPE soft03 = PH::madd(hard2, soft12, soft13);
		TYPE hard3 = PH::signum(soft03);
		*(*msg)++ = hard3;
		hard2 = PH::qmul(hard2, hard3);
		hard0 = PH::qmul(hard0, hard2);
		hard1 = PH::qmul(hard1, hard3);
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
#else
		TYPE hard0 = PH::signum(soft[(1<<U)+0]);
		TYPE hard1 = PH::signum(soft[(1<<U)+1]);
		TYPE hard2 = PH::signum(soft[(1<<U)+2]);
		TYPE hard3 = PH::signum(soft[(1<<U)+3]);
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
		hard0 = PH::qmul(hard0, hard1);
		hard2 = PH::qmul(hard2, hard3);
		hard0 = PH::qmul(hard0, hard2);
		hard1 = PH::qmul(hard1, hard3);
		*(*msg)++ = hard0;
		*(*msg)++ = hard1;
		*(*msg)++ = hard2;
		*(*msg)++ = hard3;
#endif
	}
	void leaf1(TYPE **msg, int index)
	{
		TYPE soft20 = soft[(1<<U)+0];
		TYPE soft21 = soft[(1<<U)+1];
		TYPE soft22 = soft[(1<<U)+2];
		TYPE soft23 = soft[(1<<U)+3];
		TYPE soft10 = PH::prod(soft20, soft22);
		TYPE soft11 = PH::prod(soft21, soft23);
		TYPE hard0 = PH::one();
		TYPE soft01 = PH::qadd(soft10, soft11);
		TYPE hard1 = PH::signum(soft01);
		*(*msg)++ = hard1;
		hard0 = PH::qmul(hard0, hard1);
		TYPE soft12 = PH::madd(hard0, soft20, soft22);
		TYPE soft13 = PH::madd(hard1, soft21, soft23);
		TYPE soft02 = PH::prod(soft12, soft13);
		TYPE hard2 = PH::signum(soft02);
		*(*msg)++ = hard2;
		TYPE soft03 = PH::madd(hard2, soft12, soft13);
		TYPE hard3 = PH::signum(soft03);
		*(*msg)++ = hard3;
		hard2 = PH::qmul(hard2, hard3);
		hard0 = PH::qmul(hard0, hard2);
		hard1 = PH::qmul(hard1, hard3);
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf2(TYPE **msg, int index)
	{
		TYPE soft20 = soft[(1<<U)+0];
		TYPE soft21 = soft[(1<<U)+1];
		TYPE soft22 = soft[(1<<U)+2];
		TYPE soft23 = soft[(1<<U)+3];
		TYPE soft10 = PH::prod(soft20, soft22);
		TYPE soft11 = PH::prod(soft21, soft23);
		TYPE soft00 = PH::prod(soft10, soft11);
		TYPE hard0 = PH::signum(soft00);
		*(*msg)++ = hard0;
		TYPE hard1 = PH::one();
		TYPE soft12 = PH::madd(hard0, soft20, soft22);
		TYPE soft13 = PH::madd(hard1, soft21, soft23);
		TYPE soft02 = PH::prod(soft12, soft13);
		TYPE hard2 = PH::signum(soft02);
		*(*msg)++ = hard2;
		TYPE soft03 = PH::madd(hard2, soft12, soft13);
		TYPE hard3 = PH::signum(soft03);
		*(*msg)++ = hard3;
		hard2 = PH::qmul(hard2, hard3);
		hard0 = PH::qmul(hard0, hard2);
		hard1 = PH::qmul(hard1, hard3);
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf3(TYPE **msg, int index)
	{
		TYPE soft20 = soft[(1<<U)+0];
		TYPE soft21 = soft[(1<<U)+1];
		TYPE soft22 = soft[(1<<U)+2];
		TYPE soft23 = soft[(1<<U)+3];
		TYPE hard0 = PH::one();
		TYPE hard1 = PH::one();
		TYPE soft12 = PH::qadd(soft20, soft22);
		TYPE soft13 = PH::qadd(soft21, soft23);
		TYPE soft02 = PH::prod(soft12, soft13);
		TYPE hard2 = PH::signum(soft02);
		*(*msg)++ = hard2;
		TYPE soft03 = PH::madd(hard2, soft12, soft13);
		TYPE hard3 = PH::signum(soft03);
		*(*msg)++ = hard3;
		hard2 = PH::qmul(hard2, hard3);
		hard0 = PH::qmul(hard0, hard2);
		hard1 = PH::qmul(hard1, hard3);
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf4(TYPE **msg, int index)
	{
		TYPE soft20 = soft[(1<<U)+0];
		TYPE soft21 = soft[(1<<U)+1];
		TYPE soft22 = soft[(1<<U)+2];
		TYPE soft23 = soft[(1<<U)+3];
		TYPE soft10 = PH::prod(soft20, soft22);
		TYPE soft11 = PH::prod(soft21, soft23);
		TYPE soft00 = PH::prod(soft10, soft11);
		TYPE hard0 = PH::signum(soft00);
		*(*msg)++ = hard0;
		TYPE soft01 = PH::madd(hard0, soft10, soft11);
		TYPE hard1 = PH::signum(soft01);
		*(*msg)++ = hard1;
		hard0 = PH::qmul(hard0, hard1);
		TYPE soft12 = PH::madd(hard0, soft20, soft22);
		TYPE soft13 = PH::madd(hard1, soft21, soft23);
		TYPE hard2 = PH::one();
		TYPE soft03 = PH::qadd(soft12, soft13);
		TYPE hard3 = PH::signum(soft03);
		*(*msg)++ = hard3;
		hard2 = PH::qmul(hard2, hard3);
		hard0 = PH::qmul(hard0, hard2);
		hard1 = PH::qmul(hard1, hard3);
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf5(TYPE **msg, int index)
	{
		TYPE soft20 = soft[(1<<U)+0];
		TYPE soft21 = soft[(1<<U)+1];
		TYPE soft22 = soft[(1<<U)+2];
		TYPE soft23 = soft[(1<<U)+3];
		TYPE soft10 = PH::prod(soft20, soft22);
		TYPE soft11 = PH::prod(soft21, soft23);
		TYPE hard0 = PH::one();
		TYPE soft01 = PH::qadd(soft10, soft11);
		TYPE hard1 = PH::signum(soft01);
		*(*msg)++ = hard1;
		hard0 = PH::qmul(hard0, hard1);
		TYPE soft12 = PH::madd(hard0, soft20, soft22);
		TYPE soft13 = PH::madd(hard1, soft21, soft23);
		TYPE hard2 = PH::one();
		TYPE soft03 = PH::qadd(soft12, soft13);
		TYPE hard3 = PH::signum(soft03);
		*(*msg)++ = hard3;
		hard2 = PH::qmul(hard2, hard3);
		hard0 = PH::qmul(hard0, hard2);
		hard1 = PH::qmul(hard1, hard3);
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf6(TYPE **msg, int index)
	{
		TYPE soft20 = soft[(1<<U)+0];
		TYPE soft21 = soft[(1<<U)+1];
		TYPE soft22 = soft[(1<<U)+2];
		TYPE soft23 = soft[(1<<U)+3];
		TYPE soft10 = PH::prod(soft20, soft22);
		TYPE soft11 = PH::prod(soft21, soft23);
		TYPE soft00 = PH::prod(soft10, soft11);
		TYPE hard0 = PH::signum(soft00);
		*(*msg)++ = hard0;
		TYPE hard1 = PH::one();
		TYPE soft12 = PH::madd(hard0, soft20, soft22);
		TYPE soft13 = PH::madd(hard1, soft21, soft23);
		TYPE hard2 = PH::one();
		TYPE soft03 = PH::qadd(soft12, soft13);
		TYPE hard3 = PH::signum(soft03);
		*(*msg)++ = hard3;
		hard2 = PH::qmul(hard2, hard3);
		hard0 = PH::qmul(hard0, hard2);
		hard1 = PH::qmul(hard1, hard3);
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf7(TYPE **msg, int index)
	{
#if 0
		TYPE soft20 = soft[(1<<U)+0];
		TYPE soft21 = soft[(1<<U)+1];
		TYPE soft22 = soft[(1<<U)+2];
		TYPE soft23 = soft[(1<<U)+3];
		TYPE hard0 = PH::one();
		TYPE hard1 = PH::one();
		TYPE soft12 = PH::qadd(soft20, soft22);
		TYPE soft13 = PH::qadd(soft21, soft23);
		TYPE hard2 = PH::one();
		TYPE soft03 = PH::qadd(soft12, soft13);
		TYPE hard3 = PH::signum(soft03);
		*(*msg)++ = hard3;
		hard2 = PH::qmul(hard2, hard3);
		hard0 = PH::qmul(hard0, hard2);
		hard1 = PH::qmul(hard1, hard3);
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
#else
		TYPE soft20 = soft[(1<<U)+0];
		TYPE soft21 = soft[(1<<U)+1];
		TYPE soft22 = soft[(1<<U)+2];
		TYPE soft23 = soft[(1<<U)+3];
		TYPE soft12 = PH::qadd(soft20, soft22);
		TYPE soft13 = PH::qadd(soft21, soft23);
		TYPE soft03 = PH::qadd(soft12, soft13);
		TYPE hard3 = PH::signum(soft03);
		*(*msg)++ = hard3;
		hard[index+0] = hard3;
		hard[index+1] = hard3;
		hard[index+2] = hard3;
		hard[index+3] = hard3;
#endif
	}
	void leaf8(TYPE **msg, int index)
	{
		TYPE soft20 = soft[(1<<U)+0];
		TYPE soft21 = soft[(1<<U)+1];
		TYPE soft22 = soft[(1<<U)+2];
		TYPE soft23 = soft[(1<<U)+3];
		TYPE soft10 = PH::prod(soft20, soft22);
		TYPE soft11 = PH::prod(soft21, soft23);
		TYPE soft00 = PH::prod(soft10, soft11);
		TYPE hard0 = PH::signum(soft00);
		*(*msg)++ = hard0;
		TYPE soft01 = PH::madd(hard0, soft10, soft11);
		TYPE hard1 = PH::signum(soft01);
		*(*msg)++ = hard1;
		hard0 = PH::qmul(hard0, hard1);
		TYPE soft12 = PH::madd(hard0, soft20, soft22);
		TYPE soft13 = PH::madd(hard1, soft21, soft23);
		TYPE soft02 = PH::prod(soft12, soft13);
		TYPE hard2 = PH::signum(soft02);
		*(*msg)++ = hard2;
		TYPE hard3 = PH::one();
		hard0 = PH::qmul(hard0, hard2);
		hard1 = PH::qmul(hard1, hard3);
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf9(TYPE **msg, int index)
	{
		TYPE soft20 = soft[(1<<U)+0];
		TYPE soft21 = soft[(1<<U)+1];
		TYPE soft22 = soft[(1<<U)+2];
		TYPE soft23 = soft[(1<<U)+3];
		TYPE soft10 = PH::prod(soft20, soft22);
		TYPE soft11 = PH::prod(soft21, soft23);
		TYPE hard0 = PH::one();
		TYPE soft01 = PH::qadd(soft10, soft11);
		TYPE hard1 = PH::signum(soft01);
		*(*msg)++ = hard1;
		hard0 = PH::qmul(hard0, hard1);
		TYPE soft12 = PH::madd(hard0, soft20, soft22);
		TYPE soft13 = PH::madd(hard1, soft21, soft23);
		TYPE soft02 = PH::prod(soft12, soft13);
		TYPE hard2 = PH::signum(soft02);
		*(*msg)++ = hard2;
		TYPE hard3 = PH::one();
		hard0 = PH::qmul(hard0, hard2);
		hard1 = PH::qmul(hard1, hard3);
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf10(TYPE **msg, int index)
	{
		TYPE soft20 = soft[(1<<U)+0];
		TYPE soft21 = soft[(1<<U)+1];
		TYPE soft22 = soft[(1<<U)+2];
		TYPE soft23 = soft[(1<<U)+3];
		TYPE soft10 = PH::prod(soft20, soft22);
		TYPE soft11 = PH::prod(soft21, soft23);
		TYPE soft00 = PH::prod(soft10, soft11);
		TYPE hard0 = PH::signum(soft00);
		*(*msg)++ = hard0;
		TYPE hard1 = PH::one();
		TYPE soft12 = PH::madd(hard0, soft20, soft22);
		TYPE soft13 = PH::madd(hard1, soft21, soft23);
		TYPE soft02 = PH::prod(soft12, soft13);
		TYPE hard2 = PH::signum(soft02);
		*(*msg)++ = hard2;
		TYPE hard3 = PH::one();
		hard0 = PH::qmul(hard0, hard2);
		hard1 = PH::qmul(hard1, hard3);
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf11(TYPE **msg, int index)
	{
		TYPE soft20 = soft[(1<<U)+0];
		TYPE soft21 = soft[(1<<U)+1];
		TYPE soft22 = soft[(1<<U)+2];
		TYPE soft23 = soft[(1<<U)+3];
		TYPE hard0 = PH::one();
		TYPE hard1 = PH::one();
		TYPE soft12 = PH::qadd(soft20, soft22);
		TYPE soft13 = PH::qadd(soft21, soft23);
		TYPE soft02 = PH::prod(soft12, soft13);
		TYPE hard2 = PH::signum(soft02);
		*(*msg)++ = hard2;
		TYPE hard3 = PH::one();
		hard0 = PH::qmul(hard0, hard2);
		hard1 = PH::qmul(hard1, hard3);
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf12(TYPE **msg, int index)
	{
		TYPE soft20 = soft[(1<<U)+0];
		TYPE soft21 = soft[(1<<U)+1];
		TYPE soft22 = soft[(1<<U)+2];
		TYPE soft23 = soft[(1<<U)+3];
		TYPE soft10 = PH::prod(soft20, soft22);
		TYPE soft11 = PH::prod(soft21, soft23);
		TYPE soft00 = PH::prod(soft10, soft11);
		TYPE hard0 = PH::signum(soft00);
		*(*msg)++ = hard0;
		TYPE soft01 = PH::madd(hard0, soft10, soft11);
		TYPE hard1 = PH::signum(soft01);
		*(*msg)++ = hard1;
		hard0 = PH::qmul(hard0, hard1);
		TYPE hard2 = PH::one();
		TYPE hard3 = PH::one();
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf13(TYPE **msg, int index)
	{
		TYPE soft20 = soft[(1<<U)+0];
		TYPE soft21 = soft[(1<<U)+1];
		TYPE soft22 = soft[(1<<U)+2];
		TYPE soft23 = soft[(1<<U)+3];
		TYPE soft10 = PH::prod(soft20, soft22);
		TYPE soft11 = PH::prod(soft21, soft23);
		TYPE hard0 = PH::one();
		TYPE soft01 = PH::qadd(soft10, soft11);
		TYPE hard1 = PH::signum(soft01);
		*(*msg)++ = hard1;
		hard0 = PH::qmul(hard0, hard1);
		TYPE hard2 = PH::one();
		TYPE hard3 = PH::one();
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf14(TYPE **msg, int index)
	{
		TYPE soft20 = soft[(1<<U)+0];
		TYPE soft21 = soft[(1<<U)+1];
		TYPE soft22 = soft[(1<<U)+2];
		TYPE soft23 = soft[(1<<U)+3];
		TYPE soft10 = PH::prod(soft20, soft22);
		TYPE soft11 = PH::prod(soft21, soft23);
		TYPE soft00 = PH::prod(soft10, soft11);
		TYPE hard0 = PH::signum(soft00);
		*(*msg)++ = hard0;
		TYPE hard1 = PH::one();
		TYPE hard2 = PH::one();
		TYPE hard3 = PH::one();
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	void leaf15(TYPE **, int index)
	{
		TYPE hard0 = PH::one();
		TYPE hard1 = PH::one();
		TYPE hard2 = PH::one();
		TYPE hard3 = PH::one();
		hard[index+0] = hard0;
		hard[index+1] = hard1;
		hard[index+2] = hard2;
		hard[index+3] = hard3;
	}
	template <int level>
	void left(TYPE **, int)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			soft[i+length/2] = PH::prod(soft[i+length], soft[i+length/2+length]);
	}
	template <int level>
	void right(TYPE **, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			soft[i+length/2] = PH::madd(hard[index+i], soft[i+length], soft[i+length/2+length]);
	}
	template <int level>
	void rate0_right(TYPE **, int)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			soft[i+length/2] = PH::qadd(soft[i+length], soft[i+length/2+length]);
	}
	template <int level>
	void combine(TYPE **, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			hard[index+i] = PH::qmul(hard[index+i], hard[index+i+length/2]);
	}
	template <int level>
	void rate0_combine(TYPE **, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			hard[index+i] = hard[index+i+length/2];
	}
	template <int level>
	void rate1_combine(TYPE **msg, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			hard[index+i] = PH::qmul(hard[index+i], hard[index+i+length/2] = PH::signum(soft[i+length/2+length]));
		for (int i = 0; i < length/2; i += 2) {
			soft[i] = PH::qmul(hard[index+i+length/2], hard[index+i+1+length/2]);
			soft[i+1] = hard[index+i+1+length/2];
		}
		for (int h = 2; h < length/2; h *= 2)
			for (int i = 0; i < length; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					soft[j] = PH::qmul(soft[j], soft[j+h]);
		for (int i = 0; i < length/2; ++i)
			*(*msg)++ = soft[i];
	}
	template <int level>
	void rep(TYPE **msg, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int h = length; h; h /= 2)
			for (int i = 0; i < h/2; ++i)
				soft[i+h/2] = PH::qadd(soft[i+h], soft[i+h/2+h]);
		TYPE hardi = PH::signum(soft[1]);
		*(*msg)++ = hardi;
		for (int i = 0; i < length; ++i)
			hard[index+i] = hardi;
	}
	TYPE soft[1U<<(MAX_M+1)];
	TYPE hard[1U<<MAX_M];
public:
	void operator()(TYPE *message, const TYPE *codeword, const uint8_t *program)
	{
		int level = *program++;
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length; ++i)
			soft[i+length] = codeword[i];
		int idx = 0;
		TYPE **msg = &message;
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
			case (3<<5)+3: rate0_right<3>(msg, idx); idx += 1<<(3-1); break;
			case (3<<5)+4: rate0_right<4>(msg, idx); idx += 1<<(4-1); break;
			case (3<<5)+5: rate0_right<5>(msg, idx); idx += 1<<(5-1); break;
			case (3<<5)+6: rate0_right<6>(msg, idx); idx += 1<<(6-1); break;
			case (3<<5)+7: rate0_right<7>(msg, idx); idx += 1<<(7-1); break;
			case (3<<5)+8: rate0_right<8>(msg, idx); idx += 1<<(8-1); break;
			case (3<<5)+9: rate0_right<9>(msg, idx); idx += 1<<(9-1); break;
			case (3<<5)+10: rate0_right<10>(msg, idx); idx += 1<<(10-1); break;
			case (3<<5)+11: rate0_right<11>(msg, idx); idx += 1<<(11-1); break;
			case (3<<5)+12: rate0_right<12>(msg, idx); idx += 1<<(12-1); break;
			case (3<<5)+13: rate0_right<13>(msg, idx); idx += 1<<(13-1); break;
			case (3<<5)+14: rate0_right<14>(msg, idx); idx += 1<<(14-1); break;
			case (3<<5)+15: rate0_right<15>(msg, idx); idx += 1<<(15-1); break;
			case (3<<5)+16: rate0_right<16>(msg, idx); idx += 1<<(16-1); break;
			case (3<<5)+17: rate0_right<17>(msg, idx); idx += 1<<(17-1); break;
			case (3<<5)+18: rate0_right<18>(msg, idx); idx += 1<<(18-1); break;
			case (3<<5)+19: rate0_right<19>(msg, idx); idx += 1<<(19-1); break;
			case (3<<5)+20: rate0_right<20>(msg, idx); idx += 1<<(20-1); break;
			case (3<<5)+21: rate0_right<21>(msg, idx); idx += 1<<(21-1); break;
			case (3<<5)+22: rate0_right<22>(msg, idx); idx += 1<<(22-1); break;
			case (3<<5)+23: rate0_right<23>(msg, idx); idx += 1<<(23-1); break;
			case (3<<5)+24: rate0_right<24>(msg, idx); idx += 1<<(24-1); break;
			case (3<<5)+25: rate0_right<25>(msg, idx); idx += 1<<(25-1); break;
			case (3<<5)+26: rate0_right<26>(msg, idx); idx += 1<<(26-1); break;
			case (3<<5)+27: rate0_right<27>(msg, idx); idx += 1<<(27-1); break;
			case (3<<5)+28: rate0_right<28>(msg, idx); idx += 1<<(28-1); break;
			case (3<<5)+29: rate0_right<29>(msg, idx); idx += 1<<(29-1); break;
			case (3<<5)+30: rate0_right<30>(msg, idx); idx += 1<<(30-1); break;
			case (4<<5)+3: idx -= 1<<(3-1); combine<3>(msg, idx); break;
			case (4<<5)+4: idx -= 1<<(4-1); combine<4>(msg, idx); break;
			case (4<<5)+5: idx -= 1<<(5-1); combine<5>(msg, idx); break;
			case (4<<5)+6: idx -= 1<<(6-1); combine<6>(msg, idx); break;
			case (4<<5)+7: idx -= 1<<(7-1); combine<7>(msg, idx); break;
			case (4<<5)+8: idx -= 1<<(8-1); combine<8>(msg, idx); break;
			case (4<<5)+9: idx -= 1<<(9-1); combine<9>(msg, idx); break;
			case (4<<5)+10: idx -= 1<<(10-1); combine<10>(msg, idx); break;
			case (4<<5)+11: idx -= 1<<(11-1); combine<11>(msg, idx); break;
			case (4<<5)+12: idx -= 1<<(12-1); combine<12>(msg, idx); break;
			case (4<<5)+13: idx -= 1<<(13-1); combine<13>(msg, idx); break;
			case (4<<5)+14: idx -= 1<<(14-1); combine<14>(msg, idx); break;
			case (4<<5)+15: idx -= 1<<(15-1); combine<15>(msg, idx); break;
			case (4<<5)+16: idx -= 1<<(16-1); combine<16>(msg, idx); break;
			case (4<<5)+17: idx -= 1<<(17-1); combine<17>(msg, idx); break;
			case (4<<5)+18: idx -= 1<<(18-1); combine<18>(msg, idx); break;
			case (4<<5)+19: idx -= 1<<(19-1); combine<19>(msg, idx); break;
			case (4<<5)+20: idx -= 1<<(20-1); combine<20>(msg, idx); break;
			case (4<<5)+21: idx -= 1<<(21-1); combine<21>(msg, idx); break;
			case (4<<5)+22: idx -= 1<<(22-1); combine<22>(msg, idx); break;
			case (4<<5)+23: idx -= 1<<(23-1); combine<23>(msg, idx); break;
			case (4<<5)+24: idx -= 1<<(24-1); combine<24>(msg, idx); break;
			case (4<<5)+25: idx -= 1<<(25-1); combine<25>(msg, idx); break;
			case (4<<5)+26: idx -= 1<<(26-1); combine<26>(msg, idx); break;
			case (4<<5)+27: idx -= 1<<(27-1); combine<27>(msg, idx); break;
			case (4<<5)+28: idx -= 1<<(28-1); combine<28>(msg, idx); break;
			case (4<<5)+29: idx -= 1<<(29-1); combine<29>(msg, idx); break;
			case (4<<5)+30: idx -= 1<<(30-1); combine<30>(msg, idx); break;
			case (5<<5)+3: idx -= 1<<(3-1); rate0_combine<3>(msg, idx); break;
			case (5<<5)+4: idx -= 1<<(4-1); rate0_combine<4>(msg, idx); break;
			case (5<<5)+5: idx -= 1<<(5-1); rate0_combine<5>(msg, idx); break;
			case (5<<5)+6: idx -= 1<<(6-1); rate0_combine<6>(msg, idx); break;
			case (5<<5)+7: idx -= 1<<(7-1); rate0_combine<7>(msg, idx); break;
			case (5<<5)+8: idx -= 1<<(8-1); rate0_combine<8>(msg, idx); break;
			case (5<<5)+9: idx -= 1<<(9-1); rate0_combine<9>(msg, idx); break;
			case (5<<5)+10: idx -= 1<<(10-1); rate0_combine<10>(msg, idx); break;
			case (5<<5)+11: idx -= 1<<(11-1); rate0_combine<11>(msg, idx); break;
			case (5<<5)+12: idx -= 1<<(12-1); rate0_combine<12>(msg, idx); break;
			case (5<<5)+13: idx -= 1<<(13-1); rate0_combine<13>(msg, idx); break;
			case (5<<5)+14: idx -= 1<<(14-1); rate0_combine<14>(msg, idx); break;
			case (5<<5)+15: idx -= 1<<(15-1); rate0_combine<15>(msg, idx); break;
			case (5<<5)+16: idx -= 1<<(16-1); rate0_combine<16>(msg, idx); break;
			case (5<<5)+17: idx -= 1<<(17-1); rate0_combine<17>(msg, idx); break;
			case (5<<5)+18: idx -= 1<<(18-1); rate0_combine<18>(msg, idx); break;
			case (5<<5)+19: idx -= 1<<(19-1); rate0_combine<19>(msg, idx); break;
			case (5<<5)+20: idx -= 1<<(20-1); rate0_combine<20>(msg, idx); break;
			case (5<<5)+21: idx -= 1<<(21-1); rate0_combine<21>(msg, idx); break;
			case (5<<5)+22: idx -= 1<<(22-1); rate0_combine<22>(msg, idx); break;
			case (5<<5)+23: idx -= 1<<(23-1); rate0_combine<23>(msg, idx); break;
			case (5<<5)+24: idx -= 1<<(24-1); rate0_combine<24>(msg, idx); break;
			case (5<<5)+25: idx -= 1<<(25-1); rate0_combine<25>(msg, idx); break;
			case (5<<5)+26: idx -= 1<<(26-1); rate0_combine<26>(msg, idx); break;
			case (5<<5)+27: idx -= 1<<(27-1); rate0_combine<27>(msg, idx); break;
			case (5<<5)+28: idx -= 1<<(28-1); rate0_combine<28>(msg, idx); break;
			case (5<<5)+29: idx -= 1<<(29-1); rate0_combine<29>(msg, idx); break;
			case (5<<5)+30: idx -= 1<<(30-1); rate0_combine<30>(msg, idx); break;
			case (6<<5)+3: rate1_combine<3>(msg, idx); break;
			case (6<<5)+4: rate1_combine<4>(msg, idx); break;
			case (6<<5)+5: rate1_combine<5>(msg, idx); break;
			case (6<<5)+6: rate1_combine<6>(msg, idx); break;
			case (6<<5)+7: rate1_combine<7>(msg, idx); break;
			case (6<<5)+8: rate1_combine<8>(msg, idx); break;
			case (6<<5)+9: rate1_combine<9>(msg, idx); break;
			case (6<<5)+10: rate1_combine<10>(msg, idx); break;
			case (6<<5)+11: rate1_combine<11>(msg, idx); break;
			case (6<<5)+12: rate1_combine<12>(msg, idx); break;
			case (6<<5)+13: rate1_combine<13>(msg, idx); break;
			case (6<<5)+14: rate1_combine<14>(msg, idx); break;
			case (6<<5)+15: rate1_combine<15>(msg, idx); break;
			case (6<<5)+16: rate1_combine<16>(msg, idx); break;
			case (6<<5)+17: rate1_combine<17>(msg, idx); break;
			case (6<<5)+18: rate1_combine<18>(msg, idx); break;
			case (6<<5)+19: rate1_combine<19>(msg, idx); break;
			case (6<<5)+20: rate1_combine<20>(msg, idx); break;
			case (6<<5)+21: rate1_combine<21>(msg, idx); break;
			case (6<<5)+22: rate1_combine<22>(msg, idx); break;
			case (6<<5)+23: rate1_combine<23>(msg, idx); break;
			case (6<<5)+24: rate1_combine<24>(msg, idx); break;
			case (6<<5)+25: rate1_combine<25>(msg, idx); break;
			case (6<<5)+26: rate1_combine<26>(msg, idx); break;
			case (6<<5)+27: rate1_combine<27>(msg, idx); break;
			case (6<<5)+28: rate1_combine<28>(msg, idx); break;
			case (6<<5)+29: rate1_combine<29>(msg, idx); break;
			case (6<<5)+30: rate1_combine<30>(msg, idx); break;
			case (7<<5)+3: rep<3>(msg, idx); break;
			case (7<<5)+4: rep<4>(msg, idx); break;
			case (7<<5)+5: rep<5>(msg, idx); break;
			case (7<<5)+6: rep<6>(msg, idx); break;
			case (7<<5)+7: rep<7>(msg, idx); break;
			case (7<<5)+8: rep<8>(msg, idx); break;
			case (7<<5)+9: rep<9>(msg, idx); break;
			case (7<<5)+10: rep<10>(msg, idx); break;
			case (7<<5)+11: rep<11>(msg, idx); break;
			case (7<<5)+12: rep<12>(msg, idx); break;
			case (7<<5)+13: rep<13>(msg, idx); break;
			case (7<<5)+14: rep<14>(msg, idx); break;
			case (7<<5)+15: rep<15>(msg, idx); break;
			case (7<<5)+16: rep<16>(msg, idx); break;
			case (7<<5)+17: rep<17>(msg, idx); break;
			case (7<<5)+18: rep<18>(msg, idx); break;
			case (7<<5)+19: rep<19>(msg, idx); break;
			case (7<<5)+20: rep<20>(msg, idx); break;
			case (7<<5)+21: rep<21>(msg, idx); break;
			case (7<<5)+22: rep<22>(msg, idx); break;
			case (7<<5)+23: rep<23>(msg, idx); break;
			case (7<<5)+24: rep<24>(msg, idx); break;
			case (7<<5)+25: rep<25>(msg, idx); break;
			case (7<<5)+26: rep<26>(msg, idx); break;
			case (7<<5)+27: rep<27>(msg, idx); break;
			case (7<<5)+28: rep<28>(msg, idx); break;
			case (7<<5)+29: rep<29>(msg, idx); break;
			case (7<<5)+30: rep<30>(msg, idx); break;
			default:
				assert(false);
			}
		}
	}
};

