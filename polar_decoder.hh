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
	void left(TYPE *, int)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			soft[i+length/2] = PH::prod(soft[i+length], soft[i+length/2+length]);
	}
	template <int level>
	void right(TYPE *, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			soft[i+length/2] = PH::madd(hard[index+i], soft[i+length], soft[i+length/2+length]);
	}
	template <int level>
	void comb(TYPE *, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			hard[index+i] = PH::qmul(hard[index+i], hard[index+i+length/2]);
	}
	template <int level>
	void rate0(TYPE *, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length; ++i)
			hard[index+i] = PH::one();
	}
	template <int level>
	void rate1(TYPE *mesg, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length; ++i)
			hard[index+i] = PH::signum(soft[i+length]);
		for (int i = 0; i < length; i += 2) {
			soft[i] = PH::qmul(hard[index+i], hard[index+i+1]);
			soft[i+1] = hard[index+i+1];
		}
		for (int h = 2; h < length; h *= 2)
			for (int i = 0; i < length; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					soft[j] = PH::qmul(soft[j], soft[j+h]);
		for (int i = 0; i < length; ++i)
			mesg[i] = soft[i];
	}
	template <int level>
	void rep(TYPE *mesg, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int h = length; h; h /= 2)
			for (int i = 0; i < h/2; ++i)
				soft[i+h/2] = PH::qadd(soft[i+h], soft[i+h/2+h]);
		TYPE hardi = PH::signum(soft[1]);
		*mesg = hardi;
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
		TYPE *msg = message;
		while (*program != 255) {
			switch (*program++) {
			case (0<<5)+0: left<0+2>(msg, idx); break;
			case (0<<5)+1: left<1+2>(msg, idx); break;
			case (0<<5)+2: left<2+2>(msg, idx); break;
			case (0<<5)+3: left<3+2>(msg, idx); break;
			case (0<<5)+4: left<4+2>(msg, idx); break;
			case (0<<5)+5: left<5+2>(msg, idx); break;
			case (0<<5)+6: left<6+2>(msg, idx); break;
			case (0<<5)+7: left<7+2>(msg, idx); break;
			case (0<<5)+8: left<8+2>(msg, idx); break;
			case (0<<5)+9: left<9+2>(msg, idx); break;
			case (0<<5)+10: left<10+2>(msg, idx); break;
			case (0<<5)+11: left<11+2>(msg, idx); break;
			case (0<<5)+12: left<12+2>(msg, idx); break;
			case (0<<5)+13: left<13+2>(msg, idx); break;
			case (0<<5)+14: left<14+2>(msg, idx); break;
			case (0<<5)+15: left<15+2>(msg, idx); break;
			case (0<<5)+16: left<16+2>(msg, idx); break;
			case (0<<5)+17: left<17+2>(msg, idx); break;
			case (0<<5)+18: left<18+2>(msg, idx); break;
			case (0<<5)+19: left<19+2>(msg, idx); break;
			case (0<<5)+20: left<20+2>(msg, idx); break;
			case (0<<5)+21: left<21+2>(msg, idx); break;
			case (0<<5)+22: left<22+2>(msg, idx); break;
			case (0<<5)+23: left<23+2>(msg, idx); break;
			case (0<<5)+24: left<24+2>(msg, idx); break;
			case (0<<5)+25: left<25+2>(msg, idx); break;
			case (0<<5)+26: left<26+2>(msg, idx); break;
			case (0<<5)+27: left<27+2>(msg, idx); break;
			case (0<<5)+28: left<28+2>(msg, idx); break;
			case (0<<5)+29: left<29+2>(msg, idx); break;
			case (1<<5)+0: right<0+2>(msg, idx); idx += 1<<(0+1); break;
			case (1<<5)+1: right<1+2>(msg, idx); idx += 1<<(1+1); break;
			case (1<<5)+2: right<2+2>(msg, idx); idx += 1<<(2+1); break;
			case (1<<5)+3: right<3+2>(msg, idx); idx += 1<<(3+1); break;
			case (1<<5)+4: right<4+2>(msg, idx); idx += 1<<(4+1); break;
			case (1<<5)+5: right<5+2>(msg, idx); idx += 1<<(5+1); break;
			case (1<<5)+6: right<6+2>(msg, idx); idx += 1<<(6+1); break;
			case (1<<5)+7: right<7+2>(msg, idx); idx += 1<<(7+1); break;
			case (1<<5)+8: right<8+2>(msg, idx); idx += 1<<(8+1); break;
			case (1<<5)+9: right<9+2>(msg, idx); idx += 1<<(9+1); break;
			case (1<<5)+10: right<10+2>(msg, idx); idx += 1<<(10+1); break;
			case (1<<5)+11: right<11+2>(msg, idx); idx += 1<<(11+1); break;
			case (1<<5)+12: right<12+2>(msg, idx); idx += 1<<(12+1); break;
			case (1<<5)+13: right<13+2>(msg, idx); idx += 1<<(13+1); break;
			case (1<<5)+14: right<14+2>(msg, idx); idx += 1<<(14+1); break;
			case (1<<5)+15: right<15+2>(msg, idx); idx += 1<<(15+1); break;
			case (1<<5)+16: right<16+2>(msg, idx); idx += 1<<(16+1); break;
			case (1<<5)+17: right<17+2>(msg, idx); idx += 1<<(17+1); break;
			case (1<<5)+18: right<18+2>(msg, idx); idx += 1<<(18+1); break;
			case (1<<5)+19: right<19+2>(msg, idx); idx += 1<<(19+1); break;
			case (1<<5)+20: right<20+2>(msg, idx); idx += 1<<(20+1); break;
			case (1<<5)+21: right<21+2>(msg, idx); idx += 1<<(21+1); break;
			case (1<<5)+22: right<22+2>(msg, idx); idx += 1<<(22+1); break;
			case (1<<5)+23: right<23+2>(msg, idx); idx += 1<<(23+1); break;
			case (1<<5)+24: right<24+2>(msg, idx); idx += 1<<(24+1); break;
			case (1<<5)+25: right<25+2>(msg, idx); idx += 1<<(25+1); break;
			case (1<<5)+26: right<26+2>(msg, idx); idx += 1<<(26+1); break;
			case (1<<5)+27: right<27+2>(msg, idx); idx += 1<<(27+1); break;
			case (1<<5)+28: right<28+2>(msg, idx); idx += 1<<(28+1); break;
			case (1<<5)+29: right<29+2>(msg, idx); idx += 1<<(29+1); break;
			case (2<<5)+0: idx -= 1<<(0+1); comb<0+2>(msg, idx); break;
			case (2<<5)+1: idx -= 1<<(1+1); comb<1+2>(msg, idx); break;
			case (2<<5)+2: idx -= 1<<(2+1); comb<2+2>(msg, idx); break;
			case (2<<5)+3: idx -= 1<<(3+1); comb<3+2>(msg, idx); break;
			case (2<<5)+4: idx -= 1<<(4+1); comb<4+2>(msg, idx); break;
			case (2<<5)+5: idx -= 1<<(5+1); comb<5+2>(msg, idx); break;
			case (2<<5)+6: idx -= 1<<(6+1); comb<6+2>(msg, idx); break;
			case (2<<5)+7: idx -= 1<<(7+1); comb<7+2>(msg, idx); break;
			case (2<<5)+8: idx -= 1<<(8+1); comb<8+2>(msg, idx); break;
			case (2<<5)+9: idx -= 1<<(9+1); comb<9+2>(msg, idx); break;
			case (2<<5)+10: idx -= 1<<(10+1); comb<10+2>(msg, idx); break;
			case (2<<5)+11: idx -= 1<<(11+1); comb<11+2>(msg, idx); break;
			case (2<<5)+12: idx -= 1<<(12+1); comb<12+2>(msg, idx); break;
			case (2<<5)+13: idx -= 1<<(13+1); comb<13+2>(msg, idx); break;
			case (2<<5)+14: idx -= 1<<(14+1); comb<14+2>(msg, idx); break;
			case (2<<5)+15: idx -= 1<<(15+1); comb<15+2>(msg, idx); break;
			case (2<<5)+16: idx -= 1<<(16+1); comb<16+2>(msg, idx); break;
			case (2<<5)+17: idx -= 1<<(17+1); comb<17+2>(msg, idx); break;
			case (2<<5)+18: idx -= 1<<(18+1); comb<18+2>(msg, idx); break;
			case (2<<5)+19: idx -= 1<<(19+1); comb<19+2>(msg, idx); break;
			case (2<<5)+20: idx -= 1<<(20+1); comb<20+2>(msg, idx); break;
			case (2<<5)+21: idx -= 1<<(21+1); comb<21+2>(msg, idx); break;
			case (2<<5)+22: idx -= 1<<(22+1); comb<22+2>(msg, idx); break;
			case (2<<5)+23: idx -= 1<<(23+1); comb<23+2>(msg, idx); break;
			case (2<<5)+24: idx -= 1<<(24+1); comb<24+2>(msg, idx); break;
			case (2<<5)+25: idx -= 1<<(25+1); comb<25+2>(msg, idx); break;
			case (2<<5)+26: idx -= 1<<(26+1); comb<26+2>(msg, idx); break;
			case (2<<5)+27: idx -= 1<<(27+1); comb<27+2>(msg, idx); break;
			case (2<<5)+28: idx -= 1<<(28+1); comb<28+2>(msg, idx); break;
			case (2<<5)+29: idx -= 1<<(29+1); comb<29+2>(msg, idx); break;
			case (3<<5)+0: rate0<0+1>(msg, idx); break;
			case (3<<5)+1: rate0<1+1>(msg, idx); break;
			case (3<<5)+2: rate0<2+1>(msg, idx); break;
			case (3<<5)+3: rate0<3+1>(msg, idx); break;
			case (3<<5)+4: rate0<4+1>(msg, idx); break;
			case (3<<5)+5: rate0<5+1>(msg, idx); break;
			case (3<<5)+6: rate0<6+1>(msg, idx); break;
			case (3<<5)+7: rate0<7+1>(msg, idx); break;
			case (3<<5)+8: rate0<8+1>(msg, idx); break;
			case (3<<5)+9: rate0<9+1>(msg, idx); break;
			case (3<<5)+10: rate0<10+1>(msg, idx); break;
			case (3<<5)+11: rate0<11+1>(msg, idx); break;
			case (3<<5)+12: rate0<12+1>(msg, idx); break;
			case (3<<5)+13: rate0<13+1>(msg, idx); break;
			case (3<<5)+14: rate0<14+1>(msg, idx); break;
			case (3<<5)+15: rate0<15+1>(msg, idx); break;
			case (3<<5)+16: rate0<16+1>(msg, idx); break;
			case (3<<5)+17: rate0<17+1>(msg, idx); break;
			case (3<<5)+18: rate0<18+1>(msg, idx); break;
			case (3<<5)+19: rate0<19+1>(msg, idx); break;
			case (3<<5)+20: rate0<20+1>(msg, idx); break;
			case (3<<5)+21: rate0<21+1>(msg, idx); break;
			case (3<<5)+22: rate0<22+1>(msg, idx); break;
			case (3<<5)+23: rate0<23+1>(msg, idx); break;
			case (3<<5)+24: rate0<24+1>(msg, idx); break;
			case (3<<5)+25: rate0<25+1>(msg, idx); break;
			case (3<<5)+26: rate0<26+1>(msg, idx); break;
			case (3<<5)+27: rate0<27+1>(msg, idx); break;
			case (3<<5)+28: rate0<28+1>(msg, idx); break;
			case (3<<5)+29: rate0<29+1>(msg, idx); break;
			case (3<<5)+30: rate0<30+1>(msg, idx); break;
			case (4<<5)+0: rate1<0+1>(msg, idx); msg += 1<<(0+1); break;
			case (4<<5)+1: rate1<1+1>(msg, idx); msg += 1<<(1+1); break;
			case (4<<5)+2: rate1<2+1>(msg, idx); msg += 1<<(2+1); break;
			case (4<<5)+3: rate1<3+1>(msg, idx); msg += 1<<(3+1); break;
			case (4<<5)+4: rate1<4+1>(msg, idx); msg += 1<<(4+1); break;
			case (4<<5)+5: rate1<5+1>(msg, idx); msg += 1<<(5+1); break;
			case (4<<5)+6: rate1<6+1>(msg, idx); msg += 1<<(6+1); break;
			case (4<<5)+7: rate1<7+1>(msg, idx); msg += 1<<(7+1); break;
			case (4<<5)+8: rate1<8+1>(msg, idx); msg += 1<<(8+1); break;
			case (4<<5)+9: rate1<9+1>(msg, idx); msg += 1<<(9+1); break;
			case (4<<5)+10: rate1<10+1>(msg, idx); msg += 1<<(10+1); break;
			case (4<<5)+11: rate1<11+1>(msg, idx); msg += 1<<(11+1); break;
			case (4<<5)+12: rate1<12+1>(msg, idx); msg += 1<<(12+1); break;
			case (4<<5)+13: rate1<13+1>(msg, idx); msg += 1<<(13+1); break;
			case (4<<5)+14: rate1<14+1>(msg, idx); msg += 1<<(14+1); break;
			case (4<<5)+15: rate1<15+1>(msg, idx); msg += 1<<(15+1); break;
			case (4<<5)+16: rate1<16+1>(msg, idx); msg += 1<<(16+1); break;
			case (4<<5)+17: rate1<17+1>(msg, idx); msg += 1<<(17+1); break;
			case (4<<5)+18: rate1<18+1>(msg, idx); msg += 1<<(18+1); break;
			case (4<<5)+19: rate1<19+1>(msg, idx); msg += 1<<(19+1); break;
			case (4<<5)+20: rate1<20+1>(msg, idx); msg += 1<<(20+1); break;
			case (4<<5)+21: rate1<21+1>(msg, idx); msg += 1<<(21+1); break;
			case (4<<5)+22: rate1<22+1>(msg, idx); msg += 1<<(22+1); break;
			case (4<<5)+23: rate1<23+1>(msg, idx); msg += 1<<(23+1); break;
			case (4<<5)+24: rate1<24+1>(msg, idx); msg += 1<<(24+1); break;
			case (4<<5)+25: rate1<25+1>(msg, idx); msg += 1<<(25+1); break;
			case (4<<5)+26: rate1<26+1>(msg, idx); msg += 1<<(26+1); break;
			case (4<<5)+27: rate1<27+1>(msg, idx); msg += 1<<(27+1); break;
			case (4<<5)+28: rate1<28+1>(msg, idx); msg += 1<<(28+1); break;
			case (4<<5)+29: rate1<29+1>(msg, idx); msg += 1<<(29+1); break;
			case (4<<5)+30: rate1<30+1>(msg, idx); msg += 1<<(30+1); break;
			case (5<<5)+0: rep<0+1>(msg, idx); ++msg; break;
			case (5<<5)+1: rep<1+1>(msg, idx); ++msg; break;
			case (5<<5)+2: rep<2+1>(msg, idx); ++msg; break;
			case (5<<5)+3: rep<3+1>(msg, idx); ++msg; break;
			case (5<<5)+4: rep<4+1>(msg, idx); ++msg; break;
			case (5<<5)+5: rep<5+1>(msg, idx); ++msg; break;
			case (5<<5)+6: rep<6+1>(msg, idx); ++msg; break;
			case (5<<5)+7: rep<7+1>(msg, idx); ++msg; break;
			case (5<<5)+8: rep<8+1>(msg, idx); ++msg; break;
			case (5<<5)+9: rep<9+1>(msg, idx); ++msg; break;
			case (5<<5)+10: rep<10+1>(msg, idx); ++msg; break;
			case (5<<5)+11: rep<11+1>(msg, idx); ++msg; break;
			case (5<<5)+12: rep<12+1>(msg, idx); ++msg; break;
			case (5<<5)+13: rep<13+1>(msg, idx); ++msg; break;
			case (5<<5)+14: rep<14+1>(msg, idx); ++msg; break;
			case (5<<5)+15: rep<15+1>(msg, idx); ++msg; break;
			case (5<<5)+16: rep<16+1>(msg, idx); ++msg; break;
			case (5<<5)+17: rep<17+1>(msg, idx); ++msg; break;
			case (5<<5)+18: rep<18+1>(msg, idx); ++msg; break;
			case (5<<5)+19: rep<19+1>(msg, idx); ++msg; break;
			case (5<<5)+20: rep<20+1>(msg, idx); ++msg; break;
			case (5<<5)+21: rep<21+1>(msg, idx); ++msg; break;
			case (5<<5)+22: rep<22+1>(msg, idx); ++msg; break;
			case (5<<5)+23: rep<23+1>(msg, idx); ++msg; break;
			case (5<<5)+24: rep<24+1>(msg, idx); ++msg; break;
			case (5<<5)+25: rep<25+1>(msg, idx); ++msg; break;
			case (5<<5)+26: rep<26+1>(msg, idx); ++msg; break;
			case (5<<5)+27: rep<27+1>(msg, idx); ++msg; break;
			case (5<<5)+28: rep<28+1>(msg, idx); ++msg; break;
			case (5<<5)+29: rep<29+1>(msg, idx); ++msg; break;
			case (5<<5)+30: rep<30+1>(msg, idx); ++msg; break;
			default:
				assert(false);
			}
		}
	}
};

