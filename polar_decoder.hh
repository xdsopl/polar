/*
Simplified successive cancellation list decoding of polar codes

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#pragma once

template <typename TYPE, int MAX_M>
class PolarDecoder
{
	typedef PolarHelper<TYPE> PH;
	typedef typename TYPE::value_type value_type;
	typedef typename PH::PATH PATH;

	template <int level>
	void left()
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			soft[i+length/2] = PH::prod(soft[i+length], soft[i+length/2+length]);
	}
	template <int level>
	void right(int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			soft[i+length/2] = PH::madd(hard[index+i], soft[i+length], soft[i+length/2+length]);
	}
	template <int level>
	void comb(int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			hard[index+i] = PH::qmul(hard[index+i], hard[index+i+length/2]);
	}
	template <int level>
	void rate0(PATH *path, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length; ++i)
			hard[index+i] = PH::one();
		for (int i = 0; i < length; ++i)
			for (int k = 0; k < TYPE::SIZE; ++k)
				if (soft[i+length].v[k] < 0)
					path[k] -= soft[i+length].v[k];
	}
	template <int level>
	void rate1(PATH *path, TYPE *mesg, int count, int index)
	{
		(void)path;
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
			mesg[count+i] = soft[i];
	}
	template <>
	void rate1<0>(PATH *path, TYPE *mesg, int count, int index)
	{
		PATH fork[2*TYPE::SIZE];
		for (int k = 0; k < TYPE::SIZE; ++k)
			fork[k] = fork[k+TYPE::SIZE] = path[k];
		for (int k = 0; k < TYPE::SIZE; ++k)
			if (soft[1].v[k] < 0)
				fork[k] -= soft[1].v[k];
			else
				fork[k+TYPE::SIZE] += soft[1].v[k];
		int perm[2*TYPE::SIZE];
		for (int k = 0; k < 2*TYPE::SIZE; ++k)
			perm[k] = k;
		std::nth_element(perm, perm+TYPE::SIZE, perm+2*TYPE::SIZE, [fork](int a, int b){ return fork[a] < fork[b]; });
		for (int k = 0; k < TYPE::SIZE; ++k)
			path[k] = fork[perm[k]];
		SIMD<typename TYPE::uint_type, TYPE::SIZE> map;
		for (int k = 0; k < TYPE::SIZE; ++k)
			map.v[k] = perm[k] % TYPE::SIZE;
		for (int i = 0; i < (1<<MAX_M); ++i)
			soft[i] = vshuf(soft[i], map);
		for (int i = 0; i < (1<<MAX_M); ++i)
			hard[i] = vshuf(hard[i], map);
		for (int i = 0; i < count; ++i)
			mesg[i] = vshuf(mesg[i], map);
		TYPE hrd;
		for (int k = 0; k < TYPE::SIZE; ++k)
			hrd.v[k] = perm[k] < TYPE::SIZE ? 1 : -1;
		mesg[count] = hrd;
		hard[index] = hrd;
	}
	template <int level>
	void rep(PATH *path, TYPE *mesg, int count, int index)
	{
		assert(level <= MAX_M);
		int length = 1 << level;
		PATH fork[2*TYPE::SIZE];
		for (int k = 0; k < TYPE::SIZE; ++k)
			fork[k] = fork[k+TYPE::SIZE] = path[k];
		for (int i = 0; i < length; ++i)
			for (int k = 0; k < TYPE::SIZE; ++k)
				if (soft[i+length].v[k] < 0)
					fork[k] -= soft[i+length].v[k];
				else
					fork[k+TYPE::SIZE] += soft[i+length].v[k];
		int perm[2*TYPE::SIZE];
		for (int k = 0; k < 2*TYPE::SIZE; ++k)
			perm[k] = k;
		std::nth_element(perm, perm+TYPE::SIZE, perm+2*TYPE::SIZE, [fork](int a, int b){ return fork[a] < fork[b]; });
		for (int k = 0; k < TYPE::SIZE; ++k)
			path[k] = fork[perm[k]];
		SIMD<typename TYPE::uint_type, TYPE::SIZE> map;
		for (int k = 0; k < TYPE::SIZE; ++k)
			map.v[k] = perm[k] % TYPE::SIZE;
		for (int i = 0; i < (1<<MAX_M); ++i)
			soft[i] = vshuf(soft[i], map);
		for (int i = 0; i < (1<<MAX_M); ++i)
			hard[i] = vshuf(hard[i], map);
		for (int i = 0; i < count; ++i)
			mesg[i] = vshuf(mesg[i], map);
		TYPE hrd;
		for (int k = 0; k < TYPE::SIZE; ++k)
			hrd.v[k] = perm[k] < TYPE::SIZE ? 1 : -1;
		mesg[count] = hrd;
		for (int i = 0; i < length; ++i)
			hard[index+i] = hrd;
	}
	TYPE soft[1U<<(MAX_M+1)];
	TYPE hard[1U<<MAX_M];
public:
	void operator()(PATH *metric, TYPE *message, const value_type *codeword, const uint8_t *program)
	{
		int level = *program++;
		assert(level <= MAX_M);
		int length = 1 << level;
		for (int i = 0; i < length; ++i)
			soft[i+length] = vdup<TYPE>(codeword[i]);
		int idx = 0, cnt = 0;
		PATH *met = metric;
		met[0] = 0;
		for (int k = 1; k < TYPE::SIZE; ++k)
			met[k] = 1000;
		TYPE *msg = message;
		while (*program != 255) {
			switch (*program++) {
			case (0<<5)+0: left<0+1>(); break;
			case (0<<5)+1: left<1+1>(); break;
			case (0<<5)+2: left<2+1>(); break;
			case (0<<5)+3: left<3+1>(); break;
			case (0<<5)+4: left<4+1>(); break;
			case (0<<5)+5: left<5+1>(); break;
			case (0<<5)+6: left<6+1>(); break;
			case (0<<5)+7: left<7+1>(); break;
			case (0<<5)+8: left<8+1>(); break;
			case (0<<5)+9: left<9+1>(); break;
			case (0<<5)+10: left<10+1>(); break;
			case (0<<5)+11: left<11+1>(); break;
			case (0<<5)+12: left<12+1>(); break;
			case (0<<5)+13: left<13+1>(); break;
			case (0<<5)+14: left<14+1>(); break;
			case (0<<5)+15: left<15+1>(); break;
			case (0<<5)+16: left<16+1>(); break;
			case (0<<5)+17: left<17+1>(); break;
			case (0<<5)+18: left<18+1>(); break;
			case (0<<5)+19: left<19+1>(); break;
			case (0<<5)+20: left<20+1>(); break;
			case (0<<5)+21: left<21+1>(); break;
			case (0<<5)+22: left<22+1>(); break;
			case (0<<5)+23: left<23+1>(); break;
			case (0<<5)+24: left<24+1>(); break;
			case (0<<5)+25: left<25+1>(); break;
			case (0<<5)+26: left<26+1>(); break;
			case (0<<5)+27: left<27+1>(); break;
			case (0<<5)+28: left<28+1>(); break;
			case (0<<5)+29: left<29+1>(); break;
			case (1<<5)+0: right<0+1>(idx); idx += 1<<0; break;
			case (1<<5)+1: right<1+1>(idx); idx += 1<<1; break;
			case (1<<5)+2: right<2+1>(idx); idx += 1<<2; break;
			case (1<<5)+3: right<3+1>(idx); idx += 1<<3; break;
			case (1<<5)+4: right<4+1>(idx); idx += 1<<4; break;
			case (1<<5)+5: right<5+1>(idx); idx += 1<<5; break;
			case (1<<5)+6: right<6+1>(idx); idx += 1<<6; break;
			case (1<<5)+7: right<7+1>(idx); idx += 1<<7; break;
			case (1<<5)+8: right<8+1>(idx); idx += 1<<8; break;
			case (1<<5)+9: right<9+1>(idx); idx += 1<<9; break;
			case (1<<5)+10: right<10+1>(idx); idx += 1<<10; break;
			case (1<<5)+11: right<11+1>(idx); idx += 1<<11; break;
			case (1<<5)+12: right<12+1>(idx); idx += 1<<12; break;
			case (1<<5)+13: right<13+1>(idx); idx += 1<<13; break;
			case (1<<5)+14: right<14+1>(idx); idx += 1<<14; break;
			case (1<<5)+15: right<15+1>(idx); idx += 1<<15; break;
			case (1<<5)+16: right<16+1>(idx); idx += 1<<16; break;
			case (1<<5)+17: right<17+1>(idx); idx += 1<<17; break;
			case (1<<5)+18: right<18+1>(idx); idx += 1<<18; break;
			case (1<<5)+19: right<19+1>(idx); idx += 1<<19; break;
			case (1<<5)+20: right<20+1>(idx); idx += 1<<20; break;
			case (1<<5)+21: right<21+1>(idx); idx += 1<<21; break;
			case (1<<5)+22: right<22+1>(idx); idx += 1<<22; break;
			case (1<<5)+23: right<23+1>(idx); idx += 1<<23; break;
			case (1<<5)+24: right<24+1>(idx); idx += 1<<24; break;
			case (1<<5)+25: right<25+1>(idx); idx += 1<<25; break;
			case (1<<5)+26: right<26+1>(idx); idx += 1<<26; break;
			case (1<<5)+27: right<27+1>(idx); idx += 1<<27; break;
			case (1<<5)+28: right<28+1>(idx); idx += 1<<28; break;
			case (1<<5)+29: right<29+1>(idx); idx += 1<<29; break;
			case (2<<5)+0: idx -= 1<<0; comb<0+1>(idx); break;
			case (2<<5)+1: idx -= 1<<1; comb<1+1>(idx); break;
			case (2<<5)+2: idx -= 1<<2; comb<2+1>(idx); break;
			case (2<<5)+3: idx -= 1<<3; comb<3+1>(idx); break;
			case (2<<5)+4: idx -= 1<<4; comb<4+1>(idx); break;
			case (2<<5)+5: idx -= 1<<5; comb<5+1>(idx); break;
			case (2<<5)+6: idx -= 1<<6; comb<6+1>(idx); break;
			case (2<<5)+7: idx -= 1<<7; comb<7+1>(idx); break;
			case (2<<5)+8: idx -= 1<<8; comb<8+1>(idx); break;
			case (2<<5)+9: idx -= 1<<9; comb<9+1>(idx); break;
			case (2<<5)+10: idx -= 1<<10; comb<10+1>(idx); break;
			case (2<<5)+11: idx -= 1<<11; comb<11+1>(idx); break;
			case (2<<5)+12: idx -= 1<<12; comb<12+1>(idx); break;
			case (2<<5)+13: idx -= 1<<13; comb<13+1>(idx); break;
			case (2<<5)+14: idx -= 1<<14; comb<14+1>(idx); break;
			case (2<<5)+15: idx -= 1<<15; comb<15+1>(idx); break;
			case (2<<5)+16: idx -= 1<<16; comb<16+1>(idx); break;
			case (2<<5)+17: idx -= 1<<17; comb<17+1>(idx); break;
			case (2<<5)+18: idx -= 1<<18; comb<18+1>(idx); break;
			case (2<<5)+19: idx -= 1<<19; comb<19+1>(idx); break;
			case (2<<5)+20: idx -= 1<<20; comb<20+1>(idx); break;
			case (2<<5)+21: idx -= 1<<21; comb<21+1>(idx); break;
			case (2<<5)+22: idx -= 1<<22; comb<22+1>(idx); break;
			case (2<<5)+23: idx -= 1<<23; comb<23+1>(idx); break;
			case (2<<5)+24: idx -= 1<<24; comb<24+1>(idx); break;
			case (2<<5)+25: idx -= 1<<25; comb<25+1>(idx); break;
			case (2<<5)+26: idx -= 1<<26; comb<26+1>(idx); break;
			case (2<<5)+27: idx -= 1<<27; comb<27+1>(idx); break;
			case (2<<5)+28: idx -= 1<<28; comb<28+1>(idx); break;
			case (2<<5)+29: idx -= 1<<29; comb<29+1>(idx); break;
			case (3<<5)+0: rate0<0>(met, idx); break;
			case (3<<5)+1: rate0<1>(met, idx); break;
			case (3<<5)+2: rate0<2>(met, idx); break;
			case (3<<5)+3: rate0<3>(met, idx); break;
			case (3<<5)+4: rate0<4>(met, idx); break;
			case (3<<5)+5: rate0<5>(met, idx); break;
			case (3<<5)+6: rate0<6>(met, idx); break;
			case (3<<5)+7: rate0<7>(met, idx); break;
			case (3<<5)+8: rate0<8>(met, idx); break;
			case (3<<5)+9: rate0<9>(met, idx); break;
			case (3<<5)+10: rate0<10>(met, idx); break;
			case (3<<5)+11: rate0<11>(met, idx); break;
			case (3<<5)+12: rate0<12>(met, idx); break;
			case (3<<5)+13: rate0<13>(met, idx); break;
			case (3<<5)+14: rate0<14>(met, idx); break;
			case (3<<5)+15: rate0<15>(met, idx); break;
			case (3<<5)+16: rate0<16>(met, idx); break;
			case (3<<5)+17: rate0<17>(met, idx); break;
			case (3<<5)+18: rate0<18>(met, idx); break;
			case (3<<5)+19: rate0<19>(met, idx); break;
			case (3<<5)+20: rate0<20>(met, idx); break;
			case (3<<5)+21: rate0<21>(met, idx); break;
			case (3<<5)+22: rate0<22>(met, idx); break;
			case (3<<5)+23: rate0<23>(met, idx); break;
			case (3<<5)+24: rate0<24>(met, idx); break;
			case (3<<5)+25: rate0<25>(met, idx); break;
			case (3<<5)+26: rate0<26>(met, idx); break;
			case (3<<5)+27: rate0<27>(met, idx); break;
			case (3<<5)+28: rate0<28>(met, idx); break;
			case (3<<5)+29: rate0<29>(met, idx); break;
			case (3<<5)+30: rate0<30>(met, idx); break;
			case (4<<5)+0: rate1<0>(met, msg, cnt, idx); cnt += 1<<0; break;
			case (4<<5)+1: rate1<1>(met, msg, cnt, idx); cnt += 1<<1; break;
			case (4<<5)+2: rate1<2>(met, msg, cnt, idx); cnt += 1<<2; break;
			case (4<<5)+3: rate1<3>(met, msg, cnt, idx); cnt += 1<<3; break;
			case (4<<5)+4: rate1<4>(met, msg, cnt, idx); cnt += 1<<4; break;
			case (4<<5)+5: rate1<5>(met, msg, cnt, idx); cnt += 1<<5; break;
			case (4<<5)+6: rate1<6>(met, msg, cnt, idx); cnt += 1<<6; break;
			case (4<<5)+7: rate1<7>(met, msg, cnt, idx); cnt += 1<<7; break;
			case (4<<5)+8: rate1<8>(met, msg, cnt, idx); cnt += 1<<8; break;
			case (4<<5)+9: rate1<9>(met, msg, cnt, idx); cnt += 1<<9; break;
			case (4<<5)+10: rate1<10>(met, msg, cnt, idx); cnt += 1<<10; break;
			case (4<<5)+11: rate1<11>(met, msg, cnt, idx); cnt += 1<<11; break;
			case (4<<5)+12: rate1<12>(met, msg, cnt, idx); cnt += 1<<12; break;
			case (4<<5)+13: rate1<13>(met, msg, cnt, idx); cnt += 1<<13; break;
			case (4<<5)+14: rate1<14>(met, msg, cnt, idx); cnt += 1<<14; break;
			case (4<<5)+15: rate1<15>(met, msg, cnt, idx); cnt += 1<<15; break;
			case (4<<5)+16: rate1<16>(met, msg, cnt, idx); cnt += 1<<16; break;
			case (4<<5)+17: rate1<17>(met, msg, cnt, idx); cnt += 1<<17; break;
			case (4<<5)+18: rate1<18>(met, msg, cnt, idx); cnt += 1<<18; break;
			case (4<<5)+19: rate1<19>(met, msg, cnt, idx); cnt += 1<<19; break;
			case (4<<5)+20: rate1<20>(met, msg, cnt, idx); cnt += 1<<20; break;
			case (4<<5)+21: rate1<21>(met, msg, cnt, idx); cnt += 1<<21; break;
			case (4<<5)+22: rate1<22>(met, msg, cnt, idx); cnt += 1<<22; break;
			case (4<<5)+23: rate1<23>(met, msg, cnt, idx); cnt += 1<<23; break;
			case (4<<5)+24: rate1<24>(met, msg, cnt, idx); cnt += 1<<24; break;
			case (4<<5)+25: rate1<25>(met, msg, cnt, idx); cnt += 1<<25; break;
			case (4<<5)+26: rate1<26>(met, msg, cnt, idx); cnt += 1<<26; break;
			case (4<<5)+27: rate1<27>(met, msg, cnt, idx); cnt += 1<<27; break;
			case (4<<5)+28: rate1<28>(met, msg, cnt, idx); cnt += 1<<28; break;
			case (4<<5)+29: rate1<29>(met, msg, cnt, idx); cnt += 1<<29; break;
			case (4<<5)+30: rate1<30>(met, msg, cnt, idx); cnt += 1<<30; break;
			case (5<<5)+0: rep<0>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+1: rep<1>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+2: rep<2>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+3: rep<3>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+4: rep<4>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+5: rep<5>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+6: rep<6>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+7: rep<7>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+8: rep<8>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+9: rep<9>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+10: rep<10>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+11: rep<11>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+12: rep<12>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+13: rep<13>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+14: rep<14>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+15: rep<15>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+16: rep<16>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+17: rep<17>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+18: rep<18>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+19: rep<19>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+20: rep<20>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+21: rep<21>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+22: rep<22>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+23: rep<23>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+24: rep<24>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+25: rep<25>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+26: rep<26>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+27: rep<27>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+28: rep<28>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+29: rep<29>(met, msg, cnt, idx); ++cnt; break;
			case (5<<5)+30: rep<30>(met, msg, cnt, idx); ++cnt; break;
			default:
				assert(false);
			}
		}
	}
};

