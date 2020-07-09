/*
Test bench for successive cancellation decoding of polar codes on the binary erasure channel

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#include <random>
#include <chrono>
#include <cassert>
#include <iostream>
#include <functional>

template <int M>
class PolarTransform
{
	static const int N = 1 << M;
public:
	void operator()(int8_t *output, const int8_t *input)
	{
		for (int i = 0; i < N; i += 2) {
			int in0 = input[i];
			int in1 = input[i+1];
			output[i] = in0 * in1;
			output[i+1] = in1;
		}
		for (int h = 2; h < N; h *= 2)
			for (int i = 0; i < N; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					output[j] *= output[j+h];
	}
};

template <int M>
class PolarEncoder
{
	static const int N = 1 << M;
	static const int U = 3;
public:
	void operator()(int8_t *codeword, const int8_t *message, const uint8_t *frozen)
	{
		for (int i = 0; i < N; i += 2) {
			int msg0 = frozen[i>>U]&(1<<(i&((1<<U)-1))) ? 1 : *message++;
			int msg1 = frozen[(i+1)>>U]&(1<<((i+1)&((1<<U)-1))) ? 1 : *message++;
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
	static const int U = 3;
public:
	void operator()(int8_t *codeword, const int8_t *message, const uint8_t *frozen)
	{
		for (int i = 0; i < N; i += 2) {
			int msg0 = frozen[i>>U]&(1<<(i&((1<<U)-1))) ? 1 : *message++;
			int msg1 = frozen[(i+1)>>U]&(1<<((i+1)&((1<<U)-1))) ? 1 : *message++;
			codeword[i] = msg0 * msg1;
			codeword[i+1] = msg1;
		}
		for (int h = 2; h < N; h *= 2)
			for (int i = 0; i < N; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					codeword[j] *= codeword[j+h];
		for (int i = 0; i < N; i += 2) {
			int msg0 = frozen[i>>U]&(1<<(i&((1<<U)-1))) ? 1 : codeword[i];
			int msg1 = frozen[(i+1)>>U]&(1<<((i+1)&((1<<U)-1))) ? 1 : codeword[i+1];
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
class PolarFreezer
{
	static const int N = 1 << M;
	static const int U = 3;
	static void freeze(uint8_t *bits, double pe, double th, int i, int h)
	{
		if (h) {
			freeze(bits, pe * (2-pe), th, i, h/2);
			freeze(bits, pe * pe, th, i+h, h/2);
		} else {
			bits[i>>U] |= (pe>th) << (i&((1<<U)-1));
		}
	}
public:
	int operator()(uint8_t *frozen_bits, double erasure_probability = 0.5, double freezing_threshold = 0.5)
	{
		for (int i = 0; i < N>>U; ++i)
			frozen_bits[i] = 0;
		freeze(frozen_bits, erasure_probability, freezing_threshold, 0, N / 2);
		int K = 0;
		for (int i = 0; i < N; ++i)
			K += !(frozen_bits[i>>U]&(1<<(i&((1<<U)-1))));
		return K;
	}
};

template<typename TYPE>
int popcnt(TYPE x)
{
	int cnt = 0;
	while (x) {
		++cnt;
		x &= x-1;
	}
	return cnt;
}

template <int M>
class PolarCompiler
{
	static const int N = 1 << M;
	static const int U = 3;
	static uint32_t leaf(int func, int index)
	{
		return (func << 23) | (index >> U);
	}
	static uint32_t node(int func, int level, int index)
	{
		return ((256 | (func << 5) | level) << 23) | (index >> U);
	}
	static uint32_t rate0(int level, int index)
	{
		return node(3, level, index);
	}
	static uint32_t rate1(int level, int index)
	{
		return node(4, level, index);
	}
	static int frozen_count(const uint8_t *frozen, int level)
	{
		int count = 0;
		for (int i = 0; i < 1<<(level-U); ++i)
			count += popcnt(frozen[i]);
		return count;
	}
	static void compile(uint32_t **program, const uint8_t *frozen, int index, int level)
	{
		if (level > U) {
			int count = frozen_count(frozen, level);
			if (count == 1<<level) {
				*(*program)++ = rate0(level, index);
			} else if (count == 0) {
				*(*program)++ = rate1(level, index);
			} else {
				*(*program)++ = node(0, level, index);
				compile(program, frozen, index, level-1);
				*(*program)++ = node(1, level, index);
				compile(program, frozen+(1<<(level-1-U)), index+(1<<(level-1)), level-1);
				*(*program)++ = node(2, level, index);
			}
		} else {
			*(*program)++ = leaf(*frozen, index);
		}
	}
public:
	void operator()(uint32_t *program, const uint8_t *frozen)
	{
		uint32_t *first = program;
		compile(&program, frozen, 0, M);
		*program++ = 0xffffffff;
		std::cerr << "program length = " << program - first << std::endl;
	}
};

template <int M>
class PolarDecoder
{
	static const int N = 1 << M;
	static const int U = 3;
	static int8_t signum(int8_t v)
	{
		return (v > 0) - (v < 0);
	}
	static int8_t qabs(int8_t a)
	{
		return std::abs(std::max<int8_t>(a, -127));
	}
	static int8_t qadd(int8_t a, int8_t b)
	{
		return std::min<int16_t>(std::max<int16_t>(int16_t(a) + int16_t(b), -128), 127);
	}
	static int8_t qmul(int8_t a, int8_t b)
	{
		return std::min<int16_t>(std::max<int16_t>(int16_t(a) * int16_t(b), -128), 127);
	}
	static int8_t prod(int8_t a, int8_t b)
	{
		return signum(a) * signum(b) * std::min(qabs(a), qabs(b));
	}
	static int8_t madd(int8_t a, int8_t b, int8_t c)
	{
		return qadd(qmul(a, b), c);
	}
	template <int level>
	void node0(int8_t **, int)
	{
		assert(level <= M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			soft[i+length/2] = prod(soft[i+length], soft[i+length/2+length]);
	}
	template <int level>
	void node1(int8_t **, int index)
	{
		assert(level <= M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			soft[i+length/2] = madd(hard[index+i], soft[i+length], soft[i+length/2+length]);
	}
	template <int level>
	void node2(int8_t **, int index)
	{
		assert(level <= M);
		int length = 1 << level;
		for (int i = 0; i < length/2; ++i)
			hard[index+i] *= hard[index+i+length/2];
	}
	template <int level>
	void rate0(int8_t **, int index)
	{
		assert(level <= M);
		int length = 1 << level;
		for (int i = 0; i < length; ++i)
			hard[index+i] = 1;
	}
	template <int level>
	void rate1(int8_t **msg, int index)
	{
		assert(level <= M);
		int length = 1 << level;
		for (int i = 0; i < length; ++i)
			hard[index+i] = signum(soft[i+length]);
		for (int i = 0; i < length; i += 2) {
			(*msg)[i] = hard[index+i] * hard[index+i+1];
			(*msg)[i+1] = hard[index+i+1];
		}
		for (int h = 2; h < length; h *= 2)
			for (int i = 0; i < length; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					(*msg)[j] *= (*msg)[j+h];
		*msg += length;
	}
#include "polar_leafs.hh"
	int8_t soft[2*N];
	int8_t hard[N];
	void decode(int8_t **msg, int func, int index)
	{
		switch (func) {
		case 0: leaf0(msg, index); break;
		case 1: leaf1(msg, index); break;
		case 2: leaf2(msg, index); break;
		case 3: leaf3(msg, index); break;
		case 4: leaf4(msg, index); break;
		case 5: leaf5(msg, index); break;
		case 6: leaf6(msg, index); break;
		case 7: leaf7(msg, index); break;
		case 8: leaf8(msg, index); break;
		case 9: leaf9(msg, index); break;
		case 10: leaf10(msg, index); break;
		case 11: leaf11(msg, index); break;
		case 12: leaf12(msg, index); break;
		case 13: leaf13(msg, index); break;
		case 14: leaf14(msg, index); break;
		case 15: leaf15(msg, index); break;
		case 16: leaf16(msg, index); break;
		case 17: leaf17(msg, index); break;
		case 18: leaf18(msg, index); break;
		case 19: leaf19(msg, index); break;
		case 20: leaf20(msg, index); break;
		case 21: leaf21(msg, index); break;
		case 22: leaf22(msg, index); break;
		case 23: leaf23(msg, index); break;
		case 24: leaf24(msg, index); break;
		case 25: leaf25(msg, index); break;
		case 26: leaf26(msg, index); break;
		case 27: leaf27(msg, index); break;
		case 28: leaf28(msg, index); break;
		case 29: leaf29(msg, index); break;
		case 30: leaf30(msg, index); break;
		case 31: leaf31(msg, index); break;
		case 32: leaf32(msg, index); break;
		case 33: leaf33(msg, index); break;
		case 34: leaf34(msg, index); break;
		case 35: leaf35(msg, index); break;
		case 36: leaf36(msg, index); break;
		case 37: leaf37(msg, index); break;
		case 38: leaf38(msg, index); break;
		case 39: leaf39(msg, index); break;
		case 40: leaf40(msg, index); break;
		case 41: leaf41(msg, index); break;
		case 42: leaf42(msg, index); break;
		case 43: leaf43(msg, index); break;
		case 44: leaf44(msg, index); break;
		case 45: leaf45(msg, index); break;
		case 46: leaf46(msg, index); break;
		case 47: leaf47(msg, index); break;
		case 48: leaf48(msg, index); break;
		case 49: leaf49(msg, index); break;
		case 50: leaf50(msg, index); break;
		case 51: leaf51(msg, index); break;
		case 52: leaf52(msg, index); break;
		case 53: leaf53(msg, index); break;
		case 54: leaf54(msg, index); break;
		case 55: leaf55(msg, index); break;
		case 56: leaf56(msg, index); break;
		case 57: leaf57(msg, index); break;
		case 58: leaf58(msg, index); break;
		case 59: leaf59(msg, index); break;
		case 60: leaf60(msg, index); break;
		case 61: leaf61(msg, index); break;
		case 62: leaf62(msg, index); break;
		case 63: leaf63(msg, index); break;
		case 64: leaf64(msg, index); break;
		case 65: leaf65(msg, index); break;
		case 66: leaf66(msg, index); break;
		case 67: leaf67(msg, index); break;
		case 68: leaf68(msg, index); break;
		case 69: leaf69(msg, index); break;
		case 70: leaf70(msg, index); break;
		case 71: leaf71(msg, index); break;
		case 72: leaf72(msg, index); break;
		case 73: leaf73(msg, index); break;
		case 74: leaf74(msg, index); break;
		case 75: leaf75(msg, index); break;
		case 76: leaf76(msg, index); break;
		case 77: leaf77(msg, index); break;
		case 78: leaf78(msg, index); break;
		case 79: leaf79(msg, index); break;
		case 80: leaf80(msg, index); break;
		case 81: leaf81(msg, index); break;
		case 82: leaf82(msg, index); break;
		case 83: leaf83(msg, index); break;
		case 84: leaf84(msg, index); break;
		case 85: leaf85(msg, index); break;
		case 86: leaf86(msg, index); break;
		case 87: leaf87(msg, index); break;
		case 88: leaf88(msg, index); break;
		case 89: leaf89(msg, index); break;
		case 90: leaf90(msg, index); break;
		case 91: leaf91(msg, index); break;
		case 92: leaf92(msg, index); break;
		case 93: leaf93(msg, index); break;
		case 94: leaf94(msg, index); break;
		case 95: leaf95(msg, index); break;
		case 96: leaf96(msg, index); break;
		case 97: leaf97(msg, index); break;
		case 98: leaf98(msg, index); break;
		case 99: leaf99(msg, index); break;
		case 100: leaf100(msg, index); break;
		case 101: leaf101(msg, index); break;
		case 102: leaf102(msg, index); break;
		case 103: leaf103(msg, index); break;
		case 104: leaf104(msg, index); break;
		case 105: leaf105(msg, index); break;
		case 106: leaf106(msg, index); break;
		case 107: leaf107(msg, index); break;
		case 108: leaf108(msg, index); break;
		case 109: leaf109(msg, index); break;
		case 110: leaf110(msg, index); break;
		case 111: leaf111(msg, index); break;
		case 112: leaf112(msg, index); break;
		case 113: leaf113(msg, index); break;
		case 114: leaf114(msg, index); break;
		case 115: leaf115(msg, index); break;
		case 116: leaf116(msg, index); break;
		case 117: leaf117(msg, index); break;
		case 118: leaf118(msg, index); break;
		case 119: leaf119(msg, index); break;
		case 120: leaf120(msg, index); break;
		case 121: leaf121(msg, index); break;
		case 122: leaf122(msg, index); break;
		case 123: leaf123(msg, index); break;
		case 124: leaf124(msg, index); break;
		case 125: leaf125(msg, index); break;
		case 126: leaf126(msg, index); break;
		case 127: leaf127(msg, index); break;
		case 128: leaf128(msg, index); break;
		case 129: leaf129(msg, index); break;
		case 130: leaf130(msg, index); break;
		case 131: leaf131(msg, index); break;
		case 132: leaf132(msg, index); break;
		case 133: leaf133(msg, index); break;
		case 134: leaf134(msg, index); break;
		case 135: leaf135(msg, index); break;
		case 136: leaf136(msg, index); break;
		case 137: leaf137(msg, index); break;
		case 138: leaf138(msg, index); break;
		case 139: leaf139(msg, index); break;
		case 140: leaf140(msg, index); break;
		case 141: leaf141(msg, index); break;
		case 142: leaf142(msg, index); break;
		case 143: leaf143(msg, index); break;
		case 144: leaf144(msg, index); break;
		case 145: leaf145(msg, index); break;
		case 146: leaf146(msg, index); break;
		case 147: leaf147(msg, index); break;
		case 148: leaf148(msg, index); break;
		case 149: leaf149(msg, index); break;
		case 150: leaf150(msg, index); break;
		case 151: leaf151(msg, index); break;
		case 152: leaf152(msg, index); break;
		case 153: leaf153(msg, index); break;
		case 154: leaf154(msg, index); break;
		case 155: leaf155(msg, index); break;
		case 156: leaf156(msg, index); break;
		case 157: leaf157(msg, index); break;
		case 158: leaf158(msg, index); break;
		case 159: leaf159(msg, index); break;
		case 160: leaf160(msg, index); break;
		case 161: leaf161(msg, index); break;
		case 162: leaf162(msg, index); break;
		case 163: leaf163(msg, index); break;
		case 164: leaf164(msg, index); break;
		case 165: leaf165(msg, index); break;
		case 166: leaf166(msg, index); break;
		case 167: leaf167(msg, index); break;
		case 168: leaf168(msg, index); break;
		case 169: leaf169(msg, index); break;
		case 170: leaf170(msg, index); break;
		case 171: leaf171(msg, index); break;
		case 172: leaf172(msg, index); break;
		case 173: leaf173(msg, index); break;
		case 174: leaf174(msg, index); break;
		case 175: leaf175(msg, index); break;
		case 176: leaf176(msg, index); break;
		case 177: leaf177(msg, index); break;
		case 178: leaf178(msg, index); break;
		case 179: leaf179(msg, index); break;
		case 180: leaf180(msg, index); break;
		case 181: leaf181(msg, index); break;
		case 182: leaf182(msg, index); break;
		case 183: leaf183(msg, index); break;
		case 184: leaf184(msg, index); break;
		case 185: leaf185(msg, index); break;
		case 186: leaf186(msg, index); break;
		case 187: leaf187(msg, index); break;
		case 188: leaf188(msg, index); break;
		case 189: leaf189(msg, index); break;
		case 190: leaf190(msg, index); break;
		case 191: leaf191(msg, index); break;
		case 192: leaf192(msg, index); break;
		case 193: leaf193(msg, index); break;
		case 194: leaf194(msg, index); break;
		case 195: leaf195(msg, index); break;
		case 196: leaf196(msg, index); break;
		case 197: leaf197(msg, index); break;
		case 198: leaf198(msg, index); break;
		case 199: leaf199(msg, index); break;
		case 200: leaf200(msg, index); break;
		case 201: leaf201(msg, index); break;
		case 202: leaf202(msg, index); break;
		case 203: leaf203(msg, index); break;
		case 204: leaf204(msg, index); break;
		case 205: leaf205(msg, index); break;
		case 206: leaf206(msg, index); break;
		case 207: leaf207(msg, index); break;
		case 208: leaf208(msg, index); break;
		case 209: leaf209(msg, index); break;
		case 210: leaf210(msg, index); break;
		case 211: leaf211(msg, index); break;
		case 212: leaf212(msg, index); break;
		case 213: leaf213(msg, index); break;
		case 214: leaf214(msg, index); break;
		case 215: leaf215(msg, index); break;
		case 216: leaf216(msg, index); break;
		case 217: leaf217(msg, index); break;
		case 218: leaf218(msg, index); break;
		case 219: leaf219(msg, index); break;
		case 220: leaf220(msg, index); break;
		case 221: leaf221(msg, index); break;
		case 222: leaf222(msg, index); break;
		case 223: leaf223(msg, index); break;
		case 224: leaf224(msg, index); break;
		case 225: leaf225(msg, index); break;
		case 226: leaf226(msg, index); break;
		case 227: leaf227(msg, index); break;
		case 228: leaf228(msg, index); break;
		case 229: leaf229(msg, index); break;
		case 230: leaf230(msg, index); break;
		case 231: leaf231(msg, index); break;
		case 232: leaf232(msg, index); break;
		case 233: leaf233(msg, index); break;
		case 234: leaf234(msg, index); break;
		case 235: leaf235(msg, index); break;
		case 236: leaf236(msg, index); break;
		case 237: leaf237(msg, index); break;
		case 238: leaf238(msg, index); break;
		case 239: leaf239(msg, index); break;
		case 240: leaf240(msg, index); break;
		case 241: leaf241(msg, index); break;
		case 242: leaf242(msg, index); break;
		case 243: leaf243(msg, index); break;
		case 244: leaf244(msg, index); break;
		case 245: leaf245(msg, index); break;
		case 246: leaf246(msg, index); break;
		case 247: leaf247(msg, index); break;
		case 248: leaf248(msg, index); break;
		case 249: leaf249(msg, index); break;
		case 250: leaf250(msg, index); break;
		case 251: leaf251(msg, index); break;
		case 252: leaf252(msg, index); break;
		case 253: leaf253(msg, index); break;
		case 254: leaf254(msg, index); break;
		case 255: leaf255(msg, index); break;
		case 256+(0<<5)+4: node0<4>(msg, index); break;
		case 256+(0<<5)+5: node0<5>(msg, index); break;
		case 256+(0<<5)+6: node0<6>(msg, index); break;
		case 256+(0<<5)+7: node0<7>(msg, index); break;
		case 256+(0<<5)+8: node0<8>(msg, index); break;
		case 256+(0<<5)+9: node0<9>(msg, index); break;
		case 256+(0<<5)+10: node0<10>(msg, index); break;
		case 256+(0<<5)+11: node0<11>(msg, index); break;
		case 256+(0<<5)+12: node0<12>(msg, index); break;
		case 256+(0<<5)+13: node0<13>(msg, index); break;
		case 256+(0<<5)+14: node0<14>(msg, index); break;
		case 256+(0<<5)+15: node0<15>(msg, index); break;
		case 256+(0<<5)+16: node0<16>(msg, index); break;
		case 256+(0<<5)+17: node0<17>(msg, index); break;
		case 256+(0<<5)+18: node0<18>(msg, index); break;
		case 256+(0<<5)+19: node0<19>(msg, index); break;
		case 256+(0<<5)+20: node0<20>(msg, index); break;
		case 256+(0<<5)+21: node0<21>(msg, index); break;
		case 256+(0<<5)+22: node0<22>(msg, index); break;
		case 256+(0<<5)+23: node0<23>(msg, index); break;
		case 256+(0<<5)+24: node0<24>(msg, index); break;
		case 256+(0<<5)+25: node0<25>(msg, index); break;
		case 256+(0<<5)+26: node0<26>(msg, index); break;
		case 256+(0<<5)+27: node0<27>(msg, index); break;
		case 256+(0<<5)+28: node0<28>(msg, index); break;
		case 256+(0<<5)+29: node0<29>(msg, index); break;
		case 256+(0<<5)+30: node0<30>(msg, index); break;
		case 256+(0<<5)+31: node0<31>(msg, index); break;
		case 256+(1<<5)+4: node1<4>(msg, index); break;
		case 256+(1<<5)+5: node1<5>(msg, index); break;
		case 256+(1<<5)+6: node1<6>(msg, index); break;
		case 256+(1<<5)+7: node1<7>(msg, index); break;
		case 256+(1<<5)+8: node1<8>(msg, index); break;
		case 256+(1<<5)+9: node1<9>(msg, index); break;
		case 256+(1<<5)+10: node1<10>(msg, index); break;
		case 256+(1<<5)+11: node1<11>(msg, index); break;
		case 256+(1<<5)+12: node1<12>(msg, index); break;
		case 256+(1<<5)+13: node1<13>(msg, index); break;
		case 256+(1<<5)+14: node1<14>(msg, index); break;
		case 256+(1<<5)+15: node1<15>(msg, index); break;
		case 256+(1<<5)+16: node1<16>(msg, index); break;
		case 256+(1<<5)+17: node1<17>(msg, index); break;
		case 256+(1<<5)+18: node1<18>(msg, index); break;
		case 256+(1<<5)+19: node1<19>(msg, index); break;
		case 256+(1<<5)+20: node1<20>(msg, index); break;
		case 256+(1<<5)+21: node1<21>(msg, index); break;
		case 256+(1<<5)+22: node1<22>(msg, index); break;
		case 256+(1<<5)+23: node1<23>(msg, index); break;
		case 256+(1<<5)+24: node1<24>(msg, index); break;
		case 256+(1<<5)+25: node1<25>(msg, index); break;
		case 256+(1<<5)+26: node1<26>(msg, index); break;
		case 256+(1<<5)+27: node1<27>(msg, index); break;
		case 256+(1<<5)+28: node1<28>(msg, index); break;
		case 256+(1<<5)+29: node1<29>(msg, index); break;
		case 256+(1<<5)+30: node1<30>(msg, index); break;
		case 256+(1<<5)+31: node1<31>(msg, index); break;
		case 256+(2<<5)+4: node2<4>(msg, index); break;
		case 256+(2<<5)+5: node2<5>(msg, index); break;
		case 256+(2<<5)+6: node2<6>(msg, index); break;
		case 256+(2<<5)+7: node2<7>(msg, index); break;
		case 256+(2<<5)+8: node2<8>(msg, index); break;
		case 256+(2<<5)+9: node2<9>(msg, index); break;
		case 256+(2<<5)+10: node2<10>(msg, index); break;
		case 256+(2<<5)+11: node2<11>(msg, index); break;
		case 256+(2<<5)+12: node2<12>(msg, index); break;
		case 256+(2<<5)+13: node2<13>(msg, index); break;
		case 256+(2<<5)+14: node2<14>(msg, index); break;
		case 256+(2<<5)+15: node2<15>(msg, index); break;
		case 256+(2<<5)+16: node2<16>(msg, index); break;
		case 256+(2<<5)+17: node2<17>(msg, index); break;
		case 256+(2<<5)+18: node2<18>(msg, index); break;
		case 256+(2<<5)+19: node2<19>(msg, index); break;
		case 256+(2<<5)+20: node2<20>(msg, index); break;
		case 256+(2<<5)+21: node2<21>(msg, index); break;
		case 256+(2<<5)+22: node2<22>(msg, index); break;
		case 256+(2<<5)+23: node2<23>(msg, index); break;
		case 256+(2<<5)+24: node2<24>(msg, index); break;
		case 256+(2<<5)+25: node2<25>(msg, index); break;
		case 256+(2<<5)+26: node2<26>(msg, index); break;
		case 256+(2<<5)+27: node2<27>(msg, index); break;
		case 256+(2<<5)+28: node2<28>(msg, index); break;
		case 256+(2<<5)+29: node2<29>(msg, index); break;
		case 256+(2<<5)+30: node2<30>(msg, index); break;
		case 256+(2<<5)+31: node2<31>(msg, index); break;
		case 256+(3<<5)+4: rate0<4>(msg, index); break;
		case 256+(3<<5)+5: rate0<5>(msg, index); break;
		case 256+(3<<5)+6: rate0<6>(msg, index); break;
		case 256+(3<<5)+7: rate0<7>(msg, index); break;
		case 256+(3<<5)+8: rate0<8>(msg, index); break;
		case 256+(3<<5)+9: rate0<9>(msg, index); break;
		case 256+(3<<5)+10: rate0<10>(msg, index); break;
		case 256+(3<<5)+11: rate0<11>(msg, index); break;
		case 256+(3<<5)+12: rate0<12>(msg, index); break;
		case 256+(3<<5)+13: rate0<13>(msg, index); break;
		case 256+(3<<5)+14: rate0<14>(msg, index); break;
		case 256+(3<<5)+15: rate0<15>(msg, index); break;
		case 256+(3<<5)+16: rate0<16>(msg, index); break;
		case 256+(3<<5)+17: rate0<17>(msg, index); break;
		case 256+(3<<5)+18: rate0<18>(msg, index); break;
		case 256+(3<<5)+19: rate0<19>(msg, index); break;
		case 256+(3<<5)+20: rate0<20>(msg, index); break;
		case 256+(3<<5)+21: rate0<21>(msg, index); break;
		case 256+(3<<5)+22: rate0<22>(msg, index); break;
		case 256+(3<<5)+23: rate0<23>(msg, index); break;
		case 256+(3<<5)+24: rate0<24>(msg, index); break;
		case 256+(3<<5)+25: rate0<25>(msg, index); break;
		case 256+(3<<5)+26: rate0<26>(msg, index); break;
		case 256+(3<<5)+27: rate0<27>(msg, index); break;
		case 256+(3<<5)+28: rate0<28>(msg, index); break;
		case 256+(3<<5)+29: rate0<29>(msg, index); break;
		case 256+(3<<5)+30: rate0<30>(msg, index); break;
		case 256+(3<<5)+31: rate0<31>(msg, index); break;
		case 256+(4<<5)+4: rate1<4>(msg, index); break;
		case 256+(4<<5)+5: rate1<5>(msg, index); break;
		case 256+(4<<5)+6: rate1<6>(msg, index); break;
		case 256+(4<<5)+7: rate1<7>(msg, index); break;
		case 256+(4<<5)+8: rate1<8>(msg, index); break;
		case 256+(4<<5)+9: rate1<9>(msg, index); break;
		case 256+(4<<5)+10: rate1<10>(msg, index); break;
		case 256+(4<<5)+11: rate1<11>(msg, index); break;
		case 256+(4<<5)+12: rate1<12>(msg, index); break;
		case 256+(4<<5)+13: rate1<13>(msg, index); break;
		case 256+(4<<5)+14: rate1<14>(msg, index); break;
		case 256+(4<<5)+15: rate1<15>(msg, index); break;
		case 256+(4<<5)+16: rate1<16>(msg, index); break;
		case 256+(4<<5)+17: rate1<17>(msg, index); break;
		case 256+(4<<5)+18: rate1<18>(msg, index); break;
		case 256+(4<<5)+19: rate1<19>(msg, index); break;
		case 256+(4<<5)+20: rate1<20>(msg, index); break;
		case 256+(4<<5)+21: rate1<21>(msg, index); break;
		case 256+(4<<5)+22: rate1<22>(msg, index); break;
		case 256+(4<<5)+23: rate1<23>(msg, index); break;
		case 256+(4<<5)+24: rate1<24>(msg, index); break;
		case 256+(4<<5)+25: rate1<25>(msg, index); break;
		case 256+(4<<5)+26: rate1<26>(msg, index); break;
		case 256+(4<<5)+27: rate1<27>(msg, index); break;
		case 256+(4<<5)+28: rate1<28>(msg, index); break;
		case 256+(4<<5)+29: rate1<29>(msg, index); break;
		case 256+(4<<5)+30: rate1<30>(msg, index); break;
		case 256+(4<<5)+31: rate1<31>(msg, index); break;
		default:
			assert(false);
		}
	}
public:
	void operator()(int8_t *message, const int8_t *codeword, const uint32_t *program)
	{
		for (int i = 0; i < N; ++i)
			soft[i+N] = codeword[i];
		while (*program != 0xffffffff) {
			int func = *program >> 23;
			int index = (*program & 0x007fffff) << U;
			decode(&message, func, index);
			++program;
		};
	}
};

int main()
{
	const int M = 20;
	const int N = 1 << M;
	const int U = 3; // unrolled at level 3
	const bool systematic = true;
	std::random_device rd;
	typedef std::default_random_engine generator;
	typedef std::uniform_int_distribution<int> distribution;
	auto data = std::bind(distribution(0, 1), generator(rd()));
	auto frozen = new uint8_t[N>>U];
	auto codeword = new int8_t[N];
	PolarFreezer<M> freeze;
	double erasure_probability = 0.5;
	double freezing_threshold = 0.5;
	int K = freeze(frozen, erasure_probability, freezing_threshold);
	std::cerr << "Polar(" << N << ", " << K << ")" << std::endl;
	auto message = new int8_t[K];
	auto decoded = new int8_t[K];
	PolarEncoder<M> encode;
	auto program = new uint32_t[N];
	PolarCompiler<M> compile;
	compile(program, frozen);
	std::cerr << "sizeof(PolarDecoder<M>) = " << sizeof(PolarDecoder<M>) << std::endl;
	auto decode = new PolarDecoder<M>;

#if 0
	auto epos = std::bind(distribution(0, N-1), generator(rd()));
	std::cerr << "errors erasures msec" << std::endl;
	for (int loop = 0; loop < 100; ++loop) {
		for (int i = 0; i < K; ++i)
			message[i] = 1 - 2 * data();
		if (systematic) {
			PolarSysEnc<M> sysenc;
			sysenc(codeword, message, frozen);
			for (int i = 0, j = 0; i < N; ++i)
				if (!(frozen[i>>U]&(1<<(i&((1<<U)-1)))))
					assert(codeword[i] == message[j++]);
		} else {
			encode(codeword, message, frozen);
		}
		for (int i = 0; i < erasure_probability * N; ++i)
			codeword[epos()] = 0;
		auto start = std::chrono::system_clock::now();
		(*decode)(decoded, codeword, program);
		auto end = std::chrono::system_clock::now();
		auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		if (systematic) {
			encode(codeword, decoded, frozen);
			for (int i = 0, j = 0; i < N; ++i)
				if (!(frozen[i>>U]&(1<<(i&((1<<U)-1)))))
					decoded[j++] = codeword[i];
		}
		int erasures = 0;
		for (int i = 0; i < K; ++i)
			erasures += !decoded[i];
		int errors = 0;
		for (int i = 0; i < K; ++i)
			errors += decoded[i] * message[i] < 0;
		std::cout << errors << " " << erasures << " " << msec.count() << std::endl;
	}
#else
	auto orig = new int8_t[N];
	auto noisy = new int8_t[N];
	auto symb = new double[N];
	double min_SNR = 20, max_mbs = 0;
	for (double SNR = -5; SNR <= 0; SNR += 0.1) {
		//double mean_signal = 0;
		double sigma_signal = 1;
		double mean_noise = 0;
		double sigma_noise = std::sqrt(sigma_signal * sigma_signal / (2 * std::pow(10, SNR / 10)));

		typedef std::normal_distribution<double> normal;
		auto awgn = std::bind(normal(mean_noise, sigma_noise), generator(rd()));

		int64_t awgn_errors = 0;
		int64_t quantization_erasures = 0;
		int64_t uncorrected_errors = 0;
		int64_t ambiguity_erasures = 0;
		int64_t loops = 10;
		double avg_mbs = 0;
		for (int l = 0; l < loops; ++l) {
			for (int i = 0; i < K; ++i)
				message[i] = 1 - 2 * data();

			if (systematic) {
				PolarSysEnc<M> sysenc;
				sysenc(codeword, message, frozen);
				for (int i = 0, j = 0; i < N; ++i)
					if (!(frozen[i>>U]&(1<<(i&((1<<U)-1)))))
						assert(codeword[i] == message[j++]);
			} else {
				encode(codeword, message, frozen);
			}

			for (int i = 0; i < N; ++i)
				orig[i] = codeword[i];

			for (int i = 0; i < N; ++i)
				symb[i] = codeword[i];

			for (int i = 0; i < N; ++i)
				symb[i] += awgn();

			// $LLR=log(\frac{p(x=+1|y)}{p(x=-1|y)})$
			// $p(x|\mu,\sigma)=\frac{1}{\sqrt{2\pi}\sigma}}e^{-\frac{(x-\mu)^2}{2\sigma^2}}$
			double DIST = 2; // BPSK
			double fact = DIST / (sigma_noise * sigma_noise);
			for (int i = 0; i < N; ++i)
				codeword[i] = std::min<double>(std::max<double>(std::nearbyint(fact * symb[i]), -128), 127);

			for (int i = 0; i < N; ++i)
				noisy[i] = codeword[i];

			auto start = std::chrono::system_clock::now();
			(*decode)(decoded, codeword, program);
			auto end = std::chrono::system_clock::now();
			auto usec = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
			double mbs = (double)K / usec.count();
			avg_mbs += mbs;

			if (systematic) {
				encode(codeword, decoded, frozen);
				for (int i = 0, j = 0; i < N; ++i)
					if (!(frozen[i>>U]&(1<<(i&((1<<U)-1)))))
						decoded[j++] = codeword[i];
			}

			for (int i = 0; i < N; ++i)
				awgn_errors += noisy[i] * orig[i] < 0;
			for (int i = 0; i < N; ++i)
				quantization_erasures += !noisy[i];
			for (int i = 0; i < K; ++i)
				uncorrected_errors += decoded[i] * message[i] < 0;
			for (int i = 0; i < K; ++i)
				ambiguity_erasures += !decoded[i];
		}

		avg_mbs /= loops;
		max_mbs = std::max(max_mbs, avg_mbs);
		double bit_error_rate = (double)(uncorrected_errors + ambiguity_erasures) / (double)(K * loops);
		if (!uncorrected_errors && !ambiguity_erasures)
			min_SNR = std::min(min_SNR, SNR);

		if (0) {
			std::cerr << SNR << " Es/N0 => AWGN with standard deviation of " << sigma_noise << " and mean " << mean_noise << std::endl;
			std::cerr << awgn_errors << " errors caused by AWGN." << std::endl;
			std::cerr << quantization_erasures << " erasures caused by quantization." << std::endl;
			std::cerr << uncorrected_errors << " errors uncorrected." << std::endl;
			std::cerr << ambiguity_erasures << " ambiguity erasures." << std::endl;
			std::cerr << bit_error_rate << " bit error rate." << std::endl;
			std::cerr << avg_mbs << " megabit per second." << std::endl;
		} else {
			std::cout << SNR << " " << bit_error_rate << " " << avg_mbs << std::endl;
		}
	}
	std::cerr << "QEF at: " << min_SNR << " SNR, speed: " << max_mbs << " Mb/s." << std::endl;
#endif
	return 0;
}
