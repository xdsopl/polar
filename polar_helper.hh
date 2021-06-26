/*
SIMD wrapper used by polar encoder and decoder

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#pragma once

template <typename TYPE>
struct PolarHelper
{
	static TYPE one()
	{
		return 1;
	}
	static TYPE zero()
	{
		return 0;
	}
	static TYPE signum(TYPE v)
	{
		return (v > 0) - (v < 0);
	}
	static TYPE decide(TYPE v)
	{
		return (v >= 0) - (v < 0);
	}
	template <typename IN>
	static TYPE quant(IN in)
	{
		return in;
	}
	static TYPE qabs(TYPE a)
	{
		return std::abs(a);
	}
	static TYPE qmin(TYPE a, TYPE b)
	{
		return std::min(a, b);
	}
	static TYPE qadd(TYPE a, TYPE b)
	{
		return a + b;
	}
	static TYPE qmul(TYPE a, TYPE b)
	{
		return a * b;
	}
	static TYPE prod(TYPE a, TYPE b)
	{
		return signum(a) * signum(b) * qmin(qabs(a), qabs(b));
	}
	static TYPE madd(TYPE a, TYPE b, TYPE c)
	{
		return a * b + c;
	}
	static TYPE flip(TYPE a, TYPE b, TYPE c, TYPE d)
	{
		return c == d ? qmul(a, b) : a;
	}
};

template <typename VALUE, int WIDTH>
struct PolarHelper<SIMD<VALUE, WIDTH>>
{
	typedef SIMD<VALUE, WIDTH> TYPE;
	static TYPE one()
	{
		return vdup<TYPE>(1);
	}
	static TYPE zero()
	{
		return vzero<TYPE>();
	}
	static TYPE signum(TYPE a)
	{
		return vsignum(a);
	}
	static TYPE decide(TYPE a)
	{
		return vcopysign(one(), a);
	}
	static TYPE qabs(TYPE a)
	{
		return vabs(a);
	}
	static TYPE qmin(TYPE a, TYPE b)
	{
		return vmin(a, b);
	}
	static TYPE qadd(TYPE a, TYPE b)
	{
		return vadd(a, b);
	}
	static TYPE qmul(TYPE a, TYPE b)
	{
		return vmul(a, b);
	}
	static TYPE prod(TYPE a, TYPE b)
	{
		return vmul(vmul(vsignum(a), vsignum(b)), vmin(vabs(a), vabs(b)));
	}
	static TYPE madd(TYPE a, TYPE b, TYPE c)
	{
		return vadd(vmul(a, b), c);
	}
	static TYPE flip(TYPE a, TYPE b, TYPE c, TYPE d)
	{
		return vreinterpret<TYPE>(vbsl(vceq(c, d), vmask(qmul(a, b)), vmask(a)));
	}
};

template <int WIDTH>
struct PolarHelper<SIMD<int8_t, WIDTH>>
{
	typedef SIMD<int8_t, WIDTH> TYPE;
	static TYPE one()
	{
		return vdup<TYPE>(1);
	}
	static TYPE zero()
	{
		return vzero<TYPE>();
	}
	static TYPE signum(TYPE a)
	{
		return vsignum(a);
	}
	static TYPE decide(TYPE a)
	{
		return vreinterpret<TYPE>(vorr(vmask(one()), vcltz(a)));
	}
	static TYPE qabs(TYPE a)
	{
		return vqabs(a);
	}
	static TYPE qmin(TYPE a, TYPE b)
	{
		return vmin(a, b);
	}
	static TYPE qadd(TYPE a, TYPE b)
	{
		return vqadd(a, b);
	}
	static TYPE qmul(TYPE a, TYPE b)
	{
#ifdef __ARM_NEON__
		return vmul(a, b);
#else
		return vsign(a, b);
#endif
	}
	static TYPE prod(TYPE a, TYPE b)
	{
#ifdef __ARM_NEON__
		return vmul(vmul(vsignum(a), vsignum(b)), vmin(vqabs(a), vqabs(b)));
#else
		return vsign(vmin(vqabs(a), vqabs(b)), vsign(vsignum(a), b));
#endif
	}
	static TYPE madd(TYPE a, TYPE b, TYPE c)
	{
#ifdef __ARM_NEON__
		return vqadd(vmul(a, vmax(b, vdup<TYPE>(-127))), c);
#else
		return vqadd(vsign(vmax(b, vdup<TYPE>(-127)), a), c);
#endif
	}
	static TYPE flip(TYPE a, TYPE b, TYPE c, TYPE d)
	{
		return vreinterpret<TYPE>(vbsl(vceq(c, d), vmask(qmul(a, b)), vmask(a)));
	}
};

template <>
struct PolarHelper<int8_t>
{
	static int8_t one()
	{
		return 1;
	}
	static int8_t zero()
	{
		return 0;
	}
	static int8_t signum(int8_t v)
	{
		return (v > 0) - (v < 0);
	}
	static int8_t decide(int8_t v)
	{
		return (v >= 0) - (v < 0);
	}
	template <typename IN>
	static int8_t quant(IN in)
	{
		return std::min<IN>(std::max<IN>(std::nearbyint(in), -128), 127);
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
	static int8_t flip(int8_t a, int8_t b, int8_t c, int8_t d)
	{
		return c == d ? qmul(a, b) : a;
	}
};

