/*
Test bench for SSCL decoding of polar codes

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#include <limits>
#include <random>
#include <chrono>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <functional>
#include "simd.hh"
#include "polar_helper.hh"
#include "polar_compiler.hh"
#include "polar_decoder.hh"
#include "polar_encoder.hh"
#include "polar_freezer.hh"

template <typename TYPE, int M>
class PolarTransform
{
	static const int N = 1 << M;
	typedef PolarHelper<TYPE> PH;
public:
	void operator()(TYPE *output, const TYPE *input)
	{
		for (int i = 0; i < N; i += 2) {
			TYPE in0 = input[i];
			TYPE in1 = input[i+1];
			output[i] = PH::qmul(in0, in1);
			output[i+1] = in1;
		}
		for (int h = 2; h < N; h *= 2)
			for (int i = 0; i < N; i += 2 * h)
				for (int j = i; j < i + h; ++j)
					output[j] = PH::qmul(output[j], output[j+h]);
	}
};

class PolarHistogram
{
	int hist[256];
public:
	void operator()(const uint8_t *program)
	{
		for (int i = 0; i < 256; ++i)
			hist[i] = 0;
		++program;
		while (*program != 255)
			++hist[*program++];
		int top[32] = { 0 };
		for (int j = 0; j < 6; ++j)
			for (int i = 0; i < 32; ++i)
				top[i] = std::max(top[i], hist[(j<<5)+i]);
		int len[32];
		for (int i = 0; i < 32; ++i)
			len[i] = 2 + std::log10(top[i]+!top[i]);
		int N = 0;
		for (int j = 0; j < 6; ++j)
			for (int i = N; i < 32; ++i)
				if (hist[(j<<5)+i])
					N = i;
		std::cerr << std::endl << "left: ";
		for (int i = 0; i <= N; ++i)
			std::cerr << std::setw(len[i]) << hist[(0<<5)+i];
		std::cerr << std::endl << "right:";
		for (int i = 0; i <= N; ++i)
			std::cerr << std::setw(len[i]) << hist[(1<<5)+i];
		std::cerr << std::endl << "comb: ";
		for (int i = 0; i <= N; ++i)
			std::cerr << std::setw(len[i]) << hist[(2<<5)+i];
		std::cerr << std::endl << "rate0:";
		for (int i = 0; i <= N; ++i)
			std::cerr << std::setw(len[i]) << hist[(3<<5)+i];
		std::cerr << std::endl << "rate1:";
		for (int i = 0; i <= N; ++i)
			std::cerr << std::setw(len[i]) << hist[(4<<5)+i];
		std::cerr << std::endl << "rep:  ";
		for (int i = 0; i <= N; ++i)
			std::cerr << std::setw(len[i]) << hist[(5<<5)+i];
		std::cerr << std::endl << std::endl;
	}
};

int main()
{
	const int M = 11;
	const int N = 1 << M;
	const bool systematic = true;
#if 1
	typedef int8_t code_type;
#else
	typedef float code_type;
#endif

#ifdef __AVX2__
	const int SIZEOF_SIMD = 32;
#else
	const int SIZEOF_SIMD = 16;
#endif
	const int SIMD_WIDTH = SIZEOF_SIMD / sizeof(code_type);
	typedef SIMD<code_type, SIMD_WIDTH> simd_type;

	int64_t loops = 320 / SIMD_WIDTH;
	std::random_device rd;
	typedef std::default_random_engine generator;
	typedef std::uniform_int_distribution<int> distribution;
	auto data = std::bind(distribution(0, 1), generator(rd()));
	auto frozen = new uint8_t[N];
	auto codeword = new code_type[N];
	auto temp = new simd_type[N];

	long double erasure_probability = 0.5;
	int K = (1 - erasure_probability) * N;
	double design_SNR = 10 * std::log10(-std::log(erasure_probability));
	std::cerr << "design SNR: " << design_SNR << std::endl;
	if (0) {
		PolarFreezer freeze;
		long double freezing_threshold = 0 ? 0.5 : std::numeric_limits<float>::epsilon();
		K = freeze(frozen, M, erasure_probability, freezing_threshold);
	} else {
		auto freeze = new PolarCodeConst0<M>;
		std::cerr << "sizeof(PolarCodeConst0<M>) = " << sizeof(PolarCodeConst0<M>) << std::endl;
		double better_SNR = design_SNR + 1.59175;
		std::cerr << "better SNR: " << better_SNR << std::endl;
		long double probability = std::exp(-pow(10.0, better_SNR / 10));
		(*freeze)(frozen, M, K, probability);
		delete freeze;
	}
	std::cerr << "Polar(" << N << ", " << K << ")" << std::endl;
	auto message = new code_type[K];
	auto decoded = new simd_type[K];
	PolarHelper<simd_type>::PATH metric[SIMD_WIDTH];
	auto program = new uint8_t[10*N];
	PolarCompiler compile;
	int length = compile(program, frozen, M);
	std::cerr << "program length = " << length << std::endl;
	std::cerr << "sizeof(PolarDecoder<simd_type, M>) = " << sizeof(PolarDecoder<simd_type, M>) << std::endl;
	if (1) {
		PolarHistogram histogram;
		histogram(program);
	}
	auto decode = new PolarDecoder<simd_type, M>;

	auto orig = new code_type[N];
	auto noisy = new code_type[N];
	auto symb = new double[N];
	double low_SNR = std::floor(design_SNR-3);
	double high_SNR = std::ceil(design_SNR+5);
	double min_SNR = high_SNR, max_mbs = 0;
	int count = 0;
	std::cerr << "SNR BER Mbit/s Eb/N0" << std::endl;
	for (double SNR = low_SNR; count <= 3 && SNR <= high_SNR; SNR += 0.1, ++count) {
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
		double avg_mbs = 0;
		for (int l = 0; l < loops; ++l) {
			for (int i = 0; i < K; ++i)
				message[i] = 1 - 2 * data();

			if (systematic) {
				PolarSysEnc<code_type, M> sysenc;
				sysenc(codeword, message, frozen);
				for (int i = 0, j = 0; i < N; ++i)
					if (!frozen[i])
						assert(codeword[i] == message[j++]);
			} else {
				PolarEncoder<code_type, M> encode;
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
				codeword[i] = PolarHelper<code_type>::quant(fact * symb[i]);

			for (int i = 0; i < N; ++i)
				noisy[i] = codeword[i];

			auto start = std::chrono::system_clock::now();
			(*decode)(metric, decoded, codeword, program);
			auto end = std::chrono::system_clock::now();
			auto usec = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
			double mbs = (double)K / usec.count();
			avg_mbs += mbs;

			if (systematic) {
				PolarEncoder<simd_type, M> encode;
				encode(temp, decoded, frozen);
				for (int i = 0, j = 0; i < N; ++i)
					if (!frozen[i])
						decoded[j++] = temp[i];
			}

			int best = 0;
			if (1) {
				for (int k = 0; k < SIMD_WIDTH; ++k)
					if (metric[k] < metric[best])
						best = k;
			} else {
				int errs[SIMD_WIDTH] = { 0 };
				for (int i = 0; i < K; ++i)
					for (int k = 0; k < SIMD_WIDTH; ++k)
						errs[k] += message[i] != decoded[i].v[k];
				for (int k = 0; k < SIMD_WIDTH; ++k)
					if (errs[k] < errs[best])
						best = k;
			}

			for (int i = 0; i < N; ++i)
				awgn_errors += noisy[i] * (orig[i] < 0);
			for (int i = 0; i < N; ++i)
				quantization_erasures += !noisy[i];
			for (int i = 0; i < K; ++i)
				uncorrected_errors += decoded[i].v[best] * message[i] < 0;
			for (int i = 0; i < K; ++i)
				ambiguity_erasures += !decoded[i].v[best];
		}

		avg_mbs /= loops;
		max_mbs = std::max(max_mbs, avg_mbs);
		double bit_error_rate = (double)(uncorrected_errors + ambiguity_erasures) / (double)(K * loops);
		if (!uncorrected_errors && !ambiguity_erasures)
			min_SNR = std::min(min_SNR, SNR);
		else
			count = 0;

		int MOD_BITS = 1; // BPSK
		double code_rate = (double)K / (double)N;
		double spectral_efficiency = code_rate * MOD_BITS;
		double EbN0 = 10 * std::log10(sigma_signal * sigma_signal / (spectral_efficiency * 2 * sigma_noise * sigma_noise));

		if (0) {
			std::cerr << SNR << " Es/N0 => AWGN with standard deviation of " << sigma_noise << " and mean " << mean_noise << std::endl;
			std::cerr << EbN0 << " Eb/N0, using spectral efficiency of " << spectral_efficiency << " from " << code_rate << " code rate and " << MOD_BITS << " bits per symbol." << std::endl;
			std::cerr << awgn_errors << " errors caused by AWGN." << std::endl;
			std::cerr << quantization_erasures << " erasures caused by quantization." << std::endl;
			std::cerr << uncorrected_errors << " errors uncorrected." << std::endl;
			std::cerr << ambiguity_erasures << " ambiguity erasures." << std::endl;
			std::cerr << bit_error_rate << " bit error rate." << std::endl;
			std::cerr << avg_mbs << " megabit per second." << std::endl;
		} else {
			std::cout << SNR << " " << bit_error_rate << " " << avg_mbs << " " << EbN0 << std::endl;
		}
	}
	std::cerr << "QEF at: " << min_SNR << " SNR, speed: " << max_mbs << " Mb/s." << std::endl;
	return 0;
}
