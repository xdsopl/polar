/*
Test bench for successive cancellation decoding of polar codes

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#include <limits>
#include <stdlib.h>
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

int main()
{
	const int M = 14;
	const int N = 1 << M;
	const bool systematic = true;
#if 1
	typedef int8_t code_type;
#else
	typedef float code_type;
#endif

#if 0
	const int SIMD_WIDTH = 1;
	typedef code_type simd_type;
#else
#ifdef __AVX2__
	const int SIZEOF_SIMD = 32;
#else
	const int SIZEOF_SIMD = 16;
#endif
	const int SIMD_WIDTH = SIZEOF_SIMD / sizeof(code_type);
	typedef SIMD<code_type, SIMD_WIDTH> simd_type;
#endif
	std::random_device rd;
	typedef std::default_random_engine generator;
	typedef std::uniform_int_distribution<int> distribution;
	auto data = std::bind(distribution(0, 1), generator(rd()));
	auto frozen = new uint8_t[N];
	auto codeword = reinterpret_cast<code_type *>(aligned_alloc(sizeof(simd_type), sizeof(simd_type) * N));

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
	auto message = reinterpret_cast<code_type *>(aligned_alloc(sizeof(simd_type), sizeof(simd_type) * K));
	auto decoded = reinterpret_cast<code_type *>(aligned_alloc(sizeof(simd_type), sizeof(simd_type) * K));
	PolarEncoder<simd_type, M> encode;
	auto program = new uint8_t[N];
	PolarCompiler compile;
	int length = compile(program, frozen, M);
	std::cerr << "program length = " << length << std::endl;
	std::cerr << "sizeof(PolarDecoder<simd_type, M>) = " << sizeof(PolarDecoder<simd_type, M>) << std::endl;
	auto decode = reinterpret_cast<PolarDecoder<simd_type, M> *>(aligned_alloc(sizeof(simd_type), sizeof(PolarDecoder<simd_type, M>)));

	auto orig = reinterpret_cast<code_type *>(aligned_alloc(sizeof(simd_type), sizeof(simd_type) * N));
	auto noisy = reinterpret_cast<code_type *>(aligned_alloc(sizeof(simd_type), sizeof(simd_type) * N));
	auto symb = new double[SIMD_WIDTH*N];
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
		int64_t loops = 0;
		while (uncorrected_errors < 1000 && ++loops < 320 / SIMD_WIDTH) {
			for (int i = 0; i < SIMD_WIDTH * K; ++i)
				message[i] = 1 - 2 * data();

			if (systematic) {
				if (1) {
					PolarSysEnc<simd_type, M> sysenc;
					sysenc(reinterpret_cast<simd_type *>(codeword), reinterpret_cast<simd_type *>(message), frozen);
				} else {
					for (int i = 0, j = 0; i < N; ++i)
						for (int k = 0; k < SIMD_WIDTH; ++k)
							if (frozen[i])
								codeword[SIMD_WIDTH*i+k] = 0;
							else
								codeword[SIMD_WIDTH*i+k] = message[j++];
					(*decode)(reinterpret_cast<simd_type *>(decoded), reinterpret_cast<simd_type *>(codeword), program);
					encode(reinterpret_cast<simd_type *>(codeword), reinterpret_cast<simd_type *>(decoded), frozen);
				}
				for (int i = 0, j = 0; i < N; ++i)
					for (int k = 0; k < SIMD_WIDTH; ++k)
						if (!frozen[i])
							assert(codeword[SIMD_WIDTH*i+k] == message[j++]);
			} else {
				encode(reinterpret_cast<simd_type *>(codeword), reinterpret_cast<simd_type *>(message), frozen);
			}

			for (int i = 0; i < SIMD_WIDTH * N; ++i)
				orig[i] = codeword[i];

			for (int i = 0; i < SIMD_WIDTH * N; ++i)
				symb[i] = codeword[i];

			for (int i = 0; i < SIMD_WIDTH * N; ++i)
				symb[i] += awgn();

			// $LLR=log(\frac{p(x=+1|y)}{p(x=-1|y)})$
			// $p(x|\mu,\sigma)=\frac{1}{\sqrt{2\pi}\sigma}}e^{-\frac{(x-\mu)^2}{2\sigma^2}}$
			double DIST = 2; // BPSK
			double fact = DIST / (sigma_noise * sigma_noise);
			for (int i = 0; i < SIMD_WIDTH * N; ++i)
				codeword[i] = PolarHelper<code_type>::quant(fact * symb[i]);

			for (int i = 0; i < SIMD_WIDTH * N; ++i)
				noisy[i] = codeword[i];

			auto start = std::chrono::system_clock::now();
			(*decode)(reinterpret_cast<simd_type *>(decoded), reinterpret_cast<simd_type *>(codeword), program);
			auto end = std::chrono::system_clock::now();
			auto usec = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
			double mbs = (double)(SIMD_WIDTH * K) / usec.count();
			avg_mbs += mbs;

			if (systematic) {
				encode(reinterpret_cast<simd_type *>(codeword), reinterpret_cast<simd_type *>(decoded), frozen);
				for (int i = 0, j = 0; i < N; ++i)
					for (int k = 0; k < SIMD_WIDTH; ++k)
						if (!frozen[i])
							decoded[j++] = codeword[SIMD_WIDTH*i+k];
			}

			for (int i = 0; i < SIMD_WIDTH * N; ++i)
				awgn_errors += noisy[i] * orig[i] < 0;
			for (int i = 0; i < SIMD_WIDTH * N; ++i)
				quantization_erasures += !noisy[i];
			for (int i = 0; i < SIMD_WIDTH * K; ++i)
				uncorrected_errors += decoded[i] * message[i] <= 0;
			for (int i = 0; i < SIMD_WIDTH * K; ++i)
				ambiguity_erasures += !decoded[i];
		}

		avg_mbs /= loops;
		max_mbs = std::max(max_mbs, avg_mbs);
		double bit_error_rate = (double)uncorrected_errors / (double)(SIMD_WIDTH * K * loops);
		if (!uncorrected_errors)
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
