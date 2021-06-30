/*
Test bench for successive cancellation decoding of polar codes

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
#include "polar_compiler.hh"
#include "polar_decoder.hh"
#include "polar_encoder.hh"
#include "polar_freezer.hh"

int main()
{
	const int M = 16;
	const int N = 1 << M;
	const bool systematic = true;

	int loops = 10;
	std::random_device rd;
	typedef std::default_random_engine generator;
	typedef std::uniform_int_distribution<int> distribution;
	auto data = std::bind(distribution(0, 1), generator(rd()));
	auto frozen = new uint8_t[N];
	auto codeword = new int8_t[N];

	long double erasure_probability = 1. / 3.;
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
	auto message = new int8_t[K];
	auto decoded = new int8_t[K];
	PolarEncoder<M> encode;
	auto program = new uint8_t[N];
	PolarCompiler compile;
	int length = compile(program, frozen, M);
	std::cerr << "program length = " << length << std::endl;
	std::cerr << "sizeof(PolarDecoder<M>) = " << sizeof(PolarDecoder<M>) << std::endl;
	auto decode = new PolarDecoder<M>();

	auto orig = new int8_t[N];
	auto noisy = new int8_t[N];
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
				if (1) {
					PolarSysEnc<M> sysenc;
					sysenc(codeword, message, frozen);
				} else {
					for (int i = 0, j = 0; i < N; ++i)
						codeword[i] = frozen[i] ? 0 : message[j++];
					(*decode)(decoded, codeword, program);
					encode(codeword, decoded, frozen);
				}
				for (int i = 0, j = 0; i < N; ++i)
					if (!frozen[i])
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
				codeword[i] = std::clamp((int)std::nearbyint(fact * symb[i]), -128, 127);

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
					if (!frozen[i])
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
