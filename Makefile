
CXXFLAGS = -std=c++17 -W -Wall -O2 -fno-exceptions -fno-rtti -ffast-math -ftree-vectorize
CXX = clang++ -stdlib=libc++ -march=native
#CXX = g++ -march=native

#CXX = armv7a-hardfloat-linux-gnueabi-g++ -static -mfpu=neon -march=armv7-a
#QEMU = qemu-arm

#CXX = aarch64-unknown-linux-gnu-g++ -static -march=armv8-a+crc+simd -mtune=cortex-a72
#QEMU = qemu-aarch64

.PHONY: all

all: testbench

.PHONY: test

test: testbench
	$(QEMU) ./testbench

testbench: testbench.cc *.hh
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -f testbench

