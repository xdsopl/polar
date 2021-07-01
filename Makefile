
CXXFLAGS = -std=c++17 -W -Wall -Ofast -fno-exceptions -fno-rtti
CXX = clang++ -stdlib=libc++ -march=native
#CXX = g++ -march=native

#CXX = armv7a-hardfloat-linux-gnueabi-g++ -static -mfpu=neon -march=armv7-a
#QEMU = qemu-arm

.PHONY: all

all: testbench

.PHONY: test

test: testbench
	$(QEMU) ./testbench

.PHONY: emit

emit: emitter
	$(QEMU) ./emitter > program.hh

emitter: polar_emitter.cc *.hh
	$(CXX) $(CXXFLAGS) $< -o $@

testbench: testbench.cc *.hh
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -f testbench emitter program.hh

