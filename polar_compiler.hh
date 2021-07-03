/*
Compiler for successive cancellation decoding of polar codes

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#pragma once

class PolarCompiler
{
	static const int left = 0, right = 1, comb = 2,
		rate0 = 3, rate1 = 4, rep = 5, spc = 6,
		rate0_right = 7, rate0_comb = 8, rate1_comb = 9;
	static int frozen_count(const uint8_t *frozen, int level)
	{
		int count = 0;
		for (int i = 0; i < (1<<level); ++i)
			count += frozen[i];
		return count;
	}
	static void compile(uint8_t **program, const uint8_t *frozen, int level)
	{
		assert(level > 0);
		int lcnt = frozen_count(frozen, level-1);
		int rcnt = frozen_count(frozen+(1<<(level-1)), level-1);
		if (lcnt == 1<<(level-1) && rcnt == 1<<(level-1)) {
			*(*program)++ = rate0;
		} else if (lcnt == 0 && rcnt == 0) {
			*(*program)++ = rate1;
		} else if (lcnt == 1<<(level-1) && rcnt == (1<<(level-1))-1 && !frozen[(1<<level)-1]) {
			*(*program)++ = rep;
		} else if (lcnt == 1 && rcnt == 0 && frozen[0]) {
			*(*program)++ = spc;
		} else if (lcnt == 1<<(level-1)) {
			*(*program)++ = rate0_right;
			compile(program, frozen+(1<<(level-1)), level-1);
			*(*program)++ = rate0_comb;
		} else if (rcnt == 0) {
			*(*program)++ = left;
			compile(program, frozen, level-1);
			*(*program)++ = rate1_comb;
		} else {
			*(*program)++ = left;
			compile(program, frozen, level-1);
			*(*program)++ = right;
			compile(program, frozen+(1<<(level-1)), level-1);
			*(*program)++ = comb;
		}
	}
public:
	int operator()(uint8_t *program, const uint8_t *frozen, int level)
	{
		uint8_t *first = program;
		*program++ = level;
		compile(&program, frozen, level);
		*program++ = 255;
		return program - first;
	}
};

