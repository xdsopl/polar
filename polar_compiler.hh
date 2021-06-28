/*
Compiler for successive cancellation decoding of polar codes

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#pragma once

class PolarCompiler
{
	static const int left = 0, right = 1, comb = 2,
		rate0 = 3, rate1 = 4, rep = 5, spc = 6;
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
		int count = frozen_count(frozen, level);
		if (count == 1<<level) {
			*(*program)++ = rate0;
		} else if (count == 0) {
			*(*program)++ = rate1;
		} else if (count == (1<<level)-1 && !frozen[(1<<level)-1]) {
			*(*program)++ = rep;
		} else if (count == 1 && frozen[0]) {
			*(*program)++ = spc;
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

