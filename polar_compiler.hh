/*
Compiler for successive cancellation decoding of polar codes

Copyright 2020 Ahmet Inan <xdsopl@gmail.com>
*/

#pragma once

class PolarCompiler
{
	static const int U = 2;
	static uint8_t leaf(int frozen)
	{
		return frozen;
	}
	static uint8_t node(int func, int level)
	{
		return (func << 5) | level;
	}
	static uint8_t left(int level)
	{
		return node(1, level);
	}
	static uint8_t right(int level)
	{
		return node(2, level);
	}
	static uint8_t combine(int level)
	{
		return node(3, level);
	}
	static uint8_t rate0(int level)
	{
		return node(4, level);
	}
	static uint8_t rate1(int level)
	{
		return node(5, level);
	}
	static uint8_t rep(int level)
	{
		return node(6, level);
	}
	template<typename TYPE>
	static int popcnt(TYPE x)
	{
		int cnt = 0;
		while (x) {
			++cnt;
			x &= x-1;
		}
		return cnt;
	}
	static int frozen_count(const uint8_t *frozen, int level)
	{
		int count = 0;
		for (int i = 0; i < 1<<(level-U); ++i)
			count += popcnt(frozen[i]);
		return count;
	}
	static void compile(uint8_t **program, const uint8_t *frozen, int level)
	{
		if (level > U) {
			int count = frozen_count(frozen, level);
			if (count == 1<<level) {
				*(*program)++ = rate0(level);
			} else if (count == 0) {
				*(*program)++ = rate1(level);
			} else if (count == (1<<level)-1 && frozen[(1<<(level-U))-1] == (1<<((1<<U)-1))-1) {
				*(*program)++ = rep(level);
			} else {
				*(*program)++ = left(level);
				compile(program, frozen, level-1);
				*(*program)++ = right(level);
				compile(program, frozen+(1<<(level-1-U)), level-1);
				*(*program)++ = combine(level);
			}
		} else {
			*(*program)++ = leaf(*frozen);
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

