#include <cstdint>
#include <algorithm>

/*
A 64-bit Xorshift PRNG combined a Weyl Generator.
From ideas of Marsgalia, Brent, and others.
Passes all the BigCrush tests.
*/
class xorshift64 {
public:
	xorshift64(uint64_t seed1 = 0, uint64_t seed2 = 0) {	
		seed(seed1,seed2);
	}

	explicit xorshift64(std::pair<uint64_t,uint64_t> p) {
		seed(p.first,p.second);
	}

	void seed(uint64_t seed1 = 0, uint64_t seed2 = 0) {
		u = (seed1 == 0) ? UINT64_C(15191868757011070976) : seed1;
		w = (seed2 == 0) ? UINT64_C(0x61C8864680B583EB) : seed2;
		get_raw();
	}

	void seed(std::pair<uint64_t,uint64_t> p) {
		seed(p.first,p.second);
	}
	
	std::pair<uint64_t,uint64_t> get_state() const {
		return std::make_pair(u,w);
	}
	
	// Xorshift + Weyl Generator
	uint64_t get_raw() {
		u ^= (u << 5); u ^= (u >> 15); u ^= (u << 27);
		w += UINT64_C(0x61C8864680B583EB);
		return u+w;
	}
	
	uint64_t get_uint64() {
		return get_raw();
	}
	
	uint32_t get_uint32() {
		return static_cast<uint32_t>(get_raw() >> 32);
	}
	
	// Uniform [0,1)
	double get_double53() {
		union { uint64_t u; double d; } a;
		a.u = get_raw();
		a.u = (a.u >> 12) | UINT64_C(0x3FF0000000000000);
		double q = (a.u&2048) ? (1.0-(DBL_EPSILON/2.0)) : 1.0;
		return a.d-q;
	}

	// Uniform (0,1)
	double get_double52() {
		union { uint64_t u; double d; } a;
		a.u = get_raw();
		a.u = (a.u >> 12) | UINT64_C(0x3FF0000000000000);
		double q = (1.0-(DBL_EPSILON/2.0));
		return a.d-q;
	}
	
	uint64_t operator()() {
		return get_raw();
	}
	
private:
	uint64_t u,w;
};
