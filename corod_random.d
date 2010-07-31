module corod_random;

import std.random;

import corod;

/// Interface of generator of random numbers
interface IRandoms(T = uint) : Generator!(T) {
	///
	T max();
	///
	uint period_log2();
}

/// Interface of parallal random strategies
/// TODO: block skip, leapfrog and random seeding
interface IParallelRandom(T = uint) {
	T max_N();
	IRandoms!(T) get(uint generator_no);
}


/** It is simple Linear Congruential Generator with formula:
 *
 * x = (69069*x+1234567) % 2^^32
 *
 * It is quite poor, but can be used as parts of other generators,
 * as initialisation or with combination with other generators.
 * Short period.
 *
 * Last half of the bits of this generator are too regular.
 */
final class RandomsCONG : FiberGenerator!(uint), IRandoms!(uint) {
public:
	/** Initialize random number generator using unpredictable seed */
	this() {
		this(unpredictableSeed());
	}

	/// Initialize random number generator using seed s
	this(uint s_) {
		assert(s_ < int.max);
		s = s_;
	}

	/// Returns maximal generated value. It is 2**31-1
	uint max() const {
		return uint.max;
	}
	/// period
	uint period_log2() const {
		return 32;
	}

private:
	uint s = 1;

protected:
	/// This is infinite generator, yielding random number
	override void iter() {
		uint x = s;
		while (true) {
			x = 69069*x + 1234567;
			yield(x);
		}
	}
}



/** Marsgalia's CMWC4096 (complimentary-multiply-with-carry) generator.
 *
 * It is outputs base b=2^32-1 expansion of k/p for some k, in reverse order.
 *
 * p is prime equal 18782*b^4096 + 1
 * persiod is p-1 which is more than 10^39460 (2^131086)
 * It requires 4096 seeds x and seed c
 *
 * Attention: Is this generator belived to be broken? There are few publications
 * which states that this (or something named Marsaglia Multicarry
 * fail on many tests).
 */
final class RandomsCMWC4096 : FiberGenerator!(uint), IRandoms!(uint) {
public:
	/** Initialize random number generator using unpredictable seed */
	this() {
		this(unpredictableSeed());
	}

	/// Initialize random number generator using seed s
	this(uint s) {
		foreach (i, ref xsi; Q) {
			xsi = (s = (69069*s+1234567));
		}
		c = (s = (69069*s+1234567));
	}

	/// Returns maximal generated value. It is probably 2**31-1
	uint max() const {
		return uint.max;
	}

	/// period
	uint period_log2() const {
		return 131086;
	}

private:
	uint Q[4096];
	uint c = 362436;


protected:
	/// This is infinite generator, yielding random number
	override void iter() {
		const ulong a = 18_782uL;
version(cmwc1) {
		const ulong b = 4_294_967_295uL;
		const uint r = b-1;
} else {
		const uint r = 0xfffffffeu;
}
		uint i = 4095;
		while (true) {
			i = (i+1) & 4095;
			ulong t = a*Q[i] + c;
			c = (t>>32);
version(cmwc1) {
			t = (t & b) + c;
			if (t > r) {
				c++;
				t = t - b;
			}
			Q[i] = r - t;
} else {
			uint x = cast(uint)(t + c);
			if (x < c) {
				x++;
				c++;
			}
			Q[i] = r - x;
}
			yield(Q[i]);
		}
	}
}

/* SWB
 *
 * x(n) = x(n-222) - x(n-237) - borrow mod 2^32.
 * borrow is 1 is computing x(n-1) caused overflow in 32-bit integer arithmetic.
 * period ~ 2^7578
 * SWB failes Birthday spacing test.
 *
 * ulong[256] t;
 * ubyte c=0;
 * ulong bro,x=0,y=0;
 * (c++,bro=(x<y),t[c]=(x=t[UC(c+34)])-(y=t[UC(c+19)]+bro))
 *
 * Use only in combination with other generators, like SWB+KISS.
 * KISS+SWB has period >2^770 and is highly recomended.
 */

/* SHR3
 *
 * (jsr^=(jsr<<17), jsr^=(jsr>>13), jsr^=(jsr<<5))
 *
 * period 2^32-1.
 * passes most tests, but fails this releated to binary rank tests.
 */

/// combination of three generators
final class RandomsKISS : FiberGenerator!(uint), IRandoms!(uint) {
private:
	uint[4096] Q;
	uint c = 362436;

public:
	/** Initialize random number generator using unpredictable seed */
	this() {
		this(unpredictableSeed());
	}

	/// Initialize random number generator using seed s
	this(uint s_) {
		s = s_;
	}

	/// Returns maximal generated value. It is probably 2**31-1
	uint max() const {
		return uint.max;
	}
	/// period 2^^123
	uint period_log2() const {
		return 123;
	}

private:
	uint s;

protected:
	/// This is infinite generator, yielding random number
	override void iter() {
		uint
			x = s+123456789u,
			y = (69069*x+1234567)+362436000u,
			z = (69069*y+1234567)+521288629u;
		uint c = (69069*z+1234567)+7654321u;
		immutable ulong a = 698769069uL;
		while (true) {
			x = 69069u*x + 12345u; // CONG  originaly with 1234567
			y ^= (y<<13); // SHR3, orginally with slightly other order
			y ^= (y>>17);
			y ^= (y<<5);
			ulong t = a*z + c; // MWC
			c = (t>>32);
			z = cast(uint)t;
			//yield(x + y + z); // other variant
			yield(x^z + y); // like in 1999 post.
		}
	}
}

/** Tausworthe generator with period 2^88.
 *
 * P. L'Ecuyer, “Maximally Equidistributed Combined Tausworthe Generators”, Mathematics of Computation, 65, 213 (1996), 203–213.
 * P. L'Ecuyer, “Tables of Maximally Equidistributed Combined LFSR Generators”, Mathematics of Computation, 68, 225 (1999), 261–269
 */
final class RandomsTaus88 : FiberGenerator!(uint), IRandoms!(uint) {
public:
	/** Initialize random number generator using unpredictable seed */
	this() {
		this(unpredictableSeed());
	}

	/// Initialize random number generator using seed s
	this(uint s_) {
		if (s_ < 128) {
			s_ <<= 16;
		}
		s1_ = s_;
		s2_ = (69069*s1_+1234567);
		s3_ = (69069*s2_+1234567);
	}

	/// Returns maximal generated value. It is 2**32-1
	uint max() const {
		return uint.max;
	}
	// period
	uint period_log2() const {
		return 88;
	}

private:
	uint s1_, s2_, s3_;

protected:
	/// This is infinite generator, yielding random number
	override void iter() {
		// this inialisation procedure assumes k+r>=L=32   (r == k-q)

		static assert(2*(12+19)-13 >= 32);
		static assert(4294967294u == (0xffffffffu << (32-12-19)));
		uint s1 = (s1_ & 4294967294u); // C= (0xffffffffu << (32-k)), mask of k 1, and L-k 0 
		s1 = s1 ^ (((s1 << 13) ^ s1) >> (12+19)/*k=31*/); // A = A ^ (((A << q) ^ A) >> k)
		assert((s1 & 4294967294u) != 0);

		static assert(2*(4+25)-2 >= 32);
		static assert(4294967288u == (0xffffffffu << (32-4-25)));
		uint s2 = (s2_ & 4294967288u);
		s2 = s2 ^ (((s2 << 2) ^ s2) >> (4+25)/*k=29*/);
		assert((s2 & 4294967288u) != 0);

		static assert(2*(17+11)-3 >= 32);
		static assert(4294967280u == (0xffffffffu << (32-17-11)));
		uint s3 = (s3_ & 4294967280u);
		s3 = s3 ^ (((s3 << 3) ^ s2) >> (17+11)/*k=28*/);
		assert((s2 & 4294967280u) != 0);

		// above initalisation procedure isn't really nacasary if L-k <= r-s=k-q-s, just skip first value of generator
		static assert(32-12-19 <= 12+19-13-12);
		static assert(32-4-25 <= 4+25-2-4);
		static assert(32-17-11 <= 17+11-3-17);
		while (true) {
			// QuickTaus // ((A & C) << s) ^ (((A << q) ^ A) >> (k-s))
			s1 = ((s1 & 4294967294u) << 12) ^ (((s1 << 13) ^ s1) >> 19);
			s2 = ((s2 & 4294967288u) << 4) ^ (((s2 << 2) ^ s2) >> 25);
			s3 = ((s3 & 4294967280u) << 17) ^ (((s3 << 3) ^ s3) >> 11);
			yield(s1 ^ s2 ^ s3);
		}
	}
}

/*
 * Four component LFSRs:
 *
 * L=32
 * 2^113
 * J=4
 * k1=31 k2=29 k3=28 k4=25
 * q1=6 q2=2 q3=13 q4=3
 * s1=18 s2=2 s3=7 s4=13
 *
 * Five components: (note this specific is for 64bits machines mainly).
 *
 * L = 64
 * 2^258
 * J=5
 * k1=63 k2=55 k3=52 k4=47 k5=41
 * q1=1 q2=24 q3=3 q4=5 q5=3
 * s1=10 s2=5 s3=29 s4=23 s5=8
 *
 */

/* Also known as GFSR4
 *
 * Robert M. Ziff, “Four-tap shift-register-sequence random-number generators”, Computers in Physics, 12(4), Jul/Aug 1998, pp 385–392.
 */
final class RandomsZiff : FiberGenerator!(uint), IRandoms!(uint) {
public:
	/** Initialize random number generator using unpredictable seed */
	this() {
		this(unpredictableSeed());
	}

	/// Initialize random number generator using seed s
	this(uint s) {
		foreach (ref x; r) {
			x = (s = 69069*s+1234567);
		}
	}

	/// Returns maximal generated value. It is 2**32-1
	uint max() const {
		return uint.max;
	}
	/// period, 2^^9689
	uint period_log2() const {
		return 9689;
	}

private:
	// 64kB of memory
	uint[16384] r; // smallest power of two enaugh
	// <del>smallest multiplicity 64 enaugh for us</dev>

protected:
	/// This is infinite generator, yielding random number
	override void iter() {
		// This generator is equivalent to 7 times decimated generator LFib(471,9689,xor),
		// decimation is to remove inherent 3-point correlations.
		// Other quit good results have R(50,103,200,250) decimation from femouse F(103,25)
		// and R(216,337,579,1279)
		// tests performed on few physical MC computations, shows
		// it is quite good, but still there are some weak corelations.
		uint i1 = r.length-471;
		uint i2 = r.length-1586;
		uint i3 = r.length-6988;
		uint i4 = r.length-9689;
		uint n = 0;
		immutable mask = (r.length - 1);
		while (true) {
			immutable new_r = r[i1++] ^ r[i2++] ^ r[i3++] ^ r[i4++];
			r[n++] = new_r;
			n &= mask;
			i1 &= mask;
			i2 &= mask;
			i3 &= mask;
			i4 &= mask;
			yield(new_r);
		}
	}
}


/*  Marsaglia's Lagged Fibonacci
 *
 * http://www.ciphersbyritter.com/NEWS4/RANDC.HTM
 *
 * KISS+LFIB4 have period 2^410 and is quite interesting.
 */
final class RandomsLFib4 : FiberGenerator!(uint), IRandoms!(uint) {
public:
	/** Initialize random number generator using unpredictable seed */
	this() {
		this(unpredictableSeed());
	}

	/// Initialize random number generator using seed s
	this(uint s) {
		foreach (ref x; t) {
			x = (s = 69069*s+1234567);
		}
	}

	/// Returns maximal generated value. It is 2^^32-1
	uint max() const {
		return uint.max;
	}
	// Period is 2^^31*(2^^256-1).
	uint period_log2() const {
		return 287; // about
	}

private:
	// 1kB of memory
	uint t[256];

protected:
	/// This is infinite generator, yielding random number
	override void iter() {
		ubyte c = 0;
		while (true) {
			c++;
			t[c] = t[c] + t[cast(ubyte)(c+58)] + t[cast(ubyte)(c+119)] + t[cast(ubyte)(c+178)];
			yield(t[c]);
		}
	}
}

/*
 Other good generators:
 LFib(2^64, a, b, *) - przechodza wszystkie testy. nie sa tu zaimplementowane,
 * bo w zasadzie wymagaja 64bitowych instrukcji.
  pary a,b: 17,5
  55,24
  607,273
  1279,861
 */

/* Knuth's GFSR using lagged Fibonaci. Version from 2002 edition of TAOCP.
 *
 * It is also used in R.
 *
 * Attention seeding and actuall implementation is quite different.
 * http://sunburn.stanford.edu/~knuth/programs/rng.c
 */
final class RandomsKnuthTAOCP2002: FiberGenerator!(uint), IRandoms!(uint) {
public:
	/** Initialize random number generator using unpredictable seed */
	this() {
		this(unpredictableSeed());
	}

	/// Initialize random number generator using seed s
	this(uint s) {
		foreach (ref xv; x) {
			xv = (s = 69069*s+1234567);
		}
	}

	/// Returns maximal generated value. It is 2^^30-1
	uint max() const {
		return ((1uL<<30) - 1);
	}

	/// Returns binary period log, which is about 2^^129
	uint period_log2() const {
		return 129;
	}

private:
	// 512B of memory
	uint[128] x;

protected:
	/// This is infinite generator, yielding random number
	override void iter() {
		uint i1 = x.length-100;
		uint i2 = x.length-37;
		uint n = 0;
		immutable mask = (x.length - 1);
		while (true) {
			immutable new_x = (x[i1++] - x[i2++]) % ((1<<30)-1);
			x[n++] = new_x;
			n &= mask;
			i1 &= mask;
			i2 &= mask;
			yield(new_x);
		}
	}
}


/** Delinearization, which uses generator x.
 *
 * Having generator giving x_i, generate:
 *  r_i = g^x_i mod m.
 *
 * Exponentation is done by table lookups.
 * Exponentation by squering is too slow.
 *
 *  r_i = g^(2^B H + L) mod m = (g^(2^B))^H * g^L mod m
 *
 * with (g^(2^B))^H and g^L in two 2^16,2^15 tables. (about 385 kB)
 *
 * Modulo operation is performed in possibly fastests way.
 *
 */
final class RandomsDelinearize : Generator!(uint), IRandoms!(uint) {
private:
	uint gH[65536];
	uint gL[32768];
	Generator!(uint) G;
	uint x;
public:
	/// IRandoms interface
	uint max() const {
		return uint.max;
	}
	/// ditto
	uint period_log2() const {
		return 32;
	}
	/// Generator interface
	bool next() {
		x = G.getNext();
		return true;
	}
	/// ditto
	uint get() const {
		return x; /// TODO: implement
	}
	/// ditto
	uint getNext() {
		next();
		return get();
	}

	mixin OpApplyMixin!(uint) opApply;
}


/** Generates random numbers of type uint, in the interval [0, max) 
 *
 *
 * It will use Mersenne Twister 19937 from Phobos.
 *
 */
final class Randoms(T = uint) : FiberGenerator!(T), IRandoms!(T) {
private:
	/// Internal random number generator engine
	//Random gen;
	Mt19937 gen; // good random generator
	//MinstdRand0 gen; /// "minimal standard" random number generator from Park-Miller, poor
	//MinstdRand gen; /// also poor

public:
	/** Initialize random number generator using unpredictable seed */
	this() {
		gen.seed(unpredictableSeed());
	}

	/// Initialize random number generator using seed s
	this(T s) {
		gen.seed(s);
	}

	/// Returns maximal generated value. It is 2^^32-1
	T max() const {
		return gen.max;
	}
	/// period, 2^^19937-1
	uint period_log2() {
		return 19937;
	}

protected:
	/// This is infinite generator, yielding random number
	override void iter() {
		while (true) {
			gen.popFront();
			auto r = gen.front;
			yield(r);
		}
	}
}


/** Combines two calls to Random(uint) to produce Random(ulong).
 *
 * Assumes that Random!(uint) produces all possible uint's.
 *
 * ulong generators is very usefull when using float or double uniforms.
 */
final class RandomsUL : FiberGenerator!(ulong), IRandoms!(ulong) {
private:
	IRandomss!(uint) g;
public:
	///
	this(IRandoms!(uint) g_) {
		g = g_;
	}
protected:
	override void iter() {
		while(true) {
			auto hi = g.getNext();
			auto lo = g.getNext();
			yield((hi<<32) | lo);
		}
	}
public:
	///
	ulong max() const {
		return ulong.max;
	}
	///
	uint period_log2() {
		return g.log2_period();
	}
}


Generator!(uint) combine_xor(T=uint,uint sa=0,uint sb=0)(Generator!(uint) a, Generator!(uint) b) {
	return new ConcurantGenerator!(uint, uint, uint)(a, b,
		delegate uint(uint x, uint y) { return (x << sa) ^ (y << sb); }
	);
}



/** Generates uniformly distributed random float number
 * in inteval [a; b) using uniform random number generator
 */
final class RandomsUniform(T = float) : FiberGenerator!(T) {
private:
	IRandoms!(uint) iterator;
	const /*T*/ double A, B;

public:
	/** Generates floats [a; b) */
	this(IRandoms!(uint) iterator_, T a_ = 0.0, T b_ = 1.0)
	in {
		assert(a_ < b_);
		assert(iterator_ !is null);
	}
	body {
		iterator = iterator_;
		A = ((b_ - a_)/(iterator.max + 1.0));
		B = a_;
	}

/+
	/** Generates floats [a; b) */
	this(Generator!(T) iterator_, T a_ = 0.0, T b_ = 1.0)
	in {
		assert(a_ < b_);
		assert(iterator_ !is null);
	}
	body {
		iterator = iterator_;
		A = (b_ - a_);
		B = a_;
	}
+/

	/// Same as above, but determinate automatically some new random number generator with unpredictable seed
	this(T a_ = 0.0, T b_ = 1.0) {
		this(new Randoms!(uint)(), a_, b_);
	}

	/// Same as above, but determine automatically some new random number generator with seed 'seed'
	this(uint seed, T a_ = 0.0, T b_ = 1.0) {
		this(new Randoms!(uint)(seed), a_, b_);
	}

protected:
	/** implementation
	 *
	 * Note: this is generic Map iterator
	 */
	override void iter() {
		foreach (uint x; new Foreach!(uint, iterator)()) {
			//debug assert(0 <= x);
			//debug assert(x <= gen.max);
			//debug assert(a_ <= A*x+B && A*x+B < b_);
			immutable T r = A*x+B;
			if (r == 0.0 || r == 1.0) continue;
			yield(r);
		}
	}
}


import std.math : sqrt, log, cos, sin, expi, PI, tan, exp, abs;

/** Cauchy distribution.
 *
 * ---
 * p(x) = \frac{1}{\pi} {1 \over 1 + x^2}
 * ---
 *
 * ---
 * u = F(X) = \int_{-\infty}^x p(x) = 1/2 + \arctan(x)/\pi
 *
 * x = \tan(\pi(u-\frac{1}{2}))
 * ---
 */
final class RandomsStandardCauchy(T = float) : FiberGenerator!(T) {
private:
	Generator!(T) G;

public:
	/// Generator G need to be uniform float generator on interval (0; 1)
	this(Generator!(T) G_) {
		G = G_;
		super();
	}

protected:
	/// implementation
	override void iter() {
		while (true) {
			yield(tan(PI*(G.getNext() - 0.5)));
		}
	}
}

/** Box-Muller's algorithm for generating Gauss distribution.
 *
 * ---
 * X_1, X_2 \in (0, 1)
 * 
 * Y_1 = \mi + \rho\sqrt{-2 \log X_1} \cos(2\pi X_1)
 * Y_2 = \mi + \rho\sqrt{-2 \log X_1} \sin(2\pi X_1)
 * ---
 *
 * Note: remember about simulatanius calculation of sin(x)/cos(x) using e^{ix}
 *
 * Eventually let V_1, V_2 \in circle(0;1)
 * then
 * ---
 * Y_1 = \sqrt{-2 \log X_1} (V_1 / R)
 * Y_2 = \sqrt{-2 \log X_1} (V_2 / R)
 * ---
 */
final class RandomsStandardNormalBoxMuller(T = float) : FiberGenerator!(T) {
private:
	Generator!(T) G;

public:
	/// Generator G need to be uniform float generator on interval (0; 1)
	this(Generator!(T) G_) {
		G = G_;
		super();
	}

	/// Automatically determine new uniform random number generator and initialize it to some unpredictable seed
	this() {
		this(new RandomsUniform!(T)());
	}

	/// Automatically determine new uniform random number generator and initlize it using seed
	this(uint seed) {
		this(new RandomsUniform!(T)(seed));
	}

protected:
	/// implementation
	override void iter() {
		while (true) {
			auto X_1 = G.getNext();
			auto X_2 = G.getNext();
			debug assert(0.0 < X_1 && X_1 < 1.0, "Generator G should generate (0;1)");
			debug assert(0.0 < X_2 && X_2 < 1.0, "Generator G should generate (0;1)");
			auto cossin = expi((2.0*PI)*X_2);
			auto s = sqrt(-2.0*log(X_1));
			yield(s*cossin.re);
			yield(s*cossin.im);
		}
	}
}

/** Same as above, but for not standard values of mi and rho */
final class RandomsGaussBoxMuller(T = float) : FiberGenerator!(T) {
private:
	Generator!(T) G;
	const T mi, rho;

public:
	/// generic constructor
	this(Generator!(T) G_, T mi_ = 0.0, T rho_ = 1.0) {
		G = G_;
		mi = mi_;
		rho = rho_;
		super();
	}

	/// Automatically determine new uniform random number generator and initialize it to some unpredictable seed
	this(T mi_ = 0.0, T rho_ = 1.0) {
		this(new RandomsUniform!(T)(), mi_, rho_);
	}

	/// Automatically determine new uniform random number generator and initlize it using seed
	this(uint seed, T mi_ = 0.0, T rho_ = 1.0) {
		this(new RandomsUniform!(T)(seed), mi_, rho_);
	}

protected:
	/// implementation
	override void iter() {
		while (true) {
			auto X_1 = G.getNext();
			auto X_2 = G.getNext();
			debug assert(0.0 < X_1 && X_1 < 1.0, "Generator G should generate (0;1)");
			debug assert(0.0 < X_2 && X_2 < 1.0, "Generator G should generate (0;1)");
			auto cossin = expi((2.0*PI)*X_2);
			auto rho_s = rho*sqrt(-2.0*log(X_1));
			yield(mi + rho_s*cossin.re);
			yield(mi + rho_s*cossin.im);
		}
	}
}

/** Generates Gaussian random nubmer, just like RandomGauss but with other algorithm.
 *
 * It doesn't compute cos/sin (or expi), but involves possible more operations
 * in underling number generator.
 *
 * It is 
 *
 * TODO: Benchmark needed.
 *
 * Note: Is also known as "Marsaglia polar method"
 *
 * See: Bell 1968, Knop 1969
 */
final class RandomsGaussBoxMullerFast(T = float) : FiberGenerator!(T) {
private:
	Generator!(T) G;
	const T mi, rho;

public:
	/// generic constructor
	this(Generator!(T) G_, T mi_ = 0.0, T rho_ = 1.0) {
		G = G_;
		mi = mi_;
		rho = rho_;
		super();
	}

	/// Automatically determine new uniform random number generator and initialize it to some unpredictable seed
	this(T mi_ = 0.0, T rho_ = 1.0) {
		this(new RandomsUniform!(T)(), mi_, rho_);
	}

	/// Automatically determine new uniform random number generator and initlize it using seed
	this(uint seed, T mi_ = 0.0, T rho_ = 1.0) {
		this(new RandomsUniform!(T)(seed), mi_, rho_);
	}

protected:
	/// implementation
	override void iter() {
		while (true) {
			auto X_1 = G.getNext();
			auto X_2 = G.getNext();
			debug assert(0.0 < X_1 && X_1 < 1.0, "Generator G should generate (0;1)");
			debug assert(0.0 < X_2 && X_2 < 1.0, "Generator G should generate (0;1)");
			auto V_1 = 2.0*X_1 - 1.0;
			auto V_2 = 2.0*X_2 - 1.0;
			auto R2 = V_1*V_1 + V_2*V_2;
			if (R2 < 1.0) {
				auto rho_s = rho*sqrt(-2.0*log(R2)/R2);
				yield(mi + rho_s*V_1);
				yield(mi + rho_s*V_2);
			}
		}
	}
}

alias RandomsStandardNormalBoxMuller RandomsStandardGaussT;
alias RandomsGaussBoxMuller RandomsGaussT;

alias RandomsStandardGaussT!(float) RandomsStandardGauss;
alias RandomsGaussT!(float) RandomsGauss;

/** TODO: The Gaussian Tail Distribution
 *
 * Standard Box-Muller algorithm is bad for tails, (log(0+epsilon) problem),
 *
 * Other method for Gauss: This function computes a Gaussian random variate using the alternative Marsaglia-Tsang ziggurat and Kinderman-Monahan-Leva ratio methods. The Ziggurat algorithm is the fastest available algorithm in most cases.
 * Other method for tail:  The method is based on Marsaglia's famous rectangle-wedge-tail algorithm (Ann. Math. Stat. 32, 894.899 (1961)), with this aspect explained in Knuth, v2, 3rd ed, p139,586 (exercise 11).
 *
 * Also method for integrated distrbuant, and 1-distribuant.
 */

/** Ziggurat methods for unimodal, monotone decressing pdf's
 *
 * Note: pdf f doesn't need to be normalized
 * Note: for T=float, generated uint is separated into two parts:
 *      24 bit for mantisa in float [0,1)
 *      and 7 bit part for choising one of 128 layers (127 + tail)
 *      and 1 bit is for sign
 *
 * See Marsaglia and Tsang, 1984a, 2000.
 *
 * "The Ziggurat Method for Generating Random Variables", Marsaglia, Tsang
 * Faster and simpler than previous version of Ziggurat method
 *
 * "Design Flaws in the Impementation of the Ziggurat and Monty Python methods ..." Boaz Nadler
 *
 * TODO: see Marsaglia, G. and Tsang, W.W., 1998. The monty python method for generating random variables. ACM Trans. Math. Software 24 3, pp. 341–350.
 */
final class RandomsZigguratGauss(T = float) : FiberGenerator!(T) {
private:
	Generator!(uint) G;

	// we are using 128, becaus we need: 24 bit for (0,1) float, 7 bits for tables wn/fn, and 1 bit for sign
	static uint[128] kn;
	static T[128] wn, fn;

public:
	/// generic constructor
	this(Generator!(uint) G_) {
		G = G_;
		super();
	}

private:
	static double f(double x) {
		return exp(-0.5*x*x); // hmm, unnormalized
	}
	static float f(float x) {
		return exp(-0.5f*x*x);
	}
	static double f_inv(double y) {
		return sqrt(-2.0*log(y));
	}


public:
	static this() {
		//writefln("konstruktor");

		// todo: initialize (static for given f) tables
		// normal distribution:
		const double m1 = 2_147_483_648.0;
		double dn = 3.442_619_855_899;
		const double vn = 9.912_563_035_262_17e-3;
		double tn = dn;
		double q = vn/f(dn);
		kn[0] = cast(uint)((dn/q)*m1);
		kn[1] = 0;
		wn[0] = q/m1;
		wn[127] = dn/m1;
		fn[0] = 1.0;
		double fv = f(dn);
		fn[127] = fv;
		for (int i = 126; i >= 1; i--) {
			dn = f_inv(vn/dn + fv);
			kn[i+1] = cast(uint)((dn/tn)*m1);
			tn = dn;
			fv = f(dn);
			fn[i] = fv;
			wn[i] = dn/m1;
		}
/+
		foreach (i; 0 .. 128) {
			writefln("%d kn[%d]= %d fn[%d]= %.10g wn[%d]= %.10g", i, i, kn[i], i, fn[i], i, wn[i]);
		}
+/
	}

private:
	static long long_abs(long x) {
		return (x >= 0 ? x : -x);
	}

	T UNI() {
		//return 0.5 + cast(int)G.getNext() * 0.232_830_6e-9f;
		return G.getNext() / (1.0 + uint.max);
	}


protected:

/+
	override void iter() {
		while (true) {
			uint u_hz = G.getNext();
			long i_hz = cast(long)(u_hz & 0x7fff_fff0u);
			long hz = ((u_hz & 0x8000_0000u) ? i_hz : -i_hz);
			ulong iz = u_hz & 0x7fuL;
			uint iz2 = cast(uint)iz;
			assert(long_abs(hz) == i_hz);
			yield(i_hz < kn[iz2] ? hz*wn[iz2] : normal_fix(hz, iz2));
		}
	}
+/

	/// implementation
	override void iter() {
		while (true) {
			auto hz = G.getNext()>>1;
			auto g2 = G.getNext();
			auto i = cast(uint)(g2 & 0x7fuL);
			auto sign = cast(uint)(g2 & 0x80uL);

			T x = hz*wn[i];
			if (hz < kn[i]) { // check if we are on the left of boxes
version (testing_acceptance) {
				fast_accept++;
				hist_fast[i]++;
}
				yield(sign ? x : -x);
			} else {
				immutable T r = 3.442_619_855_899; // my modification, taken from dn in static ctor
				immutable T minus_inv_r = -1.0/r;
				if (i == 0) { // base strip (tail distribution)
version (testing_acceptance) {
					tail_test++;
}
					T y;
					do {
version (testing_acceptance) {
						tail_iter++;
}
						// x = r - log(UNI()); // exponential
						x = minus_inv_r*log(UNI()); // normal
						y = -log(UNI());
					} while (y+y < x*x);
					//yield( (hz > 0) ? r+x : -r-x );
					yield(sign ? r+x : -(r+x));
				} else { // we are in the boxes, perform rejection testing, or repeat.
					if ((fn[i-1]-fn[i])*UNI() < f(x) - fn[i]) {
version (testing_acceptance) {
						slow_accept++;
						hist_slow_accept[i]++;
}
						yield(sign ? x : -x);
					} else {
version (testing_acceptance) {
						slow_reject++;
						hist_slow_reject[i]++;
}
					}
				}
			}
		}
	}

version (testing_acceptance) {
public:
	ulong fast_accept, slow_accept, slow_reject, tail_test, tail_iter;
	ulong hist_fast[128];
	ulong hist_slow_accept[128];
	ulong hist_slow_reject[128];

	~this() {
		writefln("fast_accept %d",fast_accept);
		writefln("slow_accept %d",slow_accept);
		writefln("slow_reject %d",slow_reject);
		writefln("tail_test %d",tail_test);
		writefln("tail_iter %d",tail_iter);
		for (int i =0; i < 128; i++) {
			auto total=hist_fast[i]+hist_slow_accept[i]+hist_slow_reject[i];
			auto total_slow=hist_slow_accept[i]+hist_slow_reject[i];
			writefln("i=%d fast_hit=%d slow_hit_accept=%d slow_hit_reject=%d tot=%d slow_tot=%d fast=%.4f%% slow=%.4f%%",i,hist_fast[i],hist_slow_accept[i],hist_slow_reject[i], total, total_slow, 1.0*hist_fast[i]/total, 1.0*hist_slow_accept[i]/total_slow);
		}
	}
}

private:
/+
	T normal_fix(long hz, uint iz2) {
		//const T r = 3.442_620f; // start of the right tail
		const T r = 3.442_619_855_899; // my modification, taken from dn in static ctor
		while (true) {
			T x = hz*wn[iz2];   //iz2==0 is for base strip */
			if (iz2 == 0) {
				T y;
				do {
					x = -log(UNI())*0.290_476_4f; /* 0.2904764 is 1/r*/
					y = -log(UNI());
				} while (y+y < x*x);
				return ( (hz > 0) ? r+x : -r-x );
			} else if (fn[iz2] + UNI()*(fn[iz2-1]-fn[iz2]) < exp(-0.5*x*x)) {
				return x;
			} else {
				uint u_hz = G.getNext();
				long i_hz = cast(long)(u_hz & 0x7fff_ffffu);
				hz = ((u_hz & 0x8000_0000u) ? i_hz : -i_hz);
				ulong iz = hz & 0x7fuL;
				iz2 = cast(uint)iz;

				assert(i_hz == long_abs(hz));

				if (i_hz < kn[iz2]) {
					return hz*wn[iz2];
				}
			}
		}
	}
+/
}

import ctfe_logexp : ctfe_log;

/** Monty Python method for generating Gaussian varietes by Marsaglia
 *
 * Idea is that pdf is partitioned into 4 rectangle regions (one of them infinite)
 * One of them is fully filled with pdf.
 * Two other can be quite easly be packed into one rectangle.
 * Forth is a tail.
 *
 * It can be seen as rejection method, but more appropriate is
 * calling it a transormation method.
 *
 * See: Marsaglia and Tsang, 1998
 */
final class RandomsMontyPythonGauss(T = float) : FiberGenerator!(T) {
private:
	Generator!(float) G;

public:
	/// generic constructor
	this(Generator!(float) G_) {
		G = G_;
		super();
	}

protected:
	/// implementation
	override void iter() {
		static immutable b = sqrt(2.0*PI); // this is allows sampling from tail only in 1.2% cases
		static immutable log4 = ctfe_log(4.0);
		static immutable a = sqrt(log4); // normal_icdf(1.0/(2.0*b));
		static immutable s = a/(b-a); // stretch factor to make area C as fit exactly above A

		// version without stretching (s == 1), is possible up to b=2.29. tail hit is then 2.2%

		// we are uing 2.0*normal_pdf, becuase we are only testing for right
		// side of distribution

		// screatched and rotated cap (area C')
		static if (s == 1.0) {
			static T g2(T x) {
				// todo: common expressions, (b-x) and 1/b
				return 2.0*(normal_pdf(b-x) - 1.0/b);
			}
		} else {
			static T g(T x) {
				// todo: common expressions, s*(b-x) and 1/b
				return 1.0/b - s*(2.0*normal_pdf(s*(b-x)) - 1.0/b);
			}
		}

		// we will be sampling rectangle b x 0.5/b
		while (true) {
			auto sign = 2*cast(int)(2.0*G.getNext()) - 1;
			auto x = b*G.getNext(); // horizontal component
			if (x < a) { // area A (about 47.9% hit probabilitiy)
				yield(sign*x);
			} else {
				auto y = G.getNext()/(b);
/+
				auto v = 2.8658 - x*(2.0213 - 0.3605*x);
				if (y < v) { // pretest for area B, should speed us up slightly because of lack of log,exp,sqrt
					yield(sign*x); // this quadratic pretest is sufficient 98% of the times
				} else
+/
				if (y < 2.0*normal_pdf(x)) { // exact test for area B, executed about 2%
					// other option is not test y < pdf(x),
					// but log(y) < log(2.0) - 0.5*x*x;
					yield(sign*x);
				} else {
					// TODO: pretest for C'
					// map area C' (packed in B) to C
					if (y > g(x)) { // area C
						yield(sign*s*(b-x)); // unmap
					} else { // area D', generate from tail (tranformation D'->D is too complex). about 1.2% hit probability
						yield(sign*normal_gen_tail!(T)(b, G));
					}
				}
			}
		}
	}
}



// FN3(q) [Wallace 1996]
// VNAL, alias method [ Ahrens and Dieter 1988]

/*
final class RandomsZigguratGauss(T = float) : FiberGenerator!(T) {
private:
	Generator!(uint) G;
	T delegate(T) f;
	T delegate(T) inv_f;

	// we are using 128, becaus we need: 24 bit for (0,1) float, 7 bits for tables wn/fn, and 1 bit for sign
	T[128] kn, wn, fn;

public:
	/// generic constructor
	this(Generator!(uint) G_, T delegate(T) f_, T delegate(T) inv_f_) {
	}
}
*/

/** Gauss distribution using Kinderman-Ramage algorithm
 * (picewise triangular) with fixes.
 *
 * Computer Generation of Normal Random Variables
 * A. J. Kinderman and J. G. Ramage
 * Journal of the American Statistical Association, Vol. 71, No. 356 (Dec., 1976), pp. 893-896
 *
 * Gunter Tirler, Peter Dalgaard, Wolfgang Hormann, Josef Leydold, An error in the Kinderman-Ramage method and how to fix it, Computational Statistics & Data Analysis, Volume 47, Issue 3, 1 October 2004, Pages 433-440, ISSN 0167-9473, DOI: 10.1016/j.csda.2003.11.019.
 * (http://www.sciencedirect.com/science/article/B6V8V-4BCWWKJ-1/2/ea68609e7439d02f398f69945a6b6ec8)
 */
final class RandomsGaussKindermanRamage(T = float) : FiberGenerator!(T) {
private:
	Generator!(T) G;
	const T mi, rho;

public:
	/// generic constructor
	this(Generator!(T) G_, T mi_ = 0.0, T rho_ = 1.0) {
		G = G_;
		mi = mi_;
		rho = rho_;
		super();
	}

	/// Automatically determine new uniform random number generator and initialize it to some unpredictable seed
	this(T mi_ = 0.0, T rho_ = 1.0) {
		this(new RandomsUniform!(T)(), mi_, rho_);
	}

	/// Automatically determine new uniform random number generator and initlize it using seed
	this(uint seed, T mi_ = 0.0, T rho_ = 1.0) {
		this(new RandomsUniform!(T)(seed), mi_, rho_);
	}

private:
	/// Just apply transformation N(0,1) -> N(mu,rho)
	override void yield(float x) {
		super.yield(mi + rho*x);
	}

protected:
	/// implementation
	override void iter() {

		static immutable eta = 2.216_035_867_1;

		static double f(double t) {
			auto temp = eta - abs(t);
			if (temp < 0.0) temp = 0.0;
			return normal_pdf(t) - 0.180_025_191_068_563*temp;
		}

		//immutable delta = 0.479_727_404_222_441;
		//immutable k = 0.053_377_549_506_886;
		//immutable gamma = 0.115_779_733_793_499_0;

		immutable half_eta_eta = 0.5*eta*eta;

		while (true) {
			auto u = G.getNext();
			debug assert(0.0 < u && u < 1.0, "Generator G should generate (0;1)");

			if (u < 0.884_070_402_298_758) { // main triangle
				auto v = G.getNext();
				yield(eta*(1.131_131_635_444_180*u + v - 1.0));
			} else if (!(u < 0.973_310_954_173_898)) { // tail
				double t, v, w;
				do {
					v = G.getNext();
					w = G.getNext();
					t = half_eta_eta - log(w);
				} while (v*v*t > half_eta_eta);
				if (u < 0.986_655_477_086_949) {
					yield(sqrt(2.0*t));
				} else {
					yield(-sqrt(2.0*t));
				}
			} else { // thre small triangles (in differential distribution)
				double z, t;
				if (!(u < 0.958_720_824_790_463)) {
					do {
						auto v = G.getNext();
						auto w = G.getNext();
						z = v-w;
						auto min_vw = (v < w ? v : w);
						//auto max_vw = (v > w ? v : w);
						t = eta - 0.630_834_801_921_960*min_vw;
						//if (max_vw <= 0.755_591_531_667_601) break; // quick acceptance
						if (0.034_240_503_750_111*abs(z) <= f(t)) break;  // accurate acceptance
					} while (true); // reject and repeat
				} else if (!(u < 0.911_312_780_288_703)) {
					do {
						auto v = G.getNext();
						auto w = G.getNext();
						z = v-w;
						auto min_vw = (v < w ? v : w);
						//auto max_vw = (v > w ? v : w);
						t = 0.479_727_404_222_441 + 1.105_473_661_022_070*min_vw;
						//if (max_vw <= 0.872_834_976_671_790) break; // quick acceptance
						if (0.049_264_496_373_128*abs(z) <= f(t)) break;  // accurate acceptance
					} while (true); // reject and repeat
				} else {
					do {
						auto v = G.getNext();
						auto w = G.getNext();
						z = v-w;
						auto min_vw = (v < w ? v : w);
						//auto max_vw = (v > w ? v : w);
						t = 0.479_727_404_222_441 - 0.595_507_138_015_940*min_vw;
						if (t < 0.0) {
							continue; // fix by Gunter Tirler, Peter Dalgaard, et al., 2004
						}
						//if (max_vw <= 0.805_577_924_423_817) break; // quick acceptance
						if (0.053_377_549_506_886*abs(z) <= f(t)) break; // accurate acceptance
					} while (true); // reject and repeat
				}

				if (z < 0.0) {
					yield(t);
				} else {
					yield(-t);
				}
			}
		}
	}
}

/** This on avarage consumes 1.01 samples to produce very accurate
 * floats near 0.
 */
T precise_float_gen(T)(Generator!(T) g) {
	assert(0);
}


/** Generetes normal variate from tail (x > b)
 *
 * It iterativly uses generator g.
 * In most cases it will use 2 or 4 samples.
 *
 * Note: normal scaling of uniform integer generator to float by multiplying by 2^-32,
 * is not very accurate because, for tails we want high precision
 * near 0 (for log or div singularity).
 * Without some transofmation on this floats, it can't be smaller than 2^-32,
 * which means that we will not generate z greater than 6.2!
 * Using double still will only allow 9.1.
 *
 * TODO: This method uses additional samples for small values from g,
 * so all posible representable float values are possibly generated (with proper probabilities)
 *
 */
T normal_gen_tail(T)(double r, Generator!(T) g) {
	immutable minus_inv_r = -1.0/r;
	T x, y;
	do {
		// x = r - log(UNI()); // exponential
		x = minus_inv_r*log(g.getNext()); // normal
		y = -log(g.getNext());
	} while (y+y < x*x);
	return (r+x);
}

/*
 * Inne:
 *  Ahrens'a-Dieter'a (tzw. trapezoid)
 *  inna metoda Ahrens'a-Dieter'a
 *  podobna metoda Brent'a
 *
 * Ahrens, Dieter, Efficient table-free sampling methods for the exponential, Cauchy, and normal distributions
 * i rozne referencje
 * J. H. AHRENS AND U. DIETER, "Computer Methods for Sampling from the Exponential and Normal Distributions," Commun. ACM15, 873-882 (1972)
 *
 * Hörmann, W. and Derflinger, G., 1990. The ACR method for generating normal random variables. OR Spektrum 12 3, pp. 181–185. MathSciNet
 *
 * Porownanie roznych metod, ale lekko nieobiektywne :) bo K-R sami jedna wymyslili :) i dosyc stare.
 * Computer Generation of Normal Random Variables, A. J. Kinderman and J. G. Ramage, Journal of the American Statistical Association, Vol. 71, No. 356 (Dec., 1976), pp. 893-896
 */




/// Normal pdf
double normal_pdf(double t) {
	static const coeff = 1.0/sqrt(2.0*PI);
	return coeff * exp(-0.5*t*t);
}

/** Cumulative density function for Normal distribution
 *
 * This function returns
 *   \int_{-\infty}^x {exp(-x^2 / 2) \over \sqrt{2\pi} } dx  == \frac{1}{2} ( 1 + erf(x \over \sqrt{2}) )
 *
 * See, West, G (2009) "Better approximations to cumulative normal functions", Wilmott Magazine, July, 70–76 http://www.wilmott.com/pdfs/090721_west.pdf
 * and refrences.
 *
 * This is invented by Hart (1968).    Hart, J.F. et al, 'Computer Approximations', Wiley 1968, Algorithm 5666
 *
 * This function is accurate to double precision of whole real line. (? i don't think so).
 * For better anything beyond 7 or -7 use tail funcions.
 *
 * Note: This function is far more accurate for negative arguments
 *       (when value is near 0.0), than for positive (when value is near 1.0).
 */
double normal_cdf(double x) {
	double xa = abs(x);
	double c;

	if (xa >= 7.0) {
		c = 0.0;
	} else {
		auto e = exp(-0.5*xa*xa);
		if (xa < 7.07106_78118_6547) {
			auto b = 3.52624_96599_8911e-2 * xa + 0.70038_30644_43688;
			b = b * xa + 6.37396_22035_3165;
			b = b * xa + 33.91286_60783_83;
			b = b * xa + 112.07929_14978_71;
			b = b * xa + 221.21359_61699_31;
			b = b * xa + 220.20686_79123_76;
			c = e * b;
			b = 8.83883_47648_3184e-2 * xa + 1.75566_71631_8264;
			b = b * xa + 16.06417_75792_07;
			b = b * xa + 86.78073_22029_461;
			b = b * xa + 296.56424_87796_74;
			b = b * xa + 637.33363_33788_31;
			b = b * xa + 793.82651_25199_48;
			b = b * xa + 440.41373_58247_52;
			c = c / b;
		} else {
			auto b = xa + 0.65;
			b = xa + 4.0 / b;
			b = xa + 3.0 / b;
			b = xa + 2.0 / b;
			b = xa + 1.0 / b;
			c = e / c / 2.50662_82746_31;
		}
	}

	return (x > 0 ? 1-c : c);
}
unittest {
	// bounduary values
	assert(normal_cdf(-1000.0) == 0.0);
	assert(normal_cdf(-40.0) == 0.0);
	assert(normal_cdf(-10.0) == 0.0);
	assert(normal_cdf(-7.0) == 0.0);
	assert(normal_cdf(-6.8) > 0.0);
	assert(normal_cdf(0.0) == 0.5);
	assert(normal_cdf(6.0) < 1.0);
	assert(normal_cdf(6.8) < 1.0);
	assert(normal_cdf(8.0) == 1.0);
	assert(normal_cdf(40.0) == 1.0);
	assert(normal_cdf(+1000.0) == 1.0);

	// monotonicity
	double prev = 0.0;
	for (double z = -2.0; z <= 2.0; z += 0.004) {
		auto x = normal_cdf(z);
		assert(x > prev);
		assert(x < 1.0);
		prev = x;
	}
	assert(prev < 1.0);
}

import std.stdio;

/** (Numerical) Inversion of gausian CDF
 *
 * Uses initial approximation and few Netwon's method interations.
 *
 * Note: this was using iterative bisection method, but was slow.
 *
 * TODO: use Halley's method which have faster convergance
 * TODO: it is possible to use direct approximation of ICDF by series or ratios.
 *
 * Note/TOOD: It would be far more accurate if this method would
 * be broken into two parts, near 0.0 and near 1.0, because of
 *
 * TODO: in Abramowitz and Stegun, Handbook of Mathematical Functions, Dover, 9th printing, formula 26.2.23 on page 933, is a direct approximation. (but is accurate only to 5-6 digits)
 */
double normal_icdf(double z_) {
	assert(z_ >= 0.0);
	assert(z_ <= 1.0);
	double f(double x) { // growing function
		return (normal_cdf(x) - z_);
	}

version (normal_icdf_bisection) {
	//static const double fleft0 = f(-7.0);
	//static const double fright = f(7.0);

	//double fleft = fleft0 // compile constant or static constant another
	//double fright = fright0;

	double left = -7.0;
	double right = 7.0;
	double half;
	do {
		half = left + 0.5*(right-left); // more accurate
		if (left >= half || half <= right) {
			break;
		}
		double fhalf = f(half);
		//writefln("l=%.25g h=%.25g r=%.25g f(h)=%.20g", left, half, right, fhalf);
		if (fhalf > 0.0) { // fleft is generally always negative
			right = half;
			//fright = fhalf;
		} else if (fhalf < 0.0) { // fritght is always positive
			left = half;
			//fleft = fhalf;
		} else {
			return half;
		}
	} while (true);

	return half;
} else {

	double z; // convert quantile to erfc^{-1} argument

	// we will be using left tail of CDF because it is much more accurate

	if (z_ < 0.5) {
		z = z_;
	} else {
		z = (1.0 - z_);
	}

	static double aprox_erf_inv(double z) { // aproximate inversion of error function
		static immutable double a = 0.14001228868666649; // 8.0*(PI-3.0)/(3.0*PI*(4.0-PI))
		static immutable double inva = 1.0/a;
		static immutable double d = 2.0*inva/PI;
		immutable double l = log(1.0 - z*z);
		immutable double k = (d + 0.5*l);
		immutable double x = sqrt(sqrt(k*k - inva*l) - k);
//writefln("erf_inv(%g)=%g", z, x);
		return x;
	}
	static double aprox_probit(double z) { // quantile function
		return sqrt(2.0)*aprox_erf_inv(2.0*z - 1.0);
	}

	// starting point
	double x = -aprox_probit(z);
//assert(x <= 0.0);
//writefln("z_=%g z=%g x0=%g cdf=%g cdf-z=%g", z_, z, x, normal_cdf(x), (normal_cdf(x) - z));
//assert(normal_cdf(x) <= 0.5);

	double xprev = x;
	double xprevprev = x;
	for (int i = 0; i < 4; i++) {
		auto fx = (normal_cdf(x) - z);
		if (fx == 0.0) {
			break;
		}
		x -= fx/normal_pdf(x);
		if (x == xprev || x == xprevprev) {
			break;
		}
		xprevprev = xprev;
		xprev = x;
//writefln("z_=%g z=%g x=%.26g cdf=%g, cdf-z=%.26g", z_, z, x, normal_cdf(x), (normal_cdf(x) - z));
	}

	if (z_ > 0.5) {
		return -x;
	} else {
		return x;
	}
}
}

/** Elimination of von Neuman
 *
 *
 * 1. generate (x_i, y_i) \in [a,b] \cross [0, sup(f)]
 *
 * 2. accept, and output x_i, only if y_i <= f(x_i)
 *
 * 3. else generate next
 *
 * Note: bad for sharp distributions, because there will be very big fraction (99.9% or more)
 *       of rejected points
 */
final class RandomsVonNeuman(T = float) : FiberGenerator!(T) {
private:
	Generator!(T) X;
	Generator!(T) Y;
	T delegate(T) pdf;

public:
	/** generates random variates described with probability densitity funtion pdf, on interval [a, b]
	 *
	 * pdf doesn't need to be normalized
	 *
	 * up is max pdf(x) \for x \in [a, b]
	 */
	this(Randoms!(uint) G, T delegate(T) pdf_, T a, T b, T up) {
		pdf = pdf_;
		X = new RandomsUniform!(T)(G, a, b);
		Y = new RandomsUniform!(T)(G, 0.0, up);
	}

	ulong discarded; /// number of discarded points

protected:
	/// implementation
	override void iter() {
		while (true) {
			auto x2 = X.getNext();
			auto y2 = Y.getNext();
			if (y2 < pdf(x2)) {
				yield(x2);
			} else {
				discarded++;
			}
		}
	}
}


/*
Reversed distributant:
 p(x) = 1/x 1/(1+x^2)
 u = F(x) = 1/x \int_0^x dt 1/(1+t^2) = 1/2 + arctan(x)/\pi
 x = tg(\pi (u-1/2))

Note: analytic inversion distribuant needed
*/


/** Generated random number described by inversed comulative density function icfg
 *
 *
 * Note: If possible use this function
 */
final class RandomsInverseCDF(T = float) : FiberGenerator!(T) {
private:
	Generator!(T) G;
	T delegate(T) icdf;

public:
	/// constructor
	this(Generator!(T) G_, T delegate(T) icdf_) {
		G = G_;
		icdf = icdf_;
	}

protected:
	/// implementation
	override void iter() {
		foreach (T z; G) {
			auto x = icdf(z);
			yield(x);
		}
	}
}

/** Generates random numbers descirbed by comulative density function cdf
 * by numerical inversion.
 *
 * Note: If possible use first RandomsInverseCDF, if you don't have inversion of CDF,
 * but have CDF (comulative density function),
 * this function can find numerically inversion for you
 */
final class RandomsCDF(T = float) : FiberGenerator!(T) {
private:
	Generator!(T) G;
	T delegate(T) cdf;
public:
	/// constructor
	this(Generator!(T) G_, T delegate(T) cdf_) {
		G = G_;
		cdf = cdf_;
	}
protected:
	/// implementation
	override void iter() {
		foreach (T x; G) {
			auto y = icdf(x);
			yield(y);
		}
	}
private:
	T icdf(T x) {
		return cdf(x);
	}
}

/** Generates random numbers descirbed by probability density function pdf
 * using numerical integration and inversion.
 *
 * Note: If you even don't have inversed CDF or CDF at all, you can use
 * PDF (probability density function), to calculate CDF and it's inversion
 * quickly.
 *
 * http://doi.acm.org/10.1145/945511.945517
 * "Continuous random variate generation by fast numerical inversion"
 *
 * http://epub.wu.ac.at/dyn/virlib/wp/eng/mediate/epub-wu-01_f41.pdf?ID=epub-wu-01_f41
 * "Random Variate Generation by Numerical Inversion when Only The Density is Known"
 */
final class RandomsPDF(T = float) : FiberGenerator!(T) {
private:
	Generator!(T) G;
	T delegate(T) pdf;

public:
	/// implementation
	this(Generator!(T) G_, T delegate(T) pdf_) {
		G = G_;
		pdf = pdf_;
	}

protected:
	/// implementation

	override void iter() {
		foreach (T x; G) {
			auto y = pdf(x);
			yield(y);
		}
	}
}

// VoigtProfile
// http://en.wikipedia.org/wiki/Voigt_profile


/+
B^{-1} = C = E{(X-mi)*(X-mi)^T }

covs(X,Y) = zaleznosc liniowa pomiedzy zmiennymi losowamy X i Y

rozkład
\phi(X) = k * exp( - (x-mi)^T B (x-mi) ) // chyba)

k = {1 \over (2\pi)^{d/2}  \det(B)^{1/d}} 


C = B^{-1}

d = 2

k = 1/(2 pi sqrt(detB))

C = ( var(x)   cov(x,y) )
    ( cov(x,y) var(y)   )

B = (var(x)var(y) - cov(x,y)^2)^-1 (

Jeśli zmienne są niezależne, to cov(x,y) == 0.

B_0 = (var(x)var(y))^{-1} = 

1. tworzymy wektor Z z n niezlaeżymi zmiennymi standardowym rozkladzie normalnym
2. rozklad choleskiego, AA^T = C
3. szukamy X = /mi + AZ
+/


/** Multidimensional Gauss generator.
 *
 * It will yield reference to vectors of variates.
 * One can change yielded value, but, it can change after next yield, so dont
 */
final class RandomGaussMultivariate(n, T = float) : FiberGenerator!(vector!(T)) {
private:
	Generator!(T) G;
	matrix A; // Cholesky decomposition of C, A is lower triangular matrix
	vec mi;
	vec z;
	vec x;

public:
	/** Creates multivariate normal distribution with parameters given by:
	 *
	 *  Params: covariances: symetric positivly defined matrix C,
	 *          means: given in vector mi
	 *
	 * Note: Generator G, shoule generate standard gauss distribution (with \mi=0 and \sigma=1)
	 *
	 * Remark: We are using Cholesky decomposition of matrix C
	 */
	this(Generator!(T) G_, vector!(T) mi_, matrix!(T) C) {
		assert(C.isSquare());
		assert(n == C.rows && n == C.cols);
		assert(n == mi.rows);
		G = G_;
		cholesky_decomp!(T)(C, A);
		x = new vector!(T)(n);
		z = new vector!(T)(n);
		mi = mi_;
	}

protected:
	/// implementation
	override void iter() {
		while (true) {
			// we are using separate loops, to have better cache locality
			foreach (i; 1 .. n) {
				z[i] = G.getNext();
			}
			foreach (i; 1 .. n) {
				x[i] = mi[i];
				foreach (j; 1 .. i) { // A is lower triangular
					x[i] += A[i,j]*z[j];
				}
			}
			foreach (i; 1 .. n) {
				yield(x);
			}
		}
	}
}


/**
 * Needed:

Discrete:
Bernoulli
Binomial
Geometric
Negative binomial

Continious:
Poisson
Exponential
Gamma
Weibull
Extreme value distribution
Normal
Log-normal
Chi squered
Cauchy
Fisher F
Student t

generators: access to real random numbers (like /dev/[u]random, openssl, hardware RNG).
 */

/** Other distributions:
 *
 * Cauchy   pdf(x) = 1/(1 + x^2)/\pi
 * Normal (Gauss) pdf(x) = 1/\sqrt(2\pi\rho^2)  \exp( - (x-x_0)^2 / \rho^2 )
 * Poison   pmf(k) = \lambda^k \exp(-\lambda) \over k!
 * Gamma    pdf(x) = x^{k-1} \exp(-x/\theta) \over \theta^k \Gamma(k), for x >= 0, k > 0, \theta > 0
            if k is integer, Gamma(k) = (k-1)!, and it is known as Erlang distribution: pdf(x) = \lambda^k x^{k-1} \exp(-\lambda x) \over (k-1)!
 * beta     pdf(x) = x^{\alpha-1} (1-x)^{\beta-1} / B(\alpha, \beta); where B = \Gamma(\alpha+\beta) \over \Gamma(\alpha)\Gamma(\beta)
 * Weibull  pdf(x) = k/\lambda (x / \lambda)^{k-1} \exp( - (x/\lambda)^k ), for x >= 0, k > 0 (shape), \lambda > 0 (scale)
 *          for k = 1, it is exponential distribution
 *          for k = 2, it is Rayleight distribution
 *          cdf(x) = 1 - \exp( - (x/\lambda)^k )
 *          There is also generalization known as "exponentiated Weibull dist." with cdf(x)^\alpha, \alpha > 0 (second shape)
 * Gumbel   pdf(x) = \exp(x' - \exp( x'/\beta)) / \beta, where x' = (x-\alpha)/\beta;  cdf(x) = 1 - exp(-exp(x')); simple inverse
 * Exponential
 */

