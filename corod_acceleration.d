module corod_acceleration;

import corod; /* FiberGenerator, Generator; */


/** Aitken's delta-squared process
 *
 * Ax = { x_{n+2}x_n - x_{n+1}^2 \over x_{n+2} - 2 x_{n+1} + x_n }
 *
 * Which can be also written as:
 *
 * Ax_n = x_n - { (\Delta x_n)^2 \over \Delta^2 x_n }
 *
 * Where \Delta x_n = x_{n+1} - x_n
 *       \Delta^2 x_n = x_n - 2 x_{n+1} + x_{n+2}
 *
 * for n = 0, 1, 2, 3, ...
 *
 * Keep attention if \Delta^2 x is near zero (sequence have repeating terms,
 * or difference sequence have repeapting terms)
 *
 * It is probably oldest, and simplest nonlinear sequence transformation.
 * It is(?) first of the Shanksen's transform.
 */
class Aitken(T = float) : FiberGenerator!(T) {
private:
	Generator!(T) g;
public:
	///
	this(Generator!(T) g_) {
		g = g_;
	}
protected:
	void iter() {
		auto x_n0 = g.getNext();
		auto x_n1 = g.getNext();
		auto delta_x_n_10 = x_n1 - x_n0;

		while (true) {
			auto x_n2 = g.getNext();
			auto delta_x_n_21 = x_n2 - x_n1;

			if (delta_x_n_21 != delta_x_n_10) {
				auto Ax_n = x_n0 - (delta_x_n_10*delta_x_n_10 / (delta_x_n_21 - delta_x_n_10));
				yield(Ax_n);
			} else {
				throw new Exception("repeting difference terms");
			}

			x_n0 = x_n1;
			x_n1 = x_n2;
			delta_x_n_10 = delta_x_n_21;
		}
	}
}

// todo: Euler transform for alternating series
// todo: Van_Wijngaarden_transformation for alternating series
// todo: e-algorithm by Pettern Wynn
// todo: Levin u-tranform
// todo: Wilf-Zeilberger-Ekhan method
// todo: WZ method
// todo: alternating series: http://www.math.utexas.edu/~villegas/publications/conv-accel.pdf
// np. Pade type approximation, bardzo prosta i szybka (5.282^-n), lepsza i nie wymagajaca dodatkowej pamieci w pornwniu do van Wijngaardena czy Eulera
// sa tez triki aby zastosowac te algorytmy do niealternaujacych szeregow i ulamkow lancuchowych

class EpsilonAlgorithm(T = float) : FiberGenerator!(T) {
	
}

import std.math : sqrt, pow;

/** Algorithm 1 from 'Convergence Acceleratiorn of Alternating Series", Henri Cohen, Fernando Rodriguez Villegas, Don Zagier
 *
 *  A = \sum_{k=0}^\infty (-1)^k a_k,   a_k > 0
 *
 * TODO: Change T to gmp, and long to gmz
 *
 * It is simple, robust and fast linear (matrix) tranformation for oscylating series (alternating elements)
 */
T alg1(T = double, LONG = long)(Generator!(T) A, LONG n) {
	static if (is(T == double) || is(T == double)) {
		static if(is(LONG == long)) {
		if (n > 20) {
			throw new Exception("n to big for representation intermidiate d in long variable");
		}
		}
	}

	alias LONG BIGLONG;

	auto d = (3 + sqrt(cast(T)8));
	d = pow(d, cast(uint)n);
	d = (d + 1/d)/2;

	auto di = cast(BIGLONG)(d+0.5);

	BIGLONG b = -1;
	BIGLONG c = -di; // d should be integer, but it only works properly up to n == 21
	T s = cast(T)(0);

	foreach (long k; 0 .. n) {
		c = b - c;
		auto ak = A.getNext();
		s += c*ak;
		//writefln("c=%d  b=%d  d=%d", c, b, di);
		b *= 2*(k+n)*(k-n);
		auto b2 = b;
		b /= ( (2*k + 1) * (k+1) );
		assert(b * (2*k + 1) * (k+1) == b2);
	}

	return s/di;
}

class PI_terms(T = float) : FiberGenerator!(T) {
protected:
	void iter() {
		T a = 0.0;
		int n21 = 1;
		while (true) {
			a += 1.0/cast(T)(n21);
			n21 += 2;
			yield(4.0*a);
			a += -1.0/cast(T)(n21);
			n21 += 2;
			yield(4.0*a);
		}
	}
}

class PI_serie(T = float) : FiberGenerator!(T) {
protected:
	void iter() {
		int n21 = 1;
		while (true) {
			T a = 1.0/cast(T)(n21);
			n21 += 2;
			yield(4.0*a);
		}
	}
}

version(acceleration_main_test) {

import std.stdio : writefln;

alias double T2;

alias PI_terms!(T2) PId;
alias Aitken!(T2) Ad;

void main() {
	auto pi = new PId();
	auto pi_aitken = new Ad(new PId());
	auto pi_aitken2 = new Ad(new Ad(new PId()));
	auto pi_aitken3 = new Ad(new Ad(new Ad(new PId())));

	auto mpi = 3.1415926535897931;

	int i;
	foreach (p; new MultiForeach!(T2, 4)([cast(Generator!(T2))pi, pi_aitken, pi_aitken2, pi_aitken3])) {
		writefln("%d\t%-011.6g\t%-011.6g\t%-11.6g\t%-011.6g", i++, p[0]-mpi, p[1]-mpi, p[2]-mpi, p[3]-mpi);
		if (i > 100) break;
	}

	void gopi(int n) {
		auto pi_serie = new PI_serie!(T2)();
		writefln("%d\t%.15f", n, alg1!(T2,long)(pi_serie, n));
	}
	gopi(1);
	gopi(2);
	gopi(3);
	gopi(4);
	gopi(5);
	gopi(6);
	gopi(8);
	gopi(10);
	gopi(12);
	gopi(14);
	gopi(16);
	gopi(18);
	gopi(20);
	//gopi(21);

	T2[10] xx;

	fill!(T2)(new PI_serie!(T2)(), xx);

	writefln("%s", xx);

	fill!(T2)(new PId(), xx);
	writefln("%s", xx);

	fill!(T2)(new Ad(new Ad(new PId())), xx);
	writefln("%s", xx);

}
}
