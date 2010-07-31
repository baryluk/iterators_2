module autocorel;

import std.stdio : writef, writeln, writefln;

/** Circular buffer which is always fully populated,
 */
final class CircularBuffer(T, uint n) {
private:
	T[n] c;
	uint m;
	uint allm;

	static assert(n > 0);

	invariant() {
		assert(0 <= m && m < n);
/++		if (m == 0) {
			assert(opIndex(0) == c[n-1]);
			assert(opIndex(-n+1) == c[0]);
		} else {
			assert(opIndex(0) == c[m]);
			assert(opIndex(-n+1) == c[n-m]);
		}
		foreach (i, y; this) {
			assert(y == opIndex(i));
		}
++/
	}

public:
	/** '~=' operations, removes oldest entry in buffer, and adds y as new
	 *
	 * Note: this is O(1)
	 */
	void opCatAssign(T y)
	out {
		assert(y == opIndex(0));
	}
	body {
		c[m] = y;
		m = (m + 1) % n;
		allm++;
	}

	/// i is negative number, describing how last data access
	T opIndex(int i) {
		assert(i <= 0 && -i < n);
		int k = (m+i-1);
		static assert(is(typeof(k) == int));
		assert(k < cast(int)n); // see dmd bug #259
		if (k < 0) {
			k += n;
			assert(k < n);
			return c[k];
		} else {
			return c[k];
		}
	}

	///
	int opApply(int delegate(ref int, ref T) dg) {
		int j = -n+1;
		foreach (i; m .. n) {
			auto ci = c[i];
			auto r = dg(j, ci);
			j++;
			if (r) return r;
		}
		foreach (i; 0 .. m) {
			auto ci = c[i];
			auto r = dg(j, ci);
			j++;
			if (r) return r;
		}
		return 0;
	}

	///
	void print() {
		foreach (ref int j, ref T y; this) {
			writef("c[%d] = %5.2f\t", j, y);
		}
		writeln();
	}

	unittest {
/++
		auto c = new CircularBuffer!(double, 7)();
		c.print();
		c ~= 4.1;
		c.print();
		c ~= 1.3;
		c.print();
		c.print();
		c ~= 6.6;
		c.print();
		c ~= 4.3;
		c.print();
		c ~= 7.8;
		c.print();
		c ~= 1.1;
		c.print();
		c ~= 6.1;
		c.print();
		c ~= 9.8;
		c.print();
		c ~= 7.7;
		c.print();
		c ~= 16.7;
		c.print();
		c ~= 51.6;
		c.print();
		c ~= 11.5;
		c.print();
		c ~= 57.2;
		c.print();
++/
	}
}


/** Generally it is just alias to CircularBuffer, it remembers last N
 * opCatAssign'ed values. Uses CircularBuffer as implementation.
 */
class RememberLastN(T, int n) {
private:
	CircularBuffer!(T, n) c;
public:
	this() {
		c = new CircularBuffer!(T, n)();
	}

	void opCatAssign(T y) {
		c ~= y;
	}

	T opIndex(int i) {
		return c[i];
	}
}

import misc : kahan_update;

version = use_kahan;

/** Calculates autocorelattion function
 *
 * C_k = 1/(N-K) \sum_{i=0}^{N-k} y_i * y_{i+k}
 *
 * for k = 0, 1, ... n
 *
 * y should be difference from mean, so it this form it is only usefull
 * for apriori known mean, and if it is 0.0
 */
class Autocorellation(T, int kmax) {
private:
	T[kmax+1] a;
version (use_kahan) {
	T[kmax+1] a_errs;
}
	RememberLastN!(T, kmax+1) l;
	int N;
	int kmax2;

public:
	static const max = kmax;
	this() {
		a[] = 0.0;
		l = new RememberLastN!(T, kmax+1)();
		kmax2 = kmax;
		version (use_kahan) {
			a_errs[] = 0.0;
		}
	}

	///
	void opCatAssign(T y) {
		N++;
		l ~= y;
		assert(y == l[0]);
		if (N > kmax) { // if RememberLastN is full
			foreach (j, yj; l.c) {
				// j is negative number
				version (use_kahan) {
					kahan_update!(T)(a[-j], y*yj, a_errs[-j]);
				} else {
					a[-j] += y*yj;
				}
			}
/+
			foreach (k; 1 .. kmax) {
				a[k] += y*l[k];
			}
+/
		} else { // else partially full (kmax > N, do only to N
			foreach (k; 0 .. N) {
				version (use_kahan) {
					kahan_update!(T)(a[k], y*l[-k], a_errs[k]);
				} else {
					a[k] += y*l[-k];
				}

			}
		}
	}

	/// value of C_k
	T opIndex(int k) {
		version (use_kahan) {
			return (a[k]+a_errs[k])/(N-k);
		} else {
			return a[k]/(N-k);
		}
	}

	/// value of C_k/C_0
	T normedIndex(int k) {
		version (use_kahan) {
			return (((a[k]+a_errs[k])/(a[0]+a_errs[0]))*N)/(N-k);
		} else {
			return ((a[k]/a[0])*N)/(N-k);
		}
	}

	/* TODO: - <x_i>^2, biased
	 * TODO: - <x_i><x_{i-k}>, unbiased
	 */
}

