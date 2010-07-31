module hist;

import std.stdio : writefln, writef, writeln;
import std.math : sqrt;

/** Equi-depth histogram, with predefined number and positions of bins
 *
 *
 * TODO: 2D histogram.
 */
final class Histogram(T = float) {
private:
	T a;
	T b;
	T c;
	uint n; // number of points in the boxes (doesn't include below and above)
	uint boxes;
	uint[] box;
	uint below;
	uint above;

	invariant() {
		assert(a < b);
		assert(box.length == boxes);
/+
		uint s;
		foreach (b; box) {
			assert(b <= n);
			s += b;
		}
		assert(s == n);
+/
	}

public:
	/** Create histogram for counting values in interval [a; b)
	 * divided into #boxes subintervals (bins) of the same width.
	 */
	this(T a_, T b_, uint boxes_) {
		a = a_;
		b = b_;
		boxes = boxes_;
		box.length = boxes;
		box[] = 0;
		c = 1.0/(b-a)*boxes;
		n = 0;
	}

	/** Adds point x to histogram
	 *
	 * Note: Points outside of [a; b) are added to variables above or below.
	 */
	void opCatAssign(T x) {
		if (a <= x) {
			if (x < b) {
				//const i = cast(uint)((x-a)/(b-a) * boxes);
				const i = cast(uint)((x-a)*c);
				assert(0 <= i && i < box.length);
				box[i]++;
				n++;
			} else {
				above++;
			}
		} else {
			below++;
		}
	}

	/// Remove point x from histogram.
	void opSubAssign(T x) {
		if (a <= x) {
			if (x < b) {
				const i = cast(uint)((x-a)*c);
				assert(0 <= i && i < box.length);
				assert(box[i] > 0);
				box[i]--;
				n--;
			} else {
				assert(above > 0);
				above--;
			}
		} else {
			assert(below > 0);
			below--;
		}
	}


	/// TODO: implement
	void opCatAssign(Histogram!(T) h2) {
		throw new Error("not implemnted");
	}

	/// TODO: implement
	void opSubAssign(Histogram!(T) h2) {
		throw new Error("not implemnted");
	}

	/** Returns histogram being combined result of two histograms.
	 *
	 * TODO: implement
	 */
	Histogram!(T) opCat(Histogram!(T) h2) const {
		throw new Error("not implemnted");
	}

	/// TODO: implement
	Histogram!(T) opSub(Histogram!(T) h2) const {
		throw new Error("not implemnted");
	}


	/** Returns number of events in bin assoscied with variable x.
	 *
	 * Note: For values x outside of [a; b) this method returns 0.
	 */
	uint opIndex(T x) const {
		if (a <= x && x < b) {
			const i = cast(uint)((x-a)*c);
			return box[i];
		} else {
			return 0;
		}
	}

	/** Sets number of events in bin assoscied with variable x
	 *
	 *
	 * Note: Values x outside of [a; b) are ignored;
	 */
	void opIndexAssign(uint v, T x) {
		if (a <= x && x < b) {
			const i = cast(uint)((x-a)*c);
			box[i] = v;
		}
	}

	/** Prints to standard output state of histogram:
	 *
	 * Format is following:
	 *   "boxes" lines in folowing format (columns separated with space and ended with new line):
	 * -----
	 *     i lower upper half b_i b_i_n b_i_n_c
	 * -----
	 * where:
	 *   i - index of bin, 0...boxes-1
	 *   lower, upper - lower and upper bound of bin [lower, upper)
	 *   half - midpoint of bin
	 *   b_i - number of events for this bin
	 *   b_i_n - number of events for this bin divided by n (number of events)
	 *   b_i_n_c - b_i_n multiplied by width of bin (normalized histogram value)
	 */
	void dump() const {
		const T alln = n+below+above;
		foreach (i, b_i; box) {
			const lower = i/c+a;
			const upper = (i+1)/c+a;
			const half = 0.5*(lower+upper);
			writefln("%d %g %g %g %d %g %.10g", i, lower, upper, half, b_i, b_i/alln, b_i/alln*c);
		}
	}

	/// return normalized value of histogram at bin b_i
	T normed_bin(uint i) const {
		const T alln = n+below+above;
		auto b_i = box[i];
		return c*b_i/alln;
	}

	
	///
	void consume(Generator!(T) g) {
		foreach (T x; g) {
			opCatAssign(x);
		}
	}

	/** This methods calculates mean and variance solaly based on the points stored in the histogram
	 *
	 * Note: It is much better to use Moments or CentralMoments class for numerically stable,
	 *       incremental/parallel algorithms for mean and covarianc, or higher moments.
	 *
	 * Note: Points outside of the histogram are ignore in computations.
	 *
	 * TODO: check numerical stability and accuracy
	 * TODO: if needed use Kahan's update scheme.
	 */
	T mean() const {
		T s = 0.0;
		foreach (i, b_i; box) {
			auto x = (i+0.5)/c+a;
			s += b_i*x;
		}
		return s/n;
	}


	/// ditto
	T variance(T mean0) const {
		T s2 = 0.0;
		foreach (i, b_i; box) {
			auto x = (i+0.5)/c+a;
			auto y = (x - mean0);
			s2 += b_i*y*y;
		}
		return s2/n;
	}

	/** Aproximated value of minimum and maximum point
	 *
	 * If there wear some points added outside of the histogram, respective limits will be returned.
	 *
	 * Remark: Class MinMax is generally more accurat.
	 */
	T min() const {
		if (below) { return a; }

		foreach (i, b_i; box) {
			if (b_i) {
				return i/c+a;
			}
		}

		if (above) { return b; }

		throw new Exception("no points in histogram");
	}

	/// ditto
	T max() const {
		if (above) { return b; }

		foreach_reverse (i, b_i; box) {
			if (b_i) {
				return (i+1)/c+a;
			}
		}

		if (below) { return a; }

		throw new Exception("no points in histogram");
	}

	/** Approximated median.
	 *
	 * Note: Points outside of histogram are ignored.
	 *
	 * TODO: It can be estimated sligthly more accurate using weighting/interpolation inside the bin,
	 *       and assuming uniform distribution (it is quite good assumption for bins narow enaugh).
	 */
	T median() const {
		int s;
		const half_n = n/2;
		foreach (i, b_i; box) {
			s += b_i;
			assert(s <= n);
			if (s >= half_n) {
				return (i+0.5)/c+a; // interpolation.
			}
		}

		throw new Exception("no points in histogram");
	}

	/** Approximated median.
	 *
	 * It is similar but it doesn't (or tries to do so) ignore points outside of histogram.
	 */
	T median_notignore() const {
		uint s = below;
		const uint all_n = below+n+above;
		const uint half_n = all_n/2;

		if (s >= half_n) { return a; } // we cant do more precise than that without knowing distribution in tails, or minmax values and uniformity assumption

		foreach (i, b_i; box) {
			s += b_i;
			assert(s <= n);
			if (s >= half_n) {
				return (i+0.5)/c+a; // here we can do sligthly better, by interpolation
			}
		}

		return b; // similary nothing better

		throw new Exception("no points in histogram");
	}

	/** Approximated location of quantile q.
	 *
	 * This code tries to take into account the fact that some points are outside of the histogram.
	 *
	 * If you you ask for quantile which will be outside of the histogram, then min or max will be returned;
	 */
	T quantile_notignore(float q)
	in {
		assert(0 <= q && q <= 1);
	}
	body {
		uint s = below;
		const uint all_n = below+n+above;
		float inv_all = 1.0/all_n;

		if (s*inv_all >= q) { return a; }

		foreach (i, b_i; box) {
			s += b_i;
			assert(s <= n);
			if (s*inv_all >= q) {
				return (i+0.5)/c+a; // here we can do sligthly better, by interpolation
			}
		}

		return b; // similary nothing better
	}

	/// Value for which there is maximal value of histogram
	T mode() {
		uint i_max = 0;
		T s_max = -1.0/0.0;
		foreach (i, b_i; box) {
			if (b_i >= i_max) {
				i_max = b_i;
				s_max = (i+0.5)/c + a;
			}
		}
		return s_max;
	}

/+
	unittest {
		auto h = new Histogram!(float)();
		assert(h.quantile_notignore(0.0) == h.min());
		assert(h.quantile_notignore(1.0) == h.max());
		assert(h.quantile_notignore(1.0) == h.median_notignore());
		foreach (i; 1 .. 10) {
			assert(h.quantile_notignore(0.1*i) <= h.quantile_notignore(0.1*(i+1)));
		}
		foreach (i; 1 .. 10) {
			assert(h.quantile_notignore(0.1*i) <= h.quantile_notignore(0.11*i+1));
		}
	}
+/
}

/** Keeps minimal and maximal element in stream of values x */
final class MinMax(T = float) {
	float min = 1.0/0.0; /// minimal value
	float max = -1.0/0.0; /// maximal value

	/// Take value x, and update min/max
	void opCatAssign(T x) {
		if (x < min) {
			min = x;
		}
		if (x > max) {
			max = x;
		}
	}
}

/+
Automaticalli adjusted Histograms:

First work, is very interesting, it construct online approximate (1+epsilon)
histogram in polylog(n) space, and runs in linear time.


Data-Streams and Histograms (2001)
by Sudipto Guha ,  Nick Koudas ,  Kyuseok Shim 
http://cs.kaist.ac.kr/~shim/stoc01.ps
http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.29.634

 Fast Incremental Maintenance of Approximate Histograms (1997) [124 citations — 22 self]
by Phillip B. Gibbons ,  Yossi Matias ,  Viswanath Poosala 
http://www.pittsburgh.intel-research.net/people/gibbons/papers/vldb97.pdf
http://www.informatik.uni-trier.de/~ley/db/conf/vldb/vldb97.html#GibbonsMP97
http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.100.9995

 Self-tuning histograms: Building histograms without looking at data (1999) [103 citations — 10 self]
by Ashraf Aboulnaga 
http://wwwiti.cs.uni-magdeburg.de/~eike/selftuning...
http://www.cs.wisc.edu/~ashraf/pubs/sigmod99sthist...
http://www.cs.uwaterloo.ca/~ashraf/pubs/sigmod99st...
http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.134.7127

Fast, Small-Space Algorithms for Approximate Histogram Maintenance (Extended Abstract) (2002) [2 citations — 0 self]
by Anna Gilbert ,  Yannis Kotidis ,  Sudipto Guha ,  S. Muthukrishnan ,  Piotr Indyk ,  Martin J. Strauss
In Proc. 34th ACM Symp. on the Theory of Computing. ACM 
http://www.cis.upenn.edu/~sudipto/mypapers/dynamic...
http://www.cs.princeton.edu/courses/archive/spring...
http://www.cs.utsa.edu/~qitian/seminar/Spring05/02...
http://www.research.att.com/~agilbert/ps.files/sud...
http://www.research.att.com/~agilbert/ps.files/ggi
http://people.csail.mit.edu/indyk/ggikms.pdf...
http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.57.7380
+/

/// Calculates simple moments
final class Moments(T) {
private:
	uint n;
	T s = cast(T)0.0;
	T s2 = cast(T)0.0;

	float Q = 0.0;
	float A = 0.0;

public:
	// http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
	// there also online algorithm
	// remember about Kahan update
	// or http://en.wikipedia.org/wiki/Standard_deviation#Rapid_calculation_methods
	// http://prod.sandia.gov/techlib/access-control.cgi/2008/086212.pdf

	/** Update Moments adding value x. See CentralMoments4 for more informations.
	 */
	void opCatAssign(T x) {
		s += x;
		s2 += x*x;

		n++;

		const double n_inv = 1.0/n;

		const double delta = (x - A);

		A += n_inv*delta;
		Q += delta*(x - A);
	}

	/// Like above, but not in-place
	auto opCat(T x) {
		auto NEW = new Moments!(T)();

		NEW.s = s + x;
		NEW.s2 = s + x*x;

		NEW.n = n+1;

		const double n_inv = 1.0/NEW.n;

		const double delta = (x - A);

		NEW.A = A + n_inv*delta;
		NEW.Q = Q + delta*(x - NEW.A);

		return NEW;
	}

	/** Remove value x, as it will apear that it was never added with opCatAssign.
	 *
	 * Don't remove value x, more than it was added. This will produce bogus behvaiour.
	 * See CentralMoments4 for more information
	 */
	void opSubAssign(T x) {
		s -= x;
		s2 -= x*x;

		const old_n = n - 1; // separate old_n to avoid confsions by in-place change

		const double old_n_inv = 1.0/old_n;

		//const old_A = (A*(n+1) - x)/n; // this is correct, but can be numerical instable

		// const old_delta = ((x-A)/old_n)*n; // this is old delta, but excesive operaions, so we put them directly into old_A
		const double sub_delta = (x - A); // this is different delta than in opCatAssign, and we will use it

		//const old_A = A - old_delta*n; // this is correct, but as old_delta have factor "n", an we here remove it, next expression is better
		const double old_A = A - sub_delta*old_n_inv; // this is exactly like line above, but one mul, and one div less, and more accurate

		const old_Q = Q - sub_delta*(x - old_A); // == Q - (x - A')*(x - A)

		n = old_n;
		A = old_A;
		Q = old_Q;
	}

	/// like above, but not in-place
	auto opSub(T x) {
		auto OLD = new Moments!(T)();

		OLD.s = s - x;
		OLD.s2 = s - x*x;

		OLD.n = n - 1;

		const double old_n_inv = 1.0/OLD.n;

		//const old_A = (A*(n+1) - x)/n; // this is correct, but can be numerical instable

		// const old_delta = ((x-A)*old_n_inv)*n; // this is old delta, but with excesive operaions, so we put them directly into old_A
		//OLD.A = A - old_delta/n; // this is correct, but as old_delta have factor "n", an we here remove it, next expression is better

		const double sub_delta = x - A; // this is different delta than in opCatAssign, and we will use it in different way
		OLD.A = A - sub_delta*old_n_inv; // this is exactly like line above, but one mul, and one div less, and more accurate

		OLD.Q = Q - sub_delta*(x - OLD.A); // == Q - (x - A')*(x - A)

		return OLD;
	}

	/// Compose two multisets into one, calculating internal statistics
	Moments!(T) opCat(Moments!(T) BB) {
		auto AA = this;
		Moments!(T) AB = new Moments!(T)();

		AB.n = AA.n + BB.n;
		AB.s = AA.s + BB.s;
		AB.s2 = AA.s2 + BB.s2;

		const double u_delta = BB.A - AA.A;
		const double u_delta_over_n = u_delta/AB.n; // Ab.n or n

		AB.A = (AA.A*AA.n + BB.A*BB.n)/AB.n;
//		AB.A = AA.A + BB.n*u_delta_over_n; // CGL79, it is mathematically equivallent to the AB.A above

		const AAnBBn = AA.n*BB.n;
		AB.Q = AA.Q + BB.Q + AAnBBn*u_delta_over_n*u_delta; // CGL79

		return AB;
	}

	/// Like previous, but in-place
	void opCatAssign(Moments!(T) BB) {
		const AAn = this.n;

		n += BB.n;
		s += BB.s;
		s2 += BB.s2;

		const AAnBBn = AAn*BB.n;

		const u_delta = BB.A - this.A;
		const u_delta_over_n = u_delta/n; // Ab.n or n

		A = (this.A*AAn + BB.A*BB.n)/n; /// TODO: check if this, and next line is equivalent, then AAn is not needed to be remembered)
//		A += BB.n*u_delta_over_n; // CGL79, it is mathematically equivallent to the AB.A above

		Q += BB.Q + AAnBBn*u_delta_over_n*u_delta; // CGL79
	}

	/// Remove multiset BB from this, and return reduced (uncomposed) multiset
	Moments!(T) opSub(Moments!(T) BB) {
		auto AB = this;
		Moments!(T) AA = new Moments!(T)();

		AA.n = AB.n - BB.n;
		AA.s = AB.s - BB.s;
		AA.s2 = AB.s2 - BB.s2;

/+
		// this is correct, but pretty complicated, and not so stable
		AA.A = (AB.n*AB.A - BB.n*BB.A) / AA.n;
		const AAnBBn = AA.n*BB.n;
		AA.Q = - (AA.n*AB.Q - AA.n*AB.Q + (BB.n*BB.n + AAnBBn)*BB.A*BB.A + (-2*BB.n*BB.n - 2*AAnBBn)*AB.A*BB.A + (BB.n*BB.n + AAnBBn)*AB.A*AB.A) / AB.n;
+/

		const double sub_delta = BB.A - AB.A;
		const double sub_delta_over_n = sub_delta/AA.n;

		AA.A = AB.A - sub_delta_over_n*BB.n;

		const double u_delta = BB.A - AA.A; // original delta
		const double u_delta_over_n = u_delta/AB.n; // Ab.n or n

		const AAnBBn = AA.n*BB.n;
		AA.Q = AB.Q - (BB.Q + AAnBBn*u_delta_over_n*u_delta);

		return AA;
	}

	/// Like above, but in-place
	void opSubAssign(Moments!(T) BB) {
		throw new Error("not implemented");
	}

	/// consume all values from generator g
	void consume(Generator!(T) g) {
		foreach (x; g) {
			opCatAssign(x);
		}
	}

	///
	void unconsume(Generator!(T) g) {
		foreach (x; g) {
			opSubAssign(x);
		}
	}

	///
	T mean() const {
		//return s/n;
		return A;
	}

	///
	T variance() const {
		//return (s2 - s*s/n)/(n-1);
		return Q/(n-1);
	}

	///
	T sample_deviation() const {
		return sqrt((Q / n)*(n-1));
	}

	///
	T m2() const {
		return Q;
	}

	///
	T deviation() const {
		return sqrt(variance);
	}

	///
	uint count() const {
		return n;
	}
}

/+

Mean:

u = u_1 + (y-u1)/n

n*u = n*u_1 + y-u1 = (n-1)*u_1 + y = \sum_1^{n-1} y_i + y
u = (\sum_1^{n-1} y_i + y)/n

Second central moment:

M = M_1 + (y-u_1)(y-u)

M = \sum_1^{n-1} (x-u_1)^2 + (y-u_1)(y-u)

u1 = u - (y-u1)/n

u1 = (n*u-y)/(n-1)



M = \sum_1^{n-1} (x- (u - (y-u1)/n))^2 + (y-u_1)(y-u)
= \sum_1^{n-1} ( (x- u) + (y-u1)/n)^2 + (y-u_1)(y-u)
= \sum_1^{n-1} [ (x- u)^2 +  2(x-u)(y-u1)/n + ((y-u1)/n)^2  ] + (y-u_1)(y-u)
= [ \sum_1^{n-1} (x- u)^2 ] + 2(n-1)(u1-u)(y-u1)/n + (n-1)*((y-u1)/n)^2  + (y-u + (y-u1)/n)(y-u)
= [ \sum_1^{n-1} (x- u)^2 ] + 2(n-1)(u1-u)(y-u1)/n + (n-1)*((y-u1)/n)^2  + (y-u)^2  + ((y-u1)/n)(y-u)
= [ \sum_1^n (x- u)^2 ] + 2(n-1)(u1-u)(y-u1)/n + (n-1)*((y-u1)/n)^2  + ((y-u1)/n)(y-u)
... // maxima, inculding substitution u1 = -(y-n*u)/(n-1)
= [ \sum_1^n (x- u)^2 ] + 0
= M

generaly:

u = u' + n'' * (u''-u')/n
M = M' + M'' + n' * n'' * (u''-u')^2 / m

Thrid central moment:

T = T' + T'' + n'*n''*(n'-n'')*(u''-u')^3 / n^2 + 3(n'*M'' - n''*M')*(u''-u')/n

Forth central moment:

F = F' + F'' + n'*n''*(n'^2 - n'*n'' + n''^2)*(u''-u')^4/n^3 + 6(n'^2*M'' + n''^2*M')*(u''-u')^2/n^2 + 4(n'*T'' - n''*T')*(u''-u')/n


\generally
for p >= 1

M_p = M_p' + M_p'' + \sum_{k=1}^{p-2} ( k \choose p ) [ (-n''/n)^k M_{p-k}' + (n'/n)^k M_{p-k}'' ] (u''-u')^k  + [ (n'*n'')/n * (u''-u') ] ^ p * [ n'' ^ {1-p} - (-n')^{1-p} ]

and Covariance:

C = \sum_{u,v} (u - u')(v - v')

C = C' + C'' n'*n'' / n * (u''-u')*(v''-v')

+/

version = paranoicchecks;
//version=opcat_old_version;

/** Stores first 4 central moments of set points in O(1) space, with O(1) update and merge operations.
 *
 * All operations should be stable numerically.
 *
 * See also:
 *
 * Chan, Tony F.; Golub, Gene H.; LeVeque, Randall J. (1979), "Updating Formulae and a Pairwise Algorithm for Computing Sample Variances.", Technical Report STAN-CS-79-773, Department of Computer Science, Stanford University.
 *   ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/79/773/CS-TR-79-773.pdf
 *
 * Terriberry, Timothy B. (2007), Computing Higher-Order Moments Online
 *   http://people.xiph.org/~tterribe/notes/homs.html
 *
 * Pébay, Philippe (2008), "Formulas for Robust, One-Pass Parallel Computation of Covariances and Arbitrary-Order Statistical Moments", Technical Report SAND2008-6212, Sandia National Laboratories
 *   http://prod.sandia.gov/techlib/access-control.cgi/2008/086212.pdf
 *
 * Last one contains all nacasarry information for calculating Central moments of any degre
 * and Covariance in one-pass parallel or online algorithm in efficient and stable way.
 */
final class CentralMoments4(T) {
private:
	/// number of data points
	uint n;

	/// mean
	T mu = 0.0;

	/// central moments
	// T M_0; == n
	// T M_1; == 0
	T M_2 = 0.0;
	/// ditto
	T M_3 = 0.0;
	/// ditto
	T M_4 = 0.0;

public:

version(opcat_old_version) {
} else {
	/** Add one point of data to this class updating all moments if nacassary.
	 *
	 * If it was already called with elements of multiset {y_1, ..., y_k} (in any order),
	 * then after this method internal statistics (and data returned by other metods)
	 * will describe multiset {y_1, ..., y_k} + y. Values can apear multiple times.
	 *
	 * This operation is O(1) in time and space. Suitable for online computation, with single pass over data, without additional storage.
	 *
	 * Algorithm should be numerically stable. So techniques like Kahan method isn't needed.
	 *
	 * This operation is O(1), but slightly costly (1 fdiv, 6 add/sub, 3-4 mul, 9 fadd/sum, 13-14 fmul).
	 *
	 * It takes, 31-52 microseconds (about 84 instructions) on Pentium M 1.7 GHz.
	 */
	void opCatAssign(T y) {
		alias real TTemp;
		const TTemp u_delta = (y-mu);
		const new_n = n+1;
		const TTemp u_delta_over_n = u_delta/new_n;

		const TTemp u_delta_over_n_2 = u_delta_over_n*u_delta_over_n;
		const TTemp new_mu = mu + u_delta_over_n;

		alias long LL;

		const long ntemp = cast(LL)n*(cast(LL)n-1);
		version(paranoicchecks) {
				const temp2 = cast(LL)new_n*cast(LL)new_n - 3*cast(LL)new_n + 3;
				assert(ntemp+1 == temp2);
		}

		version(paranoicchecks) assert(ntemp+1 == (cast(LL)new_n*cast(LL)new_n - 3*cast(LL)new_n + 3));
		version(paranoicchecks) assert(ntemp+1 == cast(LL)n*(cast(LL)n-1) + 1);
		version(paranoicchecks) assert(cast(LL)(new_n-1)*(cast(LL)new_n*cast(LL)new_n - 3*cast(LL)new_n + 3) == cast(LL)n*(ntemp+1));

		M_4 +=
			- 4.0*M_3*u_delta_over_n
			+ 6.0*M_2*u_delta_over_n_2
			+ u_delta_over_n_2*u_delta_over_n*u_delta*n*(ntemp+1);

		version(paranoicchecks) assert((cast(LL)new_n-1)*(cast(LL)new_n-2) == cast(LL)n*(cast(LL)n-1));
		version(paranoicchecks) assert(cast(LL)n*(cast(LL)n-1) == (cast(LL)new_n*cast(LL)new_n - 3*cast(LL)new_n + 2));
		version(paranoicchecks) assert((cast(LL)new_n*cast(LL)new_n - 3*cast(LL)new_n + 2) == ntemp);

		M_3 +=
			- 3.0*M_2*u_delta_over_n
			+ u_delta_over_n_2*u_delta*ntemp; //(new_n-1)*(new_n-2) == n*(n-1) == new_n^2 - 3*new_n + 2 == (ntmp+1) - 1 == ntmp
		M_2 +=
			+ (y-mu)*(y-new_mu);
		n = new_n;
		mu = new_mu;
	}

	// New version (above) is essential same, but rearanged operations, so few less variables. Also new_n-1 changed to n (because new_n=n+1)
	// and most operations are in-place, because we are calculating moments in reversed order.
	// paranoicchecks, checks if old and new variables are correct. they still need to be checked if are OK for big number of points. n^2 can overflow uint, so maybe use ulong somewhere.

	// Restuls:

/+
n = 100000
u = 3.19926
M_2 = 485447
M_3 = 17546.1
M_4 = 142008
var = 4.85452
skew = 0.0164047
kurtosis = 0.06026
Mean: 3.199256  Var:  4.854517  Std-Dev:  2.203297  N:  100000
+/

}

version(opcat_old_version) {
	// This operation is O(1), but slightly costly (1 fdiv, 6 add/sub, 3-4 mul, 9 fadd/sum, 13-14 fmul).
	// It takes, 49-64 microseconds (about 90 instruction) on Pentium M 1.7 GHz.
	void opCatAssign(T y) {
		const u_delta = (y-mu);
		const new_n = n+1;
		const u_delta_over_n = u_delta/new_n;

		const new_mu = mu + u_delta_over_n;

		const new_M_2 = M_2
			+ (y-mu)*(y-new_mu);
		const u_delta_over_n_2 = u_delta_over_n*u_delta_over_n;
		const new_M_3 = M_3
			- 3.0*M_2*u_delta_over_n
			+ u_delta_over_n_2*u_delta*n*(new_n-2);  // TODO: new_n-1==n, so: n*(n-1): or expand it, so it will be (new_n^2 - 3*new_n + 2), like in new_M_4
		alias long LL;
		long temp2 = cast(LL)new_n*cast(LL)new_n - 3*cast(LL)new_n + 3;
		//assert(new_n*new_n - 3*new_n + 3 == temp2);
		const new_M_4 = M_4
			- 4.0*M_3*u_delta_over_n
			+ 6.0*M_2*u_delta_over_n_2
			+ u_delta_over_n_2*u_delta_over_n*u_delta*n*temp2; //(new_n*new_n - 3*new_n + 3);
		n = new_n;
		mu = new_mu;
		M_2 = new_M_2;
		M_3 = new_M_3;
		M_4 = new_M_4;
	}

/+
n = 100000
u = 3.19926
M_2 = 485447
M_3 = 17546.1
M_4 = 142008
var = 4.85452
skew = 0.0164047
kurtosis = 0.06026
Mean: 3.199256  Var:  4.854517  Std-Dev:  2.203297  N:  100000
+/

}

	/** Reverse operation than opCatAssign, remove point from a sequence.
	 *
	 * This operation is O(1), but not optimalised heavly like opCatAssign.
	 *
	 * Essentially it calculates internal statistics identical to the situation
	 * when one value y was never in original sequence added in opCatAssign.
	 *
	 * Note: Don't call this metho with y which wasn't added using opCatAssign,
	 *       because it will give nonsensical values.
	 *
	 * Note: You should call this method only on values which was added by opCatAssign,
	 *         and only such number of times as it was added by opCatAssign.
	 *
	 * Remark: You can remove points from multiset in any order (not nacasarly reversed opCatAssign order).
	 */
	void opSubAssign(T y) {
		throw new Error("not implmented");
	}

	/** Merge two CentralMoments4 data structures, producing CentralMoments of all datapoints.
	 *
	 * This operation is O(1), but slightly costly (1 fdiv, 4 add/sub, 5 mul, 13 fadd/sum, 21 fmul), so one can first partition data, compute Moments in parallel and then quickly merge them
	 */
	auto opCat(CentralMoments4!(T) B) { // const
		auto A = this;
		auto AB = new CentralMoments4!(T)();
		AB.n = A.n + B.n;
		const u_delta = B.mu-A.mu;
		const u_delta_over_n = u_delta/AB.n; // Ab.n or n
		AB.mu = A.mu + B.n*u_delta_over_n; // CGL79

		alias long LL;

		const AnBn = cast(LL)A.n*cast(LL)B.n;
		AB.M_2 = A.M_2 + B.M_2
			+ AnBn*u_delta_over_n*u_delta; // CGL79
		const u_delta_over_n_2 = u_delta_over_n*u_delta_over_n;
		AB.M_3 = A.M_3 + B.M_3
			+ u_delta_over_n_2*u_delta*AnBn*(A.n-B.n)
			+ 3.0*(A.n*B.M_2 - B.n*A.M_2)*u_delta_over_n; // Ter08
		const An2 = cast(LL)A.n*cast(LL)A.n;
		const Bn2 = cast(LL)B.n*cast(LL)B.n;
		AB.M_4 = A.M_4 + B.M_4
			+ u_delta_over_n_2*u_delta_over_n*u_delta*AnBn*(An2 - AnBn + Bn2)
			+ 6.0*(An2*B.M_2 + Bn2*A.M_2)*u_delta_over_n_2
			+ 4.0*(A.n*B.M_3 - B.n*A.M_3)*u_delta_over_n; // Ter08
		return AB;
	}

	/// Inplace opCat
	void opCatAssign(CentralMoments4!(T) B) {
		throw new Error("not implmented");
	}

	/// Reverse operation than opCat
	auto opSub(CentralMoments4!(T) B) { // const
		return B;
	}

	///
	void consume(Generator!(T) g) {
		foreach (x; g) {
			opCatAssign(x);
		}
	}

	/// Moments
	T zeroth_central_moment() const { return n; }
	///
	T first_central_moment() const { return 0.0; }
	///
	T second_central_moment() const { return M_2; }
	///
	T third_central_moment() const { return M_3; }
	///
	T forth_central_moment() const { return M_4; }

	/// sample moments
	T m0() const { return 1.0; }
	///
	T m1() const { return 0.0; }
	///
	T m2() const { return M_2/n; }
	///
	T m3() const { return M_3/n; }
	///
	T m4() const { return M_4/n; }


	/// mean
	T mean() const { return mu; }
	/// unbiased estimator of variance
	T variance() const { return M_2/(n-1); }
	/// sqrt(var)
	T deviation() const { return sqrt(variance); }

	/// g_1 (sample skewness)
	//T skewness() const { return sqrt(cast(T)n)*M_3 / (sqrt(M_2)*M_2); }
	//T skewness() const { auto k2 = k_2(); return k_3() / (sqrt(k2)*k2); }
	T skewness() const { auto m2_ = m2(); return m3 / (sqrt(m2_)*m2_); }
	/// g_1
	T g_1() const { return skewness(); }
	/// population skewness, G_1 = k_3 / k_2^{3/2}
	T G_1() const { return sqrt(cast(T)(n)*cast(T)(n-1))*g_1()/cast(T)(n-2); }
	/// ditto
	T population_skewness() const { return G_1(); }

	/// Kurtosis. g_2
	//T kurtosis() const { return (cast(T)n)*M_4 / (M_2*M_2); }
	//T kurtosis() const { auto k2 = k_2(); return k_4() / (k2*k2); }
	T kurtosis() const { auto m2_ = m2(); return m4() / (m2_*m2_); }
	/// g_2
	T g_2() const { return kurtosis(); }
	/// population kurtosis, G_2 = k_4 / k_2^2
	T G_2() const { return (n-1)*((n+1)*g_2() + 6) / (n-2) / (n-3); }
	/// ditto
	T population_kurtosis() const { return G_2(); }
	/// Excess kurtosis
	T excess_kurtosis() const { return kurtosis() - 3.0; }
	/// ditto
	T population_excess_kurtosis() const { return G_2() - 3.0; }

	/// Cumulants calculated from momemnts
	T k_1() const { return mu; }
	///
	T k_2() const { return M_2; }
	///
	T k_3() const { return M_3; }
	///
	T k_4() const { auto k2 = k_2(); return M_4 - 3.0*k2*k2; } // M_4 = k_4 + 3k_2^2

	/// number of data points
	uint count() const { return n; }

	///
	void dump() const {
		writefln("n = %d", n);
		writefln("mu = %.8g", mu);
		writefln("M_2 = %.8g", M_2);
		writefln("M_3 = %.8g", M_3);
		writefln("M_4 = %.8g", M_4);
		writefln("var = %.8g", variance());
		writefln("dev = %.8g", deviation());
		writefln("skew (g1) = %.8g", skewness());
		writefln("pop. skew (G1) = %.8g", G_1());
		writefln("kurtosis (g2) = %.8g", kurtosis());
		writefln("excess kurtosis = %.8g", excess_kurtosis());
		writefln("pop. kurtosis (G2) = %.8g", G_2());
	}
}


/+
Benchmarks:

!float and all float:

1000000    52325269    52325269          52     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv //  wersja przyopyumalizowana
 100000     3329632     3329632          33     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv
 100000     4713355     4713355          47     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv
1200000    59042145    59042145          49     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv
1000000    43503207    43503207          43     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv
1000000    39674069    39674069          39     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv
1000000    44485670    44485670          44     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv
 100000     3192523     3192523          31     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv

1000000    26065973    26065973          26     _D4hist9Histogram11opCatAssignMFfZv
 100000     1625359     1625359          16     _D4hist9Histogram11opCatAssignMFfZv
 100000     1283064     1283064          12     _D4hist9Histogram11opCatAssignMFfZv
1200000    24171976    24171976          20     _D4hist9Histogram11opCatAssignMFfZv
1000000    22676347    22676347          22     _D4hist9Histogram11opCatAssignMFfZv

1000000    22265460    22265460          22     _D4hist14__T7MomentsTfZ7Moments11opCatAssignMFfZv
 100000     2014933     2014933          20     _D4hist14__T7MomentsTfZ7Moments11opCatAssignMFfZv
 100000     1714127     1714127          17     _D4hist14__T7MomentsTfZ7Moments11opCatAssignMFfZv
1200000    27080748    27080748          22     _D4hist14__T7MomentsTfZ7Moments11opCatAssignMFfZv
1000000    24887533    24887533          24     _D4hist14__T7MomentsTfZ7Moments11opCatAssignMFfZv


1000000    46942741    46942741          46     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv // 1.7 ghz
1000000    34645162    34645162          34     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv // 1.7 ghz
11000000   561457037   561457037          51     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv // 1.7 ghz

wersja starsza:
1000000    52433586    52433586          52     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv
1000000    62209959    62209959          62     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv
1000000    51220722    51220722          51     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv
 100000     4968509     4968509          49     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv
1000000    64327111    64327111          64     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv
 200000    13845377    13845377          69     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv

1000000    49648119    49648119          49     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv // 1.7 ghz
1000000    53855740    53855740          53     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv // 1.7 ghz


// przypomptymalizowana z real na wazne zmienne chwilowe:
1000000    70040455    70040455          70     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv
1000000    67945573    67945573          67     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv

// double
1000000    36358382    36358382          36     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv
1000000    52299268    52299268          52     _D4hist23__T15CentralMoments4TfZ15CentralMoments411opCatAssignMFfZv


!double:
1000000    56388447    56388447          56     _D4hist23__T15CentralMoments4TdZ15CentralMoments411opCatAssignMFdZv
1000000    86565933    86565933          86     _D4hist23__T15CentralMoments4TdZ15CentralMoments411opCatAssignMFdZv
1000000    62852809    62852809          62     _D4hist23__T15CentralMoments4TdZ15CentralMoments411opCatAssignMFdZv


!double
n = 1000000
u = 5.3399634
M_2 = 22778723
M_3 = 6362915.6
M_4 = 2299537.9
var = 22.778746
dev = 4.7727084
skew = 0.058527787
pop. skew = 0.058527875
kurtosis = 0.004431817
pop. kurtosis = 0.0044378392

!float with few double
n = 1000000
u = 5.3400612
M_2 = 22766598
M_3 = 6361920.5
M_4 = 2298040.2
var = 22.766621
dev = 4.7714381
skew = 0.058565389
pop. skew = 0.058565475
kurtosis = 0.0044336496
pop. kurtosis = 0.0044396715

!float all
n = 1000000
u = 5.3400612
M_2 = 22766598
M_3 = 6361917.5
M_4 = 2298033.8
var = 22.766621
dev = 4.7714381
skew = 0.05856536
pop. skew = 0.058565449
kurtosis = 0.004433637
pop. kurtosis = 0.004439659

!double with few real
n = 1000000
u = 5.3399634
M_2 = 22778723
M_3 = 6362915.6
M_4 = 2299537.9
var = 22.778746
dev = 4.7727084
skew = 0.058527787
pop. skew = 0.058527875
kurtosis = 0.004431817
pop. kurtosis = 0.0044378392

+/

/** Covariance estimation
 *
 * Uses this update formula.
 *
 * ---
 * C = \sum_{u,v} (u - u')(v - v')
 *
 * C = C' + C'' + n'*n'' / n * (u''-u')*(v''-v') // n = n'+n''
 *
 * C = C' + n' / (n'+1) * (u-u')*(v-v') // for n''=1
 * ---
 *
 * TODO: first we note that Covariance class calculates also mean value of u and v, something which is already calculated in Moments
 * 
 */
final class Covariance(T) {
	T C = 0.0;
	T mu = 0.0;
	T mv = 0.0;
	uint n = 0;

	this() {
	}

	/** Add pair (s,t) to set of number for which we calculate covariance.
	 *
	 * This operation is O(1)
	 *
	 * TODO: This calculates also mean value of u, and v (numbers s, t)
	 */
	void update(T s, T t) {
		alias real TTemp;
		const TTemp u_delta = (s-mu);
		const TTemp v_delta = (t-mv);
		const new_n = n+1;

		C +=  u_delta*v_delta * n / new_n; // for n''=1

		const TTemp u_delta_over_n = u_delta/new_n;
		const TTemp u_delta_over_n_2 = u_delta_over_n*u_delta_over_n;
		mu += u_delta_over_n;

		const TTemp v_delta_over_n = v_delta/new_n;
		const TTemp v_delta_over_n_2 = v_delta_over_n*v_delta_over_n;
		mv += v_delta_over_n;

		n = new_n;
	}

	/** Merge two Covariance structures, and returns Covariance structure describing union of underlaing sets
	 *
	 * This operation is O(1)
	 */
	auto opCat(Covariance!(T) B) { // const
		auto A = this;
		auto AB = new Covariance!(T)();
		AB.n = A.n + B.n;

		// new mean of u
		const u_delta = B.mu-A.mu;
		const u_delta_over_n = u_delta/AB.n; // Ab.n or n
		AB.mu = A.mu + B.n*u_delta_over_n; // CGL79

		// new mean of v
		const v_delta = B.mv-A.mv;
		const v_delta_over_n = v_delta/AB.n; // Ab.n or n
		AB.mv = A.mv + B.n*v_delta_over_n; // CGL79

		// new covariance
		AB.C = A.C + B.C + (u_delta * v_delta * A.n * B.n / AB.n);

		return AB;
	}


	/// Unbiased estimator of covariance
	T cov() const {
		return C/(n-1);
	}

	/// Mean value of u
	T mean_u() const {
		return mu;
	}

	/// Mean value of v
	T mean_v() const {
		return mv;
	}

	/// Number of elements added to this structure
	uint count() const {
		return n;
	}

	/// Print all internall data
	void dump() const {
		writefln("n = %d", n);
		writefln("mu = %.8g", mu);
		writefln("mv = %.8g", mv);
		writefln("C = %.8g", C);
		writefln("cov(u,v) = %.8g", cov());
	}
}


//final class CovarianceMatrix!(T) {
//	
//}

/+
	TODO:
		współczynnik korelacji Pearsona (znormalizowana kowariancja)
		korelacje nielinowe (rangowe), lepsze do modeli nieliowych i bardziej odporne na odstające dane
				Współczynnik korelacji rang Spearmana aka Korelacja rang Spearmana (lub: korelacja rangowa Spearmana, rho Spearmana)
					http://pl.wikipedia.org/wiki/Rho_Spearmana
				# tau Kendalla
					http://pl.wikipedia.org/wiki/Tau_Kendalla
				# gamma Kruskala
				# d Sommersa
 +/


import corod : Generator,
	Array,
	FiberGenerator,
	GeneratorEndException,
	Foreach;

/** 'jack-knife' method for calculating general estimator and their approximated errors
 *
 * estimates[i] = estimator(x less Block_i)
 *
 * avg_estimator = avarage(estimates);
 * err_estimator = \sqrt( N-1 / N \sum_i (estimates[i] - avg_estimator)^2 )
 *
 * TODO: composible estimators: given big set A, and B (eventually of same size),
 *       return small number of variables v_A = f(A), v_B = f(B) (preferebly in O(|A|)), from which
 *       estimators can be computed fast (ie. in O(1)). Also make possible to generate
 *       (possibly in O(1)), v_C = compose(v_A, v_B), which is the same as f(A+B).
 *       Such algorithms are called parallel, or if |B|=1, online algorithms.
 *       Many estimators can be callculated in such ways.
 *
 * TODO: for parallel estimators it is sumtimes possible to incrementally not only add single point or set B,
 *       but reversibly remove it back from the sum,
 *       this enables simple scheme: compute (possibly in parallel) v_X = f(A)
 *       for each i, compute esimates[i] = remove(v_X, x_i)
 *
 * TODO: Calculate multiple estimators in run.
 *
 * Note: generator returned by g_dg_ is traversed many times
 *
 * TODO: divide into blocks in other way: exclude elements with index = i (mod n/blocks)
 *       or for each n/blocks continues values, generate random permutation p, and permute to different blocks
 */
class JackKnife(T, T2) {
public:
	///
	this(Generator!(T) delegate() g_dg_, uint blocks_ = 20, uint n_ = 1_000_000)
	in {
		assert(blocks_ <= n_);
		assert(blocks_ > 1);
	}
	body {
		g_dg = g_dg_;
		n = n_;
		blocks = blocks_;
	}

private:
	Generator!(T) delegate() g_dg;
	uint blocks;
	uint n;

protected:
	class Runer : FiberGenerator!(T) {
		private int block_excluded;
		public this(int block_excluded_) { block_excluded = block_excluded_; }
		protected override void iter() {
			auto G = g_dg();
			foreach (block_i; 0 .. blocks) {
				if (block_i != block_excluded) {
					foreach (i; block_i*n/blocks .. (block_i+1)*n/blocks) {
						yield(G.getNext());
					}
				} else {
					// TODO: skip generator by N values
					foreach (i; block_i*n/blocks .. (block_i+1)*n/blocks) {
						G.getNext();
					}
				}
			}
		}
	}

public:
	///
	void run(T2 delegate(Generator!(T), uint) estimator_dg) {
		// M_full = estimator_dg(g_dg());

		M = new Moments!(T2)();
		foreach (block_excluded; 0 .. blocks) {
			auto x = estimator_dg(new Runer(block_excluded), n - n/blocks);
			writef("%.15g ", x);
			M ~= x;
		}
		writeln();
	}

	Moments!(T2) M_full;

	Moments!(T2) M;

	///
	T2 estimator() {
		return M.mean;
	}

	///
	T2 error() {
		return M.sample_deviation(); // blad sredniej
	}
}

///
interface Estimtion(T, T2) {
	T2 estimate(Generator!(T));
	T2 compose(T2, T2);
	T finalize(T2);
}


/** It computes estimator_dg for each block, and then compose them in jack-knife way
 *
 * generator returned by g_dg_ is traversed only once.
 *
 * TODO: Check if mentioned above technique is numerically stable.
 *
 *       If note it can be applied on slightly smaller scale with blocks. (So only 2*O(b) passes is needed not O(b^2).
 *
 * TODO: Calculate multiple estimators in run (preferably be mean of muliple finalizers?).
 *
 * TODO: remove this ugly parameters to interfacing class.
 *
 * Similar works: Fast Estimators of the Jackknife, Journal article by J.S. Buzas; The American Statistician, Vol. 51, 1997
 * http://www.jstor.org/pss/
 */
class JackKnifeComposable(T, T2, T3) {
public:
	///
	this(Generator!(T) delegate() g_dg, uint blocks_ = 20, uint n_ = 1_000_000)
	in {
		assert(blocks_ <= n_);
		assert(blocks_ > 1);
	}
	body {
		n = n_;
		blocks = blocks_;
		block_estimator.length = blocks;

		G = g_dg();
	}

private:
	T2[] block_estimator;
	Generator!(T) G;
	uint blocks;
	uint n;

protected:
	class Runer : FiberGenerator!(T) {
		protected override void iter() {
			foreach (i; 0 .. n/blocks) {
				yield(G.getNext());
			}
		}
	}

public:
	///
	void run(T2 delegate(Generator!(T), uint) estimator_dg) {
		//auto all = new Moments!(T3)();
		foreach (i, ref b; block_estimator) {
			b = estimator_dg(new Runer(), n/blocks); // compute estimator for each block
			//all ~= b;// wrong type, todo
		}
		//M_full = all;// wrong type, todo
	}

	///
	void run2(T2 delegate(T2, T2) compose_dg, T3 delegate(T2) final_dg) {
		M = new Moments!(T3)();

		foreach (block_excluded; 0 .. blocks) { // for each possible excluded block
			T2 temp;
			bool initialized = false;

			foreach (block_i, block; block_estimator) {
				if (block_i != block_excluded) {
					if (initialized) {
						temp = compose_dg(temp, block); // perform composition
					} else {
						temp = block;
						initialized = 1;
					}
					//temp = (initalized ? compose_dg(temp, block) : (initialized = 1, block));
				} else {
					continue;
				}
			}

			auto x = final_dg(temp);

			writef("%.15g ", x);

			M ~= x; // and avarage over exlusion, and calculate error
		}
		writeln();
	}

	/** For some estimator it is possible to remove point from set so restoring previous value.
	 * In this way we can create simple two pass algorithm with full (non-blocked) jack-knife.
	 * first compose everything, then for each point, uncompose it, estimate, and compose it back.
	 */
	void run2(T2 delegate(T2, T2) compose_dg, T2 delegate(T2, T2) uncompose_dg, T3 delegate(T2) final_dg) {
		M = new Moments!(T3)();

		T2 all;
		bool initialized = false;
		foreach (block_i, block; block_estimator) { // this loop generally can be performed in run(...), and even without constructing subGenerators!
			if (initialized) {
				all = compose_dg(all, block); // perform composition // todo: how about in-plcae composition
			} else {
				all = block;
				initialized = 1;
			}
			//all = (initalized ? compose_dg(all, block) : (initialized = 1, block));
		}

		//M_full = all; // wrong type, todo

		foreach (block_i, block; block_estimator) {
			auto temp = uncompose_dg(all, block); // uncompose this data
			auto x = final_dg(temp);

			writef("%.15g ", x);

			M ~= x; // and avarage over exlusion, and calculate error
		}
		writeln();

	}

	Moments!(T3) M_full;

	Moments!(T3) M;

	///
	T3 estimator() {
		return M.mean;
	}

	///
	T3 error() {
		return M.sample_deviation(); // blad sredniej
	}
}

/** Bootstrap estimation
 *
 * Original genrator is iterated once, and data copied to temporary array.
 *
 * Note: If possible n should be power of 2, because of the way how random generator works,
 *       eventually n << uint.max, so '%' operation have good uniformity.
 */
class Bootstrap(T, T2) {
public:
	/// R should be random number generator on [0, n) or [0, uint.max) is possible.
	this(Generator!(T) delegate() g_dg_, Generator!(uint) R_, uint samples_ = 1_000, uint samples_size_ = 100_000, uint n_ = 1_048_576)
	in {
		assert(samples_ >= 10);
		assert(samples_size_ >= 10);
		//assert(samples_size_/100 <= n);
	}
	body {
		auto G = g_dg_();
		R = R_;
		n = n_;
		samples = samples_;
		samples_size = samples_size_;

		temp.length = n;

		M_full = new Moments!(T2)();
		MM = new MinMax!(T2)();

		foreach (i, ref t; temp) {
			t = G.getNext();
			M_full ~= t;
			MM ~= t;
		}
	}

private:
	Generator!(uint) R;
	uint n;
	uint samples;
	uint samples_size;
	T[] temp;

protected:
	class Runer : FiberGenerator!(T) {
		protected override void iter() {
			foreach (i; 0 .. samples_size) {
					auto random_element = temp[R.getNext() % n]; // random sampling
					yield(random_element);
			}
		}
	}

public:
	///
	void run(T2 delegate(Generator!(T), uint) estimator_dg) {
		M = new Moments!(T2)();
		//M_full.mean + 4.0*M_full.deviation
		H = new Histogram!(T2)(MM.min, MM.max, samples/20);
		foreach (sample_i; 0 .. samples) {
			auto x = estimator_dg(new Runer(), samples_size);
			writef("%.15g ", x);
			M ~= x;
			H ~= x;
		}
		writeln();
	}

	Moments!(T2) M_full;
	MinMax!(T2) MM;

	Moments!(T2) M;
	Histogram!(T2) H;

	///
	T2 estimator() {
		return M.mean;
	}

	///
	T2 error() {
		return M.sample_deviation(); // blad sredniej
	}

	///
	void histogram_dump() {
		H.dump();
	}
}

version(jk_main_test) {

//debug = subfunc_calls;

import corod_random : RandomsUniform, RandomsGaussT, Randoms;
import std.conv : to;

/// average
T avg(T)(Generator!(T) h) {
	auto M = new Moments!(T)();
	M.consume(h);
	return M.mean;
}

alias double MT;

void main(string[] args) {
	uint seed = 128317;
//	uint seed = 12831;

	Generator!(T) g(T)() {
		debug(subfunc_calls) writefln("new g");
		return new RandomsGaussT!(T)(new RandomsUniform!(T)(seed), 1.0, 1.0);
	}

	// simple average
	T est(T)(Generator!(T) h, uint n) {
		debug(subfunc_calls) writefln("avg");
		return avg!(T)(h);
	}

	uint blocks = 100;
	uint samples = 10_000;

	if (args.length > 1) {
		samples = to!(uint)(args[1]);
	}
	if (args.length > 2) {
		blocks = to!(uint)(args[2]);
	}

	auto jk = new JackKnife!(MT, MT)(&g!(MT), blocks, samples);

/*
b=20

n=100
avg = 0.132221
err = 0.106734

n=1_000
avg = -0.003090
err = 0.026954

n=10_000
avg = 0.015995
err = 0.010105

n=100_000
avg = -0.001149
err = 0.003241

n=1_000_000
avg = 0.000653
err = 0.001007

b=100
avg = 0.000653
err = 0.001037

seed2:

n=1_000_000
avg = -0.000346
err = 0.001113

*/

	debug writeln("starting run");
	jk.run(&est!(MT));

	writefln("avg = %.15g", jk.estimator);
	writefln("err = %.15g", jk.error);

	writeln();

	Moments!(T) moments_estimate(T)(Generator!(T) h, uint n) {
		debug(subfunc_calls) writeln("estimate");
		auto M = new Moments!(T)();
		M.consume(h);
		return M;
	}
	Moments!(T) moments_compose(T)(Moments!(T) A, Moments!(T) B) {
		debug(subfunc_calls) writeln("compose");
		return (A ~ B);
	}
	Moments!(T) moments_uncompose(T)(Moments!(T) A, Moments!(T) B) {
		debug(subfunc_calls) writeln("compose");
		return (A - B);
	}
	T moments_finalize_mean(T)(Moments!(T) X) {
		debug(subfunc_calls) writeln("finalize");
		return X.mean;
	}


	auto jk2 = new JackKnifeComposable!(MT, Moments!(MT), MT)(&g!(MT), blocks, samples);

	debug writeln("starting run");
	jk2.run(&moments_estimate!(MT));

	debug writeln("starting compose");

	//jk2.run2(&moments_compose!(MT), &moments_finalize_mean!(MT));

	jk2.run2(&moments_compose!(MT), &moments_uncompose!(MT), &moments_finalize_mean!(MT));

	debug writeln("end compose");

	writefln("avg = %.15g", jk2.estimator);
	writefln("err = %.15g", jk2.error);


	uint seed2 = 85149751;
	auto boot_samples = 65536;
	debug writeln("reading data");
	auto boot = new Bootstrap!(MT, MT)(&g!(MT), new Randoms!(uint)(seed2), 1000, boot_samples, boot_samples);

	debug writeln("starting bootstrap");
	boot.run(&est!(MT));

	writefln("avg = %.15g", boot.estimator);
	writefln("err = %.15g", boot.error);

	boot.histogram_dump();

}

}


/** This class computes epsilon-approximated q-quantile for at most
 * n point data set in one pass, with as little of additionall memory as possible.
 *
 * Presented algorithm is deterministic. It always produces correct results.
 *
 * This class uses algorithms which ensures that e-approximated q-quantile,
 * And elements from sequence G is said to be and epsilon-approximate
 * of q-quantile if its rank is beetween ceil((q-e)N) and ceil((q+e)N).
 *
 * Space complexity O(\epsilon^{-1} \log^2 (\epsilon N))
 *
 * Work1: "Approximate Medians and other Quantiles in One Pass and with Limited Memory",
 *       Gurmeet Singh Manku, Sridhar Rajagopalan, Bruce G. Lindsay
 *       SIGMOD 1998, Vol 27, No 2, p 426-35, June 1998
 * Patent: Single Pass Space Efficient System and Method for Generating Approximate Quantiles
 *         Satisfying an Apriori User-Defined Approximation Error by B G Lindsay, G S Manku,
 *         S Rajagopalan, US Patent #06108658, Issued: Aug 22, 2000.
 *
 * There is more general algorithm, which doesn't need a-priory knowledge about size of dataset (n):
 *
 * Work2: "Random Sampling Techniqus for Space Efficient Online Computation of Order Statistics of Large Datasets",
 *       Gurmeet Singh Manku, Sridhar Rajagopalan, Bruce G. Lindsay
 *
 * Patent: Single Pass Space Efficient System and Method for Generating an Approximate Quantile
 *         in a Data Set Having an Unknown Size by B G Lindsay, G S Manku, S Rajagopalan,
 *         US Patent #06343288, Issued: Jan 29, 2002.
 *
 * TODO: "Space-Efficient Online Computation of Quantile Summaries", Michael Greenwald, Sanjeev Khanna
 *       Jest to algorytm ktory ma asymptotycznie mniejsze zuzycie pamieci oraz nie potrzebuje znac rozmiaru danych.
 *       W porownaniu do Work2, jest to algorytm deterministyczny. (W Work2 mamy jedyie 1-delta szansy na to ze mamy,
 *       e-aproskymacje). O(\epsilon^{-1} \log (\epsilon N)).
 *
 * Note: All mentioned here algorihtm can be easly parallelized.
 */
class Quantile(T = float) {
private:
	// b - number of buffers
	uint b;
	// k - number of elements of type T in each buffer
	uint k;

	final class Buffer {
		T[] d; // data in buffer
		uint w; // weight
		bool full_; // full or empty
		int l_; // level;
		bool aux = false;

		// properties
		bool full() { return full_; }
		bool full(bool full2) {
			if (!aux) {
				if (full2 && !full_) { // update counts of full/empty buffers on change
					count_full_++;
					count_empty_--;
				} else if (!full2 && full_) { // ditto
					count_full_--;
					count_empty_++;
					if (smallest_l_ == l_) {
						smallest_l_ = -1; // invalidate cache for safty
					}
				}
			}
			full_ = full2;
			return full2;
		}

		int l() { return l_; }
		int l(int l2) {
			if (!aux) {
				if (full_ && smallest_l_ >= 0) {
					if (l2 < smallest_l_) {
						smallest_l_ = l2; // this is new smallest value
					} else if (l_ == smallest_l_ && l2 > smallest_l_) {
						smallest_l_ = -1; // invalidate cache, and perform full search in smallest_l() for safty
					}
				}
			}
			l_ = l2;
			return l2;
		}
	}

	Buffer[] buffers;

	bool ended;

	uint maxN;
	float epsilon;

public:
	/// prepears nacassary structures for estimation of quantiles with epsilon-approximation
	this(uint maxN_, float epsilon_ = 0.001) {
		ended = false;

		maxN = maxN_;
		epsilon = epsilon_;

		mem(maxN, epsilon, b, k);

		//writefln("b=%d k=%d", b, k);

		buffers.length = b;

		foreach (ref buf; buffers) {
			buf = new Buffer();
			buf.d.length = k;
			buf.w = 0;
			buf.l = 0;
			buf.full = false;
		}

		multi_i.length = b;
		to_collapse.length = b;

		temp = new Buffer();
		temp.aux = true;
		temp.d.length = k;

		count_full_ = 0;
		count_empty_ = b;

		smallest_l_ = -1;

		cbuf.length = k;
	}

	/** Returns estimated number of bytes which will b used by algorithm */
	static size_t mem(uint maxN, float epsilon, ref uint b, ref uint k) {
		static const float[] rows = [
			0.100,
			0.050,
			0.010,
			0.005,
			0.001,
			0.0005, // computed by witold baryluk
			0.0001 // computed by witold baryluk
		];

		static const cols = [100_000, 1_000_000, 10_000_000, 100_000_000, 1_000_000_000 /*, 10_000_000_000*/];

		static const bs = [
			[5, 7, 10, 15, 12, 17],
			[6, 6, 8, 7, 8, 19],
			[7, 12, 9, 10, 10, 12],
			[3, 8, 8, 8, 7, 8],
			[3, 5, 5, 9, 10, 15],
			[3, 3, 8, 8, 8, 7],
			[2, 3, 5, 12, 9, 10]
		];
		static const ks = [
			[55, 54, 60, 51, 77, 69],
			[78, 117, 129, 211, 235, 116],
			[217, 229, 412, 596, 765, 767],
			[953, 583, 875, 1290, 2106, 2341],
			[2778, 3031, 5495, 4114, 5954, 5099],
			[4762, 9524, 5828, 12900, 21052],
			[16667, 27778, 30304, 22894, 41136, 59538]
		];

		uint r;
		for (r = 0; r < rows.length; r++) {
			if (rows[r] <= epsilon) {
				break;
			}
		}
		if (r == rows.length) {
			throw new Exception("too small epsilon");
		}

		uint c;
		for (c = 0; c < cols.length; c++) {
			if (cols[c] >= maxN) {
				break;
			}
		}
		if (c == cols.length) {
			throw new Exception("too big n");
		}

		b = bs[r][c];
		k = ks[r][c];

		return (b+1)*k*T.sizeof + Buffer.sizeof*(b+1);
	}



	invariant() {
		assert(b >= 2);

		assert(count_empty_ + count_full_ == b);

		uint countW = 0;
		foreach (bi, ref buf; buffers) {
			if (buf.full) {
				countW += buf.w;
			}
		}
		assert(W == countW);
	}


private:
	uint N; // counts how many we read from input generator
	uint infs; // counts how many -inf/+inf we additionally added to buffers

	uint W = 0; // sum of weights in full buffers, it is sufficient to increment it by one in NEW,
				// because in COLLAPSE we destroy wY=\sum_w(Xi) of weights and create one buffer with wY

	/// consumes at most k values from G, put them into buf.d and assign weight=1, and level=l
	void NEW(Buffer buf, Generator!(T) G, uint l)
	in {
		assert(count_empty() >= 1);
	}
	body {
		assert(buf.full == false);
		uint i;
		try {
			for (i = 0; i < k; i++) {
				buf.d[i] = G.getNext();
				N++;
			}
		} catch (GeneratorEndException gee) {
			if (gee.who is G) {
				for (; i < k; i++) {
					if (i & 1) {
						buf.d[i] = -1.0/0.0;
					} else {
						buf.d[i] = 1.0/0.0;
						}
					infs++;
				}
				throw gee;
			} else {
				throw gee;
			}
		} finally {
			buf.full = true;
			buf.w = 1;
			buf.l = l;
			buf.d.sort; // so sorting would is not needed in COLLAPSE
			W += 1;
			assert(N <= maxN);
		}
	}

	/// cachine of smallest l
	int smallest_l_;

	/// This retrivies smallest l from cache, or do search if
	/// cache was somehow invalidated.
	/// search is fast, because we have at most 30 buffers.
	public int smallest_l()
	out(ret) {
		assert(ret >= 0);
		assert(smallest_l_ >= 0);
		assert(ret == smallest_l_);
		assert(ret < int.max);
	}
	body {
		if (smallest_l_ >= 0) { // caching
			return smallest_l_;
		}
		int min_l = int.max;
		foreach (bi, ref buf; buffers) {
			if (buf.full && buf.l < min_l) {
				min_l = buf.l;
			}
		}
		//assert(min_l < int.max);
		smallest_l_ = min_l;
		return min_l;
	}

	/// remembers how many empty and full buffers we have (not counting temporary?)
	uint count_empty_;
	/// ditto
	uint count_full_;

	/// return number of empty buffers
	public uint count_empty()
	out(ret) {
		assert(ret == count_empty_);

		assert(ret >= 0);
		assert(ret <= b);
		uint count = 0;
		foreach (bi, ref buf; buffers) {
			if (!buf.full) {
				count++;
			}
		}
		//writefln("%d = %d = %d", count, ret, count_empty_);
		assert(count == ret);
	}
	body {
		return count_empty_;
	}

	/// remembers which one of two possible positions for wY even perform
	bool had_even_wY = false;

	/// temporary buffer for COLLAPSE
	Buffer temp;
	/// multiple indexes for multi-merge algorithm
	uint[] multi_i;

	/// collapse set of buffers pointed by pointers in slice bufs,
	/// collapse it to temp, and set it level l
	void COLLAPSE(Buffer[] bufs, uint l) {
		assert(bufs.length >= 2);
		uint wY = 0;
		foreach (ref buf; bufs) {
			assert(buf.full == true);
			buf.full = false; // after this method all buffers will be empty
			wY += buf.w;
		}

		uint j = 0;

		int jw(int jj) {
			if (wY % 2 == 1) {
				return jj*wY + (wY + 1) / 2;
			} else {
				if (had_even_wY) {
					return jj*wY + (wY + 2) / 2;
				} else {
					return jj*wY + wY / 2;
				}
			}
		}

		uint next = jw(j);

		// multi-merge
		uint counter = 0;
		multi_i[] = 0;
		while (true) {
			uint min_bi = bufs.length;
			T min_v = 1.0/0.0;
			uint taken = false;
			foreach (bi, ref buf; bufs) {
				if (multi_i[bi] < k) { // if our index is still not hitted end of buffer
					if (buf.d[multi_i[bi]] <= min_v) {
						min_bi = bi;
						min_v = buf.d[multi_i[bi]];
						taken = true;
					}
				}
			}
			//writefln("min_bi = %d", min_bi);
			if (taken == false) {
				//writeln("taken false");
				break;
			}
			assert(min_bi < bufs.length);
			//writefln("min_v = %f", min_v);
			/+ if (!(min_v < 1.0/0.0)) {
				writefln("min_bi=%s", min_bi);
				writefln("buf[bi][]=%s", bufs[min_bi].d);
				writefln("multi_i[]=%s", multi_i);
				foreach (bi, ref buf; bufs) {
					if (multi_i[bi] < k) { // if our index is still not hitted end of buffer
						writefln("multi_i[%d]=%d, val=%f", bi, multi_i[bi], buf.d[multi_i[bi]]);
					} else {
						writefln("multi_i[%d]=%d(end), val=none", bi, multi_i[bi]);
					}
				}
				writeln();
			}
			+/
			debug if (!(min_v < 1.0/0.0)) {
				assert(min_v == bufs[min_bi].d[multi_i[min_bi]]);
			} else {
				assert(min_v < 1.0/0.0);
			}
			if (taken == false) {
				break;
			}

			counter += bufs[min_bi].w;
			assert(counter <= wY*k);
			if (counter >= next) {
				assert(j < k);
				temp.d[j] = bufs[min_bi].d[multi_i[min_bi]];
				j += 1;
				next = jw(j);
			}
			multi_i[min_bi]++;
		}
		assert(j == k);
		temp.full = true;
		temp.w = wY;
		temp.l = l;

		foreach (ref buf; bufs) {
			assert(buf.full == false);
			buf.l = 0;
			buf.w = 0;
		}

		auto target = bufs[0];

		// swap data
		auto t = target.d;
		target.d = temp.d;
		temp.d = t;

		temp.full = false;
		target.full = true;
		target.w = wY;
		target.l = l;

		if (wY % 2 == 0) {
			had_even_wY = !had_even_wY;
		}
	}

	T OUTPUT(uint p)
	in {
		assert(p < W*k);
		assert(count_full_ >= 1);
	}
	body {
		auto bufs = buffers;

		// multi-merge
		uint counter = 0;
		multi_i[] = 0;
		while (true) {
			uint min_bi = 0;
			T min_v = 1.0/0.0;
			uint taken = false;
			foreach (bi, ref buf; bufs) {
				if (buf.full) { // only perform search on full buffers
					if (multi_i[bi] < k) { // if our index is still not hitted end of buffer
						if (buf.d[multi_i[bi]] <= min_v) {
							min_bi = bi;
							min_v = buf.d[multi_i[bi]];
							taken = true;
						}
					}
				}
			}
			assert(min_v < 1.0/0.0);
			if (taken == false) {
				break;
			}
			counter += bufs[min_bi].w;
			if (counter >= p) { // stop when needed
				return bufs[min_bi].d[multi_i[min_bi]];
			}
			multi_i[min_bi]++;
		}
		assert(0);
	}


	Buffer[] to_collapse;

public:

	void opCatAssign(Generator!(T) G) {
		while (true) {
			next_faza();
			consume_step(G);
		}
	}

private:
	/// calculate what next faze we should peform
	void next_faza() {
		if (!ended) {
			try {
				while (/*true*/ in_faza == 0) {
					uint how_many_empty = count_empty();
					if (how_many_empty == 1) {
						int l = smallest_l();
						Buffer to_new;
						foreach (b, ref buf; buffers) {
							if (buf.full == false) {
								to_new = buf;
								break;
							}
						}
						assert(to_new !is null);
						in_faza = 1;
						faza_1_l = l;
						faza_1_to_new = to_new;
					} else if (how_many_empty >= 2) {
						in_faza = 2;
						faza_2_buffers = buffers;
						faza_2_c = how_many_empty;
					} else {
						int l = smallest_l();
						assert(how_many_empty == 0);
						uint j = 0;
						foreach (b, ref buf; buffers) {
							if (buf.full && buf.l == l) {
								to_collapse[j++] = buf;
							}
						}
						COLLAPSE(to_collapse[0 .. j], l+1);
					}
				}
			} catch (GeneratorEndException gee) {
				/+// todo: move this to opCatAssign(G)
					if (gee.who is G) {
					; // generator G ended, no problem.
					// all relevant operations (like flags, and auxilary
					// variables was setup in safe order
					return;
				} else {
					throw gee;
				}
				+/
			} finally {
				//ended = true;
			}
		} else {
			throw new Exception("already called get()");
		}
	}

	/// in what phase we are
	int in_faza = 0;
	
	/// if in phase 1, what level new buffer should have, and where it is
	int faza_1_l = 0;
	/// ditto
	Buffer faza_1_to_new;

	/// if in phase 2, how many buffer do we need to read (c), and in what slice we should search for free bufs
	int faza_2_c = 0;
	/// ditto
	Buffer[] faza_2_buffers;

	/// consume some elements from g according to the phase (at most k for 1, or at most c*k for 2)
	/// if g is long enaugh exactly k or c*k will be read from it
	void consume_step(Generator!(T) G) {
		assert(!ended);
		switch (in_faza) {
			case 1:
				assert(faza_1_to_new !is null);
				NEW(faza_1_to_new, G, faza_1_l);
				break;
			case 2:
				assert(faza_2_buffers.length > 0);
				int bc = 0;
				foreach (b, ref buf; faza_2_buffers) {
					if (buf.full == false) {
						NEW(buf, G, 0);
						bc++;
					}
				}
				assert(bc >= 2);
				break;
			default:
				throw new Exception("unfazed");
		}
	}

	/// consume at most k elements from g
	/// if g is more than k elements, exactly k elements will be read.
	void consume_small_step(Generator!(T) G) {
		assert(!ended);
		switch (in_faza) {
			case 1:
				assert(faza_1_to_new !is null);
				NEW(faza_1_to_new, G, faza_1_l);
				in_faza = 0;
				break;
			case 2:
				assert(faza_2_buffers.length > 0);
				assert(faza_2_c > 0);
				foreach (b, ref buf; faza_2_buffers) {
					if (buf.full == false) {
						NEW(buf, G, 0);
						faza_2_buffers = faza_2_buffers[b .. $];
						faza_2_c--;
						break;
					}
				}
				if (faza_2_c == 0) {
					in_faza = 0;
				}
				break;
			default:
				throw new Exception("unfazed");
		}
	}

	/// buffer for opCatAssign(T)
	T[] cbuf;
	/// index in this buffer
	int cbuf_i = 0;

public:
	void begin() {
		next_faza();
	}

	/// add value x
	void opCatAssign(T x) {
		assert(!ended);
		assert(cbuf !is null);
		assert(cbuf.length == k);
		assert(cbuf_i >= 0 && cbuf_i < k);
		cbuf[cbuf_i++] = x;
		if (cbuf_i >= k) {
			consume_small_step(new Array!(T)(cbuf));
			cbuf_i = 0;
			next_faza();
		}
	}

	/// finalize all data
	void end() {
		if (cbuf_i > 0) {
			auto tempG = new Array!(T)(cbuf[0 .. cbuf_i]);
			try {
				consume_small_step(tempG);
			} catch (GeneratorEndException gee) {
				if (gee.who is tempG) {
					;
				} else {
					throw gee;
				}
			}
		}
		ended = true;
	}

	/// This function should be only called once, after all points was added by opIndex
	/// it is safe to call this function multiple times with different values q
	T quantile(float q) {
		assert(epsilon <= q && q <= 1.0-epsilon);
		ended = true;
		assert(N > 0);
		assert(W > 0);
		auto beta = cast(float)(N + infs) / N;
		auto q2 = (2*q + beta - 1.0) / (2.0*beta);
		//writefln("q=%f  beta=%f q2=%f  W=%d", q, beta, q2, W);
		return OUTPUT(cast(uint)(q2*k*W));
	}
}


/** Moving avarage of order l
 *
 * First l point aren't exactly moving avarage
 */
class MovingAvarage(T) : Generator!(T) {
private:
	Generator!(T) g;
	T[] last;
public:
	this(Generator!(T) g_, int l) {
		last.length = l;
	}
protected:
	override void iter() {
		foreach (ref x; last) {
			x = g.getNext();
		}
		
	}
}

version(q_main_test) {

//debug = subfunc_calls;

import corod_random : RandomsUniform, RandomsGaussT, Randoms;
import std.conv : to;

alias double MT;

void main(string[] args) {
	uint seed = 128317;
//	uint seed = 12831;

	auto g = new RandomsGaussT!(MT)(new RandomsUniform!(MT)(seed), 0.0, 1.0);

	auto maxN = 10_000_000; // points
	auto N = 20; // how many quantiles to output

	auto eps = 0.0001;

	auto q = new Quantile!(MT)(maxN, eps);

//	writefln("b=%d  k=%d", q.b, q.k);

	q.opCatAssign(new Foreach!(MT, g, maxN)());

	multiq(q, N, eps);
}

}


/// computes N equidistant quantiles, and prints equidepth histogram
void multiq(T)(Quantile!(T) q, int N, T eps = 0.0) {
	T[] hh = new T[N];

	if (eps == 0.0) {
		eps = q.epsilon;
	}

	foreach (i; 1 .. N) {
		auto qi = i*(1.0/N);
		auto qqi = q.quantile(qi);
		hh[i] = qqi;
		//writefln("%.4f %.8f", qi, q.quantile(qi));
	}
	// gnuplot> plot "testowanie_histogramu_2.txt" u 2:1 w p, 0.5*erf(x/sqrt(2.0)) + 0.5

	foreach (i; 1 .. N-1) {
		auto qi = i*(1.0/N);
		writefln("%d %.4f %.8f %.8f %.8f %.8f %.8f",
			i, qi, hh[i], hh[i+1],
			(1.0/(hh[i+1]-hh[i])) / N,
			(1.0/(hh[i+1]-hh[i] + 2.0*eps)) / N,
			(1.0/(hh[i+1]-hh[i] - 2.0*eps)) / N);
	}
	// plot "z12c.q" u (($3+$4)/2):7 w boxes, exp(-abs(x))/2 w l

	// gnuplot>
	// plot "testowanie_histogramu_4.hist" u 4:7 w boxes, "testowanie_histogramu_4.hist" u 4:5 w boxes, "testowanie_histogramu_4.hist" u 4:6 w boxes, exp(-x*x / 2) / sqrt(2*pi)
	// plot "testowanie_histogramu_4.hist" u (0.5*($3+$4)):7 w l, "testowanie_histogramu_4.hist" u (0.5*($3+$4)):5 w l, "testowanie_histogramu_4.hist" u (0.5*($3+$4)):6 w l, exp(-x*x / 2) / sqrt(2*pi)
}
