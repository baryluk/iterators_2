module series;


import std.stdio : writefln, writeln;
import std.cstream : derr;

import corod;
import corod_random;
import autocorel : Autocorellation;
import hist : Histogram;

import core.memory : GC;
import core.thread : Thread; // sleep

/** Differencing of first degree:
 *
 * y'_n = x_n - x_{n-1}
 */
class Differencing : FiberGenerator!(float) {
private:
	/*final*/ Generator!(float) g;
public:
	///
	this(Generator!(float) g_) {
		g = g_;
	}
protected:
	///
	void iter() {
		auto prev = g.getNext();
		while (true) {
			auto current = g.getNext();
			yield(current - prev);
			prev = current;
		}
	}
}

/** Differencing of second degree:
 *
 * y''_n = y'_n - y'_{n-1}
 *
 * where y'_n = x_n - x_{n-1}
 *
 * so
 *
 * y''_n = x_n - x_{n-1} - x_{n-1} + x_{n-2}
 *
 * y''_n = x_n - 2*x_{n-1} + x_{n-2}
 *
 *
 * This is because this process is can be described by lag operators:
 *
 * y''_n = (1 - L)^2 x_n
 *
 * similary for higher orders,
 *
 * y^{(d)}_n = (1 - L)^d x_n = \sum_{i=0}^d ({i \choice d}) x_{n-i}
 *
 */
class DifferencingN(uint N) : FiberGenerator!(float) {
	static assert(0);
}


/** Autoregressive model
 *
 * ---
 * y_n = c + \alpha_0 \eta_n + \sum_{i=1}^p \beta_i y_{n-i}
 * ---
 *
 * Analitical autocovariance is given for AR(1) by:
 * ---
 * \sigma^2 / (1 - \phi^2) * \phi^{|n|}
 * ---
 */
class AR(uint p) : FiberGenerator!(float) {
private:
	/*final*/ Generator!(float) g;
	/*const*/ float c;
	/*const*/ float alpha_0;
	/*const*/ float[p] beta;
	float[p] ys;

public:
	///
	this(Generator!(float) g_, float c_, float alpha_0_, float[/*p*/] beta_, float y0) {
		c = c_;
		g = g_;
		alpha_0 = alpha_0_;
		//beta = beta_;
		assert(beta_.length == p);
		beta[] = beta_[];
		ys[] = y0;
	}

protected:
	/// implementation
	override void iter() {
		while (true) {
			float eta_n = g.getNext();
			float y_n = c + alpha_0 * eta_n;
			foreach (i, b_i; beta) {
				y_n += b_i * ys[p-i-1]; // or [i] or [p-i-1]
			}
			ys[0 .. $-1] = ys[1 .. $]; /// TODO: circular buffer for large p
			ys[p-1] = y_n;
			yield(y_n);
		}
	}
}
/+
Fitting AR(p) model to data:

0. Select order of the model, p.

1. Compute Autocorellation function gamma_m

2. Solve Yule-Walker equation for \phi_k and \sigma_\epsilon^2:
  \gamma_m = \sum_{k=1}^p \phi_k \gamma_{m-k} + \sigma_\epsilon^2 \delta_m

where
  \delta_m is Kronecker's delta
  \sigma_\epsilon^2 is variance of underling distribution

How?

2a. First for m > 1, solve Linear system (symetric?)
2b. then compute \sigma_\epsilon^2 = \phi_0 - \sum_{k=1}^p \phi_k \gamma_{-k}

+/

T[] yule_walker(T)(uint p, Generator!(T) gamma) {
	if (p == 0) {
		return new T[0];
	}

	T c0 = gamma.getNext();
	T inv_c0 = 1.0/c0;

	// direct inversions using least-squares method
	if (p == 1) {
		T r1 = inv_c0*gamma.getNext();
		return [r1];
	}
	if (p == 2) {
		T r1 = inv_c0*gamma.getNext();
		T r2 = inv_c0*gamma.getNext();
		auto temp = 1.0/(1.0 - r1*r1);
		return [r1*(1.0-r2)*temp, (r2-r1 *r1)*temp];
	}

/*
	vector b;
	// symetric
	// Toeplitz-style matrix
	// use Levinson-Durbin algorithm for fast solution of such system
	matrix R;
	return solve_cholesky(R, b);
*/

	return null;
}

/+
T[] burg(uint p, Generator!(T) gamma) {
}


/// Akaike's Information Criterion
uint aic() {
}

/// Rissanen/Schwartz Minimum Description Length (MDL)
/// similar to AIC
+/

/** Brownian motion is particular instance of AR(1) model
 */
class Brownian : AR!(1) {
	this(Generator!(float) g_, float c_, float alpha_0_, float beta_0_, float y0) {
		super(g_, c_, alpha_0_, [beta_0_], y0);
	}
}

/** Moving average
 *
 * ---
 * y_n = \mu + \theta_0 \eta_n + \sum_{i=1}^q \theta_i eta_{n-i}
 * ---
 */
class MA(uint q) : FiberGenerator!(float) {
private:
	/*final*/ Generator!(float) g;
	/*const*/ float mu;
	float[q+1] etas;
	/*const*/ float[q+1] theta;

public:
	///
	this(Generator!(float) g_, float mu_, float[/*q+1*/] theta_) {
		g = g_;
		mu = mu_;
		//theta = theta_;
		assert(theta_.length == q+1);
		theta[] = theta_[];
		etas[] = 0.0f;
	}

protected:
	/// implementation
	override void iter() {
		while (true) {
			//etas[0 .. $-1] = etas[1 .. $]; ///  object.Exception: overlapping array copy
			foreach (i; 0 .. q) {
				etas[i] = etas[i+1]; /// TODO: circular buffer for large p
			}
			etas[q] = g.getNext();
			float y_n = mu;
			foreach (i, theta_i; theta) {
				y_n += theta_i * etas[q-i];
			}
			yield(y_n);
		}
	}
}

/** ARMA = AR+MA ( aka Box-Jenkins )
 *
 * y_t = \sum_{i=1}^p \rho_i y_{t-i} + \theta_0 \eta_t + \sum_{i=1}^q \theta_i \eta_{t-i}
 *
 */
class ARMA(uint p, uint q) : FiberGenerator!(float) {
private:
	/*final*/ Generator!(float) g;
//	/*const*/ float mu;
	float[q+1] etas;
	float[q+1] theta;

	/*const*/ float alpha_0;
	/*const*/ float[p] beta;
	/*const*/ float[p] ys;

public:
	///
	this(Generator!(float) g_, float[/*p*/] beta_, float[/*q+1*/] theta_, float y0) {
		g = g_;
//		mu = mu_;
		theta = theta_;
		beta = beta_;
		etas[] = 0.0f;
		ys[] = y0;
	}

protected:
	/// implementation
	override void iter() {
		while (true) {
			//etas[0 .. $-1] = etas[1 .. $]; ///  object.Exception: overlapping array copy
			foreach (i; 0 .. q) {
				etas[i] = etas[i+1]; /// TODO: circular buffer for large p
			}
			etas[q] = g.getNext();
			float y_n = 0.0;
			foreach (i, theta_i; theta) {
				y_n += theta_i * etas[q-i];
			}
			foreach (i, b_i; beta) {
				y_n += b_i * ys[p-i-1];
			}
			//ys[0 .. $-1] = ys[1 .. $]; /// TODO: circular buffer for large p
			// overlaping array copy :/
			foreach (i; 0 .. p-1) { // this loop can be incorporated directly in the previous one
				ys[i] = ys[i+1];
			}
			ys[p-1] = y_n;

			yield(y_n);
		}
	}
}

/** ARMA with exogenous inputs
 *
 * y_t = \eta_t + \sum_{i=1}^p \rho_i y_{t-i} + \sum_{i=1}^q \theta_i \eta_{t-i} + \sum_{i=1}^b \zeta_i d_{t-i}
 *
 * where $d_i$ is known and external (independed) time serie.
 *
 */
class ARMAX(uint p, uint q, uint b) : FiberGenerator!(float) {
}

/** ARIMA(p,d,q) is particular case of ARMA(p+d,q)
 */
class ARIMA(uint p, uint d, uint q) : FiberGenerator!(float) {
	/*final*/ ARMA!(p,q) arma;
public:
	this(Generator!(float) g_, float[/*p*/] beta_, float[/*q+1*/] theta_, float y0) {
		arma = new ARMA!(p,q)(g_, beta_, theta_, y0);
	}
protected:
	override void iter() {
		static if (d > 0) {
			float prev = arma.getNext();
		}
		while (true) {
			static if (d == 0) {
				yield(arma.getNext());
			} else static if (d == 1) {
				float cur = arma.getNext();
				prev += cur;
				yield(prev);
			} else {
				static assert(0, "not implemented yet (but it is simple");
			}
		}
	}
}

/* ARCH(p)
 *
 * y_i = \sqrt{h_t} \eta_i
 *
 * h_t = \alpha_0 + \sum_{i=1}^p \alpha_i h_{t-i}
 *
 * Where alpha_i > 0
 *
 * h_t == \sigma^2_t
 */
class ARCH(uint p) : FiberGenerator!(float) {
}

/** GARCH(p,q)
 *
 * x_i = \sqrt{\sigma^2_i} \eta_i
 *
 * \sigma^2_t = \sum_{i=1}^p b_i \sigma^2_{t-i} + c_0 + \sum_{i=1}^q c_i x_{t-i}^2
 *
 * Where c_i > 0, b_i > 0
 *
 * Which is equivalent to:
 *
 * \sigma^2_t = \sum_{i=1}^p b_i \sigma^2_{t-i} + c_0 + \sum_{i=1}^q c_i \sigma^2_{t-i} \eta_{t-i}^2
 *
 *
 */
class GARCH(uint p, uint q) : FiberGenerator!(float) {
public:
	this(Generator!(float) g_, float[/*p*/] b_, float[/*q+1*/] c_, float s0 = 0.0) {
		b = b_;
		c0 = c_[0];
		c = c_[1 .. $];
		sigma2s[] = s0;
		x2s[] = 0.0;
		g = g_;
	}
private:
	/*final*/ Generator!(float) g;
	float c0;
	/*const*/ float[q] c;
	/*const*/ float[p] b;
	float[q] x2s;
	float[p] sigma2s;

//	float[p] etas;

protected:

	override void iter() {
		double sigma2_prev = 1.0/0.6;
		double x2_prev = 0.0;
		while (true) {
			double sigma2_i = 1.0 + 0.2 * sigma2_prev + 0.2*x2_prev;
			double eta_i = g.getNext();
			double x_i = sqrt(sigma2_i) * eta_i;
			yield(x_i);
			sigma2_prev = sigma2_i;
			//x2_prev = sigma2_i*eta_i*eta_i;
			x2_prev = x_i*x_i;
		}
	}

/+
	override void iter() {
		while (true) {
			float sigma2_i = c0;
			foreach (i, b_i; b) {
				sigma2_i += b_i * sigma2s[q-1-i];
			}
			foreach (i, c_i; c) {
				sigma2_i += c_i * x2s[p-1-i];
			}

/*
 * x_i = \sqrt{\sigma^2_i} \eta_i
 *
 * \sigma^2_t = \sum_{i=1}^p b_i \sigma^2_{t-i} + c_0 + \sum_{i=1}^q c_i x_{t-i}^2
 */
			assert(sigma2_i >= 0);

			float eta_i = g.getNext();
			float x_i = sqrt(sigma2_i) * eta_i;

			yield(x_i);

			foreach (i; 0 .. q-1) {
				x2s[i] = x2s[i+1];
			}
			x2s[q-1] = x_i*x_i;
			foreach (i; 0 .. p-1) {
				sigma2s[i] = sigma2s[i+1];
			}
			sigma2s[p-1] = sigma2_i;
/*
			foreach (i; 0 .. p-1) {
				etas[i] = etas[i+1];
			}
			etas[p-1] = eta_i;
*/
		}
	}
+/
}

/** Logarithmic interest rates.
 *
 * This is multiplicative version of a Differencing
 *
 * Note: This is only usefull for strictly positive time serie g.
 *
 * Note: First observation can be lost.
 */
class LIR(T) : FiberGenerator!(T) {
private:
	/*final*/ Generator!(T) g;
public:
	///
	this(Generator!(T) g_) {
		g = g_;
	}
protected:
	///
	void iter() {
		auto prev = g.getNext();
		assert(prev > 0.0);
		foreach (current; g) {
			assert(current > 0.0);
			yield(log(current/prev));
			prev = current;
		}
	}
}

/** inverse Logarithmic interest rates.
 *
 * This is multiplicative version of a Integration,
 * it just returns prices given logarithmic interest rates.
 */
class InvLIR(T) : FiberGenerator!(T) {
private:
	/*final*/ Generator!(T) g;
	/*const*/ T y0;
	/*const*/ T d;
public:
	///
	this(Generator!(T) g_, T y0_ = 1.0, T divfactor = 1.0) {
		g = g_;
		y0 = y0_;
		d = 1.0/divfactor;
	}
protected:
	///
	void iter() {
		auto y = y0;
		while (true) {
			y *= exp(g.getNext()*d);
			yield(y);
		}
	}
}


/** This is approximated and fast exponential smoothing.
 *
 * It essentially computed online O(1), sum:
 *
 * \alpha := 1/c
 *
 * s_n := \sum_{i=0}^n y_{n-i} \alpha^i / norm_n
 *
 * where norm_n := \sum_{i=0}^n \alpha^i
 */
class ExponentialSmoothing(T) : FiberGenerator!(T) {
private:
	/*final*/ Generator!(T) g;
	/*const*/ T c;
public:
	///
	this(Generator!(T) g_, T c_ = 2.0) {
		assert(c_ > 0.0);
		g = g_;
		c = 1.0/c_;
	}
protected:
	///
	void iter() {
		auto y = 0.0;
		auto norm = 0.0;
		while (true) {
			y = (y*c + g.getNext());
			norm = (norm*c + 1.0);
			yield(y/norm);
		}
	}
}

/+
Gladzenie:

MA  (sredni ruchome, rozne)
Exponential moving avarage (latwe do obliczenia)
Price oscilator (roznica pomiedzy srednia ruchoma krotka i dluga)
MACD (MACD 12/62)  (srednie kroaczace)  (MACD is a moving avg. of difference beetwen long term and short term moving avg.)
Wstegi Bollingera (trzy linie: n-okresowa srednia ruchoma ceny, +- k odchylen standardowych (a dokladniej srednich ruchomych z odchylen standardowych))
DEMA (ciekawe hpodwojnie eksponencjalna srednia ruchoma, todo)


// floats aren't the best for financial data, but for simulations they are good
struct Price {
	int timestamp;
	uint volumen;
	float open;
	float minimum;
	float maximum;
	float close;
}

wskazniki:
Oscylator stochastyczny; %K = (close - minimum) / (maximum - minimum); %D = srednia 3 dniowa z %K (tzw. szybki oscyl. stoch.)
Wskaznik sily wzglednej: RSI = 1 - 1 / (1 + RS);  RS = (srednia wartosc wzrostu cen zamkniecia z y dni / sr. wart. spadku cen zamkn. z y dni)
    i mniejsze "y" tym czulszy oscylator
+/

/+
// dokladniej zobaczyc w Box and Jenkins (1976)
// oraz w Granger and Newbold (1977)
// obie ksiazki maja nowsze wydania

+/

