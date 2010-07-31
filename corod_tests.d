module corod_tests;

import corod;
import corod_random;

import hist;

import std.stdio;

import std.conv;

void main(string[] args) {

	ulong l = 1000;
	if (args.length >= 1) {
		l = to!(int)(args[2]);
	}

	switch (to!(int)(args[1])) {
		case 1:
			simpletests();
			break;
		case 2:
			arrays();
			break;
		case 3:
			rs();
			break;
		case 4:
			t();
			break;
		case 5:
			tg(l);
			break;
		case 6:
			tg_covs(l);
			break;
		case 7:
			bufki1(l);
			break;
		case 8:
			int buf_size = to!(int)(args[3]);
			bufki2(l, buf_size);
			break;
	}
}

void simpletests() {
	//auto f = new Test();
	auto f = new Array!(int)([2,6,1,7,2,4,8,2,3,7,2,3,5,8,6,2,2,4]);

	foreach (x; new Foreach!(int, f)()) {
		writefln("foreach test in main: x=%d", x);
	}

	f.reset();
	foreach (int x; new Foreach!(int, f)) {
		writefln("foreach in main: x=%d", x);
	}
}

void arrays() {
	auto a1 = new Array!(int)([4,1,2,5,7,1,2]);
	auto a2 = new Array!(float)([1.2, 3.4, 5.3, 5.1, 2.1, 3.2]);

	concurant_foreach!(int, float)(a1, a2, (int x, float y) {
		writefln("concurant_foreach %d %f", x, y);
	});

	auto a3 = new Array!(int)([4,1,2,5,7,1,2]);
	auto a4 = new Array!(float)([1.2, 3.4, 5.3, 5.1, 2.1, 3.2]);

	auto F = new ConcurantGenerator!(int, float, float)(a3, a4, (int x, float y) {
		return x*y;
	});

	foreach (xy; new Foreach!(float, F)()) {
		writefln("foreach xy: ", xy);
	}

	auto r3 = new Range!(int)(10);

	foreach (x; new Foreach!(int, r3)()) {
		writefln("range: ", x);
	}

	a1.reset();
	writefln("sum: %d", sumint(a1));

	a2.reset();
	writefln("sum: %f", sumfloat(a2));
}

void rs(int l = 1000) {
	auto rr = new Randoms!()();
	foreach (x; new Foreach!(uint, rr, l)()) {
		writeln(x);
	}

	t();
}

void t(int l = 1_000_000) {
	auto rr = new Randoms!(uint)();
	auto ru = new RandomsUniform!(float)(rr, 0.0, 100.0);

	auto h = new Histogram!(float)(0.0, 100.0, 100);
	auto m = new Moments!(float)();
	foreach (x; new Foreach!(float, ru, l)()) {
		//writeln(x);
		h ~= x;
		m ~= x;
	}

	h.dump();
	writefln("Mean: %f  Var:  %f  Std-Dev:  %f  N:  %d", m.mean, m.variance, m.deviation, m.count);
}

void tg(ulong l) {
	auto rr = new Randoms!(uint)(123);
	auto ru = new RandomsUniform!(float)(rr);
	//auto ru = new RandomsUniform!(float)(rr, 0.0, 4.0);
	//auto rg = new RandomsStandardGauss!(float)(ru);
	auto rg = new RandomsGauss(ru, 3.2, 2.2);
	//auto rg = new RandomsGaussFast!(float)(ru, 3.2, 2.2);
	//auto rg = ru;


	//auto rg = new RandomsInverseCDF!(float)(ru, (float x) { return x*x; });

	auto h = new Histogram!(float)(-0.0, 20.0, 100);
	auto m = new Moments!(float)();
	auto cm = new CentralMoments4!(double)();
	auto c = new Covariance!(float)();
	foreach (x; new Foreach!(float, rg, l)()) {
		//writeln(x);
		h ~= x;
		m ~= x;
		cm ~= x;
		c.update(x, x);
	}

	//h.dump();

	cm.dump();
	c.dump();

	//writefln("Mean: %f  Var:  %f  Std-Dev:  %f  N:  %d", m.mean, m.variance, m.deviation, m.count);
}

void tg_covs(ulong l) {
	auto a_rg = new RandomsGauss(8216512, 0.0, 1.0);
	auto b_rg = new RandomsGauss(5772813, 0.0, 1.0);

	//auto a_rg = new RandomsUniform(8216512, 0.0, 1.0);
	//auto b_rg = new RandomsUniform(5772813, 0.0, 1.0);

	auto c = new Covariance!(float)();

	/+concurant_foreach!(float, float)(a_rg, b_rg, (float x, float y) {
		writeln(x, ' ', y);
		c.update(x, y);
	});
	+/

	auto a_m = new Moments!(float)();
	auto b_m = new Moments!(float)();


	auto A_11 = 2.0;
	auto A_21 = 2.0;
	auto A_22 = 5.0;

	// cholesky 2x2 decomposition
	auto L_11 = sqrt(A_11);
	auto L_21 = A_21/L_11;
	auto L_22 = sqrt(A_22 - L_21*L_21);

	auto m1 = 1.0;
	auto m2 = 0.5;

int i = 0;
	foreach (float x1, float y1; new ConcurantForeach!(float, float)(a_rg, b_rg)) {
		// L * (x,y), matrix*vector
		auto x = L_11*x1 + m1;
		auto y = L_21*x1 + L_22*y1 + m2;
		a_m ~= x;
		b_m ~= y;
		c.update(x, y);
		writeln(x, ' ', y);

if (i++ > l) {
	break;
}
	};

	c.dump();
	writeln();

	writefln("corr(x,y) = %.8g", c.cov/sqrt(a_m.variance*b_m.variance));

	writeln();
	writefln("var(x)^2 = %.8f", a_m.variance);
	writefln("var(y)^2 = %.8f", b_m.variance);
	writefln("cov(x,y) = %.8g", c.cov);
}


void bufki1(ulong l) {
	auto rg1 = new Fib1();

	auto aa = 1;

	foreach (x; new Foreach!(int, rg1, l)()) {
		//writefln("rg1 : %d", x);
		aa += x;
	}
	writefln("rg1 : %d", aa);
}

void bufki2(ulong l, int buf_size = 32) {
	auto rg2 = new Fib1Buf();
	rg2.setBuf(buf_size);

	auto aa = 1;

	foreach (x; new Foreach!(int, rg2, l)()) {
		//writefln("rg2 : %d", x);
		aa += x;
	}
	writefln("rg2 : %d", aa);
}

