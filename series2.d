module series2;

import gnuplot : Gnuplot;

import corod;
import corod_random;
import series;

import hist : MinMax, Histogram, Moments, CentralMoments4, JackKnife, JackKnifeComposable;

import autocorel : Autocorellation;
//import std.stdio : writefln;

void go() {
	auto g = new Gnuplot();
	g.set("terminal x11");
}


enum ProcessType {
	MA,
	AR,
	ARMA,
	ARIMA,
	GARCH,
	Brown
}

import std.conv : to;

void gen(Gnuplot g,
	size_t length = 120,
	uint seed = 1231,
	bool shownoise = true,
	bool showpoints = true,
	ProcessType pt = ProcessType.MA,
	bool showratio = false,
	bool showgrid = false,
	bool separateplots = false,
	bool showratio_logscale = false,
	size_t forecast = 0,
	bool show_invlir = false,
	float invlir_factor = 1.0f,
	bool show_expsmooth = false,
	float expsmooth_factor = 2.0f,
	bool show_autocorel = false,
	size_t autocorel_length = 100,
	float[] ps = null, float[] qs = null,
	bool plots = true,
	bool noise_stats = true,
	float[] save_buffer = null) {


	const bool acki = true; //false;
	const bool teski = false;


	if (plots) {
		g.set("style line 1 lw 2");
	}

	//size_t M = 100000; // number of monte carlo experiments
	size_t M = 1;

	//auto h_end = new Histogram!(float)(-60.0, 60.0, 100); // for Brownian, l=200
	auto h_end = new Histogram!(float)(-4.0, 16.0, 100); // for AR1 a, AR1 b, l = 200

	Generator!(float) rg_org = new RandomsGauss(seed, 0.0f, 1.0f);


	CentralMoments!(T) moments_estimate(T)(Generator!(T) h, size_t n) {
		auto M = new Moments!(T)();
		M.consume(h);
		return M;
	}
	CentralMoments!(T) moments_compose(T)(CentralMoments!(T) A, CentralMoments!(T) B) {
		return (A ~ B);
	}
	CentralMoments!(T) moments_uncompose(T)(CentralMoments!(T) A, CentralMoments!(T) B) {
		return (A - B);
	}
	T moments_finalize_mean(T)(CentralMoments!(T) X) {
		return X.mean;
	}

/+
	alias float MT;
	auto jk2 = new JackKnifeComposable!(MT, Moments!(MT), MT)(&g!(MT), 20, length);
	jk2.run(&moments_estimate!(MT));

	jk2.run2(&moments_compose!(MT), &moments_uncompose!(MT), &moments_finalize_mean!(MT));
+/

	auto rg_minmax = new MinMax!(float)();
	auto rg_moments = new CentralMoments4!(float)();
	auto rg_hist = new Histogram!(float)(-10.0, 10.0, 40);
	if (noise_stats) {
		rg_org = new GeneratorProxy!(float)(rg_org, delegate void(float x) {
				rg_minmax ~= x;
				rg_hist ~= x;
				rg_moments ~= x;
			});
	}

	if (forecast > 0) {
		rg_org = new Concat!(float)([
			new Limit!(float)(rg_org, length),
			new Limit!(float)(new Constant!(float)(0.0f), forecast)
		]);
	}

	Generator!(float) rg = rg_org;

	GeneratorProxy!(float) rg_proxy;
	if (shownoise || showratio) {
		rg_proxy = new GeneratorProxy!(float)(rg_org);
		rg = rg_proxy;
	}

	auto percent = 0;

	//Generator!(float) ar;

	foreach (Mi; 0 .. M) {
		auto p = 100*Mi/M;
		if (p > percent) {
			percent = p;
			derr.writefln("%02d%%", percent);
			if (percent % 10 == 0) {
				GC.collect();
				derr.writefln("GC.collect");
			}
		}

		bool ratio = false;

		Generator!(float) genwraped_ar = new GenWrap!(float, FiberGenerator!(float))(
				delegate FiberGenerator!(float) () { switch (pt) {


				case ProcessType.Brown: return new Brownian(rg, 0.0, 1.0, 1.0, 0.0);

				// examples from lecture06 and lecture07 of PFG's TS lecture
				//case ProcessType.AR: return new AR!(1)(rg, 0.0, 1.0, [/*beta0=*/0.75f], 0.0);
				//new AR!(1)(rg, 0.0, 1.0, [/*beta0=*/-0.75f], 0.0);
				//new AR!(2)(rg, 0.0, 1.0, [/*beta0=*/0.75f, -0.05], 0.0);
				//new AR!(2)(rg, 0.0, 1.0, [/*beta0=*/0.75f, -0.5], 0.0);
				//new AR!(3)(rg, 0.0, 1.0, [/*beta0=*/0.5f, -0.125, 0.5], 0.0);
				//case ProcessType.MA: return new MA!(2)(rg, 0.0, [0.25f, 0.5, 0.25]);
				//new MA!(2)(rg, 0.0, [-0.25f, 0.5, -0.25]);
				//new ARMA!(1,2)(rg, [0.75f], [0.25f, 0.5, 0.25], 0.0f);
				//case ProcessType.ARMA: return new ARMA!(3,0)(rg, [1.0f/6.0f, 2.0f/3.0f, 1.0f/6.0f], [0.25f], 0.0f); // this ARMA is VERY similar to Brown
				//new ARMA!(3,0)(rg, [1.0f/6.0f, 2.0f/3.0f, 1.0f/6.0f], [0.25f], 0.0f);
				//case ProcessType.ARIMA: return new ARIMA!(3,1,0)(rg, [1.0f/6.0f, 2.0f/3.0f, 1.0f/6.0f], [0.25f], 0.0f);


				case ProcessType.AR: return new AR!(1)(rg, 0.0, 1.0, qs[0 .. 1], 0.0);
				case ProcessType.MA: return new MA!(2)(rg, 0.0, ps[0 .. 3]);
				case ProcessType.ARMA: return new ARMA!(3,0)(rg, qs[0 .. 3], ps[0 .. 1], 0.0f);
				case ProcessType.ARIMA: return new ARIMA!(3,1,0)(rg, qs[0 .. 3], ps[0 .. 1], 0.0f);

				//new AR!(1)(rg, 10.0, 1.0, [/*beta0=*/0.8f], 0.0); // converges to 50, and oscilates beetwen 45 and 55
				//new AR!(2)(rg, 10.0, 1.0, [/*beta0=*/0.8f, 0.1f], 0.0);
				//new AR!(1)(rg, 10.0, 1.0, [/*beta0=*/0.99f], 0.0); // converges to 1/(1-beta0)*y0
				//new AR!(1)(rg, 10.0, 1.0, [/*beta0=*/-0.8f], 0.0); // oscilates near 5
				// Brownian motion

				//new AR!(1)(rg, 0.0, 1.0, [0.999f], 0.0);
				//new MA!(1)(rg, 10.0, [0.1f, 0.1f]);
				//new MA!(2)(rg, 10.0, [0.1f, 0.2f, -0.3f]);
				//new MA!(2)(rg, 10.0, [0.2f, 0.1f, 0.05f]);
				//new MA!(2)(rg, 10.0, [1.0f, 0.2f, 0.1f]);
				//new ARMA!(2,1)(rg, [0.8f, 0.1f], [1.0f, 0.2f], 10.0);
//				case ProcessType.GARCH: ratio = true; return new GARCH!(1,1)(rg, [0.2f], [1.0f, 0.1], 0.0);
				case ProcessType.GARCH: ratio = true;
					//return new GARCH!(1,1)(rg, [0.8f], [1.0f, 0.17], 0.0);
					return new GARCH!(1,1)(rg, [0.2f], [1.0f, 0.2], 0.0);
				}
			});

		// the same generator for comparission/drawing in background
		//auto rg2 = new RandomsGauss(seed, 0.0f, 1.0f);

		static if (acki) {
			auto ac1 = new Autocorellation!(float, 100)();
		}

		//float yn;
		//size_t i;

		float total_length = length+forecast;

		float yy;
		genwraped_ar = new GeneratorProxy!(float)(genwraped_ar, delegate void(float x) {
				yy = x;
			});

		float y_before_invlir;

		if (show_invlir) {
			if (invlir_factor > 0.0) {
				assert(invlir_factor > 0.0);
				genwraped_ar = new GeneratorProxy!(float)(genwraped_ar, delegate void(float x) {
					y_before_invlir = x;
				});
				genwraped_ar = new InvLIR!(float)(genwraped_ar, 1.0, invlir_factor);
			}
			if (plots) {
				g.set("yrange [0:]");
			}
		} else {
			if (plots) {
				g.set("yrange [*:*]");
			}
		}

		float y_before_smooth;

		if (show_expsmooth) {
			genwraped_ar = new GeneratorProxy!(float)(genwraped_ar, delegate void(float x) {
				y_before_smooth = x;
			});
			assert(expsmooth_factor > 1.0);
			genwraped_ar = new ExponentialSmoothing!(float)(genwraped_ar, expsmooth_factor);
		}

		if (plots) {
//			g.set("title 'points of serie'");
			if (showgrid) {
				g.set("grid");
			} else {
				g.unset("grid");
			}

			g.set("xrange [0:"~to!(string)(total_length)~"]");

//			g.plot(["'-' w l ls 1"]);
		}

		float gwn;
		//if (shownoise || showratio) {
		if (rg_proxy !is null) {
			rg_proxy.set_monitor(delegate void(float x) {
				gwn = x;
			});
		}

		scope ff = new Foreach!(float, genwraped_ar, total_length)();

		size_t buffer_i = 0;

//		g.emit_start();
		g.shm_emit_start(0);
		g.shm_emit_start(1);
		g.shm_emit_start(2);
		g.shm_emit_start(4);

		size_t i = 0;
		float yn;
		foreach (y; ff) {
			yn = y;

			if (!show_expsmooth) {
				y_before_smooth = y; // set it for conviniance
			}

			static if (teski) {
				writeln(i++, " ", y);
			}

			static if (acki) {
				if (show_autocorel) {
					float for_ac;
					if (show_invlir) {
						for_ac = y_before_invlir;
					} else {
						for_ac = y_before_smooth; // y or y_before_smooth
					}
					if (pt == ProcessType.GARCH) {
						ac1 ~= for_ac*for_ac;
					} else {
						ac1 ~= for_ac;
					}
				}
			}

			if (save_buffer !is null && buffer_i < save_buffer.length) {
				save_buffer[buffer_i++] = y_before_smooth;
			}

			if (plots) {
				if (show_expsmooth) {
					g.shm_emit_data(0, "%d %f %f", i++, y, y_before_smooth);
				} else {
					//g.emit_data("%f", y);
					g.shm_emit_data(0, "%d %f", i++, y);
				}

				if (show_invlir) {
					g.shm_emit_data(4, "%f", y_before_invlir);
				}

				if (shownoise) {
					g.shm_emit_data(1, "%f", gwn);
				}
				if (ratio && showratio) {
					g.shm_emit_data(2, "%f", (gwn == 0.0 ? 1.0 : yy/gwn));
				}
			}
		}
		//g.emit_end();
		g.shm_emit_end(0);
		g.shm_emit_end(1);
		g.shm_emit_end(2);
		g.shm_emit_end(4);

		static if (acki) {
			//writeln("Corellations:");
			if (plots) {
				//g.set("title 'Corellations'");
				//g.plot(["'-' w p"]);
				//g.emit_start();
				g.shm_emit_start(3);
				if (show_autocorel) {
					foreach (k; 0 .. ac1.max) {
						//writefln("C1 %d = %f", k, ac1[k]);
						//g.emit_data("%d %f", k, ac1.normedIndex(k));
						g.shm_emit_data(3, "%d %f", k, ac1.normedIndex(k));
					}
				}
				//g.emit_end();
				g.shm_emit_end(3);
				//g.plot(["'/dev/shm/gnuplot_d_3' w p"]);
			}
		}

		if (plots) {
			auto plot0_points = "'/dev/shm/gnuplot_d_0' with "~(showpoints ? "lp" : "l")~" ls 1 title \"process\"";
			if (show_expsmooth) {
				plot0_points = "'/dev/shm/gnuplot_d_0' using 1:3 with "~(showpoints ? "lp" : "l")~" ls 3 title \"before smoothing\" ," ~ plot0_points;
			}
			auto plot1_noise = "'/dev/shm/gnuplot_d_1' with l title \"noise\"";
			auto plot2_sigma = "'/dev/shm/gnuplot_d_2' with l title \"sigma_i\"";
			auto plot3_autocorel = "'/dev/shm/gnuplot_d_3' with lp title \"autocorel_k\"";

			auto plot4_lir = "'/dev/shm/gnuplot_d_4' with "~(showpoints ? "lp" : "l")~" ls 6 title \"process (LIR)\"";
			//if (show_expsmooth) {
			//	plot4_lir = "'/dev/shm/gnuplot_d_4' using 1:3 with "~(showpoints ? "lp" : "l")~" ls 3 title \"before smoothing\" ," ~ plot0_points;
			//}


			g.set("key center top");

			g.unset("origin");
			g.unset("size");
			g.unset("bmargin");
			g.unset("tmargin");

			g.set("xrange [0:"~to!(string)(total_length)~"]");
			g.set("x2range [0:"~to!(string)(total_length)~"]");


			if (separateplots) {
				int rows = 1;
				if (shownoise) { rows++; }
				if (ratio && showratio) { rows++; }

				//g.set("multiplot layout "~to!(string)(rows)~",1");

				g.set("lmargin 9");
				g.set("rmargin 9");

				g.set("multiplot");

				g.set("mxtics 10");

				g.set("format x \"\"");
				g.unset("xtics");
				g.set("x2tics");
				g.set("mx2tics 10");
				g.set("format x2 \"%g\"");
				if (show_invlir) {
					g.set("origin 0,0.65");
					g.set("size 1,0.3");
				} else {
					g.set("origin 0,0.35");
					g.set("size 1,0.6");
				}
				g.set("bmargin 0");
				g.set("tmargin 0");
				g.set("ylabel \"y_i\" offset 1");
				g.plot([plot0_points]);
				g.unset("ylabel");
				g.set("format x2 \"\"");
				g.unset("x2tics");
				g.set("xtics");

				if (shownoise) {
					g.set("origin 0,0.25");
					g.set("size 1,0.1");
					g.set("bmargin 0");
					g.set("tmargin 0");
					g.set("format x \"\"");
					g.set("ylabel \"eta_i\" offset 1");
					g.set("yrange [*:*]");

					g.plot([plot1_noise]);
					g.unset("ylabel");
				}

				if (ratio && showratio) {
					g.set("format x");
					if (show_autocorel) {
						g.set("origin 0,0.15");
						g.set("size 1,0.1");
					} else {
						g.set("origin 0,0.05");
						g.set("size 1,0.2");
					}
					g.set("bmargin 0");
					g.set("tmargin 0");
					g.set("format x \"\"");
					g.set("ylabel \"sigma_i\" offset 1");
					g.unset("ytics");
					g.set("y2tics");
					g.set("yrange [:*]");
					if (showratio_logscale) {
						g.set("logscale y");
					}

					g.plot([plot2_sigma]);
					g.unset("logscale y");
					g.unset("ylabel");
					g.unset("y2tics");
					g.set("ytics");
				}


				if (show_invlir) {
//					g.set("format x");
					g.set("origin 0,0.35");
					g.set("size 1,0.3");
					g.set("bmargin 0");
					g.set("tmargin 0");
					g.set("format x \"\"");
					g.set("ylabel \"y_i\" offset 1");
					g.unset("ytics");
					g.set("y2tics");
					g.set("yrange [*:*]");
					g.set("xzeroaxis");

					g.plot([plot4_lir]);
					g.unset("xzeroaxis");
					g.unset("ylabel");
					g.unset("y2tics");
					g.set("ytics");

				}
				g.unset("mx2tics");
				g.unset("mxtics");

				if (show_autocorel) {
					g.set("format x");
					if (ratio && showratio) {
						g.set("origin 0,0.05");
						g.set("size 1,0.08");
					} else {
						g.set("origin 0,0.05");
						g.set("size 1,0.18");
					}
					g.set("bmargin 0");
					g.set("tmargin 0");
					g.set("format x");
					g.set("ylabel \"autocorel_k\" offset 1");
					g.set("yrange [:*]");
					g.set("xzeroaxis");

					g.set("xrange [0:"~to!(string)(autocorel_length)~"]");

					g.plot([plot3_autocorel]);
					g.unset("xzeroaxis");
					g.unset("ylabel");
					g.unset("xrange");
				}


				g.unset("multiplot");
			} else {
				auto what_plot = [plot0_points];
				if (shownoise) {
					what_plot ~= plot1_noise;
				}
				if (ratio && showratio) {
					what_plot ~= plot2_sigma;
				}
				if (show_autocorel) {
					what_plot ~= plot3_autocorel;
				}
				g.plot(what_plot);
			}
		} // plots

		h_end ~= yn;


	if (noise_stats) {
		writefln("min=%f max=%f", rg_minmax.min, rg_minmax.max);
		rg_moments.dump();
		//rg_hist.dump();
/+
		jk2.run2(&moments_compose!(MT), &moments_uncompose!(MT), &moments_finalize_mean!(MT));
		writefln("avg = %.15g", jk2.estimator);
		writefln("err = %.15g", jk2.error);
+/
	}


	if (show_autocorel) {
		writefln("original qs = %s", qs);
		writefln("original ps = %s", ps);
	}

	} // for M

	//h_end.dump();
}

void anal(Generator!(float) serie, Gnuplot g, bool lir = false, bool show_noise = false, bool show_autocorel = false) {


	if (lir) {
		serie = new LIR!(float)(serie);
	}

	auto ac2 = new Autocorellation!(float, 100)();

	if (0) {
		g.shm_emit_start(5);
		size_t i;
		foreach (y; serie) {
			if (show_noise) {
				g.shm_emit_data(5, "%d %f", i++, y);
			}
			ac2 ~= y;
		}
		g.shm_emit_end(5);


		if (show_noise) {
			g.plot(["'/dev/shm/gnuplot_d_5' u 1:2 w l"]);
		}

		if (show_autocorel) {
			g.shm_emit_start(6);
			foreach (k; 0 .. ac2.max) {
				g.shm_emit_data(6, "%d %f", k, ac2.normedIndex(k));
			}
			g.shm_emit_end(6);
			g.plot(["'/dev/shm/gnuplot_d_6' u 1:2 w l"]);
		}
	}


	writefln("yule walker");
	auto yw = yule_walker!(float)(1, ac_gen!(100)(ac2));
	writefln("yw(1) = %s", yw);
	yw = yule_walker!(float)(2, ac_gen!(100)(ac2));
	writefln("yw(2) = %s", yw);
	yw = yule_walker!(float)(2, ac_gen!(100)(ac2));
	writefln("yw(3) = %s", yw);

}


Generator!(float) ac_gen(uint N)(Autocorellation!(float, 100) ac) {
	return new class Generator!(float) {
		private int k = -1;
		override float getNext() { next(); return get(); }
		override float get() { return ac.normedIndex(k); }
		override bool next() { return (k++ < N); }
		mixin OpApplyMixin!(float) opApply;
	};
}


version (series_main_test) {
void main() {
	go();
	Thread.sleep(500_000_000);
}
}
