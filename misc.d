module misc;

/** One step of Kahan's sumation algorithm, initially s, and c should be 0 */
void kahan_update(T)(ref T s, T x_i, ref T c) {
	T y = x_i - c;
	//volatile;
	synchronized {};
	T t = s + y;
	//volatile;
	synchronized {};
	c = (t-s) - y;
	s = t;
}
