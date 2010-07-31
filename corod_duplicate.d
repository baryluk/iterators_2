module corod_duplicate;

import corod;

/** This is simplified case of 2 consumers */
class Remember2limit(T) {
private:
	T[] buffer;
	int i, j;
	int p1, p2;



	invariant() {
		assert(i != j);
		assert(p1 == i || p2 == i);
		if (p1 == p2) {
			assert((i+1) % buffer.length == j);
		}
	}

public:
	this(int n) {
		buffer = new T[n];
		i = 0;
		j = 0;
		p1 = 0;
		p2 = 0;
	}

	int S() const { return buffer.length; }

	/// rank=1|2
	T getNext(int rank, T delegate(int) takeNewDG)
	in {
		assert(rank == 1 && rank == 2);
	}
	body {
		T v = buffer[(rank == 1 ? p1 : p2)];
		if (i < j) {
			if (rank == 1) {
				if (p1 == i) { // p1 is oldest
					p1 = (p1 + 1) % S;
					if (p2 == i) {
						; //nothing p2 is still referencing i
					} else {
						i = (i + 1) % S;
					}
				} else { // p2 is oldest
					assert(p2 == i);
					p1 = (p1 + 1) % S;
				}
			} else {
				assert(rank == 2);
				if (p2 == i) { // p2 is oldest
					p2 = (p2 + 1) % S;
					if (p1 == i) {
						; //nothing p1 is still referencing i
					} else {
						i = (i + 1) % S;
					}
				} else { // p1 is oldest
					assert(p1 == i);
					p2 = (p2 + 1) % S;
				}
			}
		} else {
			
		}
		return v;
	}
}

/** Remembers history of values, but not indefinitly,
 * but deffinied by multiple index pointers.
 *
 * It remember as many values as needed from stream of values,
 * depending how many stream readers have already readed.
 * If so values are not needed anymore (all readers are after this value)
 * it is removed from buffer.
 *
 * It is essentially just like having multiple processes
 * have indexes i_rank >= 0, indexint table normal_buffer;
 * They can only do one operation: normal_buffer[++i_rank].
 * If i_rank+1 is valid index for normal_buffer,
 * this operations returns value at index i_rank+1, and advances i_rank by 1,
 * If it's not valid, first normal_buffer is enlarged if needed,
 * to have space for new values, then external function getNext()
 * is called, and its return value is saved in normal_buffer[i_rank+1],
 * this value is also returned to initial process and i_rank is advanced by 1.
 *
 * This way multiple processes/consumers can consume from single stream.
 *
 * What additionally does this class? It uses as small space
 * as possible, utilizing dynamically resized circular buffer,
 * and keeping values of i_rank for each rank, as also
 * it's minimal value (and coresponding rank). If no process
 * can consume position p (because min{i_rank} >= p, so \all i_rank >= p),
 * then p is removed from buffer.
 *
 * This ensures all operations are O(1) (eventually finding new min{i_rank}
 * can be O(max{rank}), but it is sufficient to maintain order on ranks, to
 * ensure O(1) aso in this case).
 *
 * In mosts cases this will also lead to constant size of buffers,
 * if all consumers read data from stream in roughly the same rate.
 *
 * Remark: If max{i_rank}-min{i_rank} have approximetly growing behaviour,
 *     (even linear) then underling buffer will grow indefinitly.
 *      Keep this difference bounded and as small as possible
 *      (in practice it rearly exceads 5).
 *
 * Important characteristic:
 *   If you have created RememberN!(T) r, with n as number of processes,
 *   and then created table:
 * -----
 *        int[n] k; T[n][S] vals;
 * -----
 *       and performed:
 * -----
 *      while (k.max < S) {
 *         auto rank = random(n).next;
 *         if (k[rank] < S) {
 *           vals[k[rank]++] = r.getNext(rank, anydelegate); // anydelegate not modifing k, rank, vals, S and n
 *         }
 *      } // this loop eventually terminates because random is uniform generator
 *        // random can be any permutation of (0, ..., 0, 1, ..., 1, ..., ..., n-1, ..., n-1) sequence
 *                                             \__.___/
 *                                              S-times
 *        // it also terminats, because each k[rank] is then incremented S times, so after processing
 *           of all elements k.max == S.
 * -----
 *     then, after loop it is true:
 * -----
 *        foreach (i; 0 .. S) {
 *           T v = vals[i][0];
 *           foreach (r; 0 .. n) {
 *               assert(vals[i][r] == v);
 *           }
 *        }
 * -----
 *  in short words: getNext(rank, dg), will return the same values in order for each valid rank independly.
 *
 * TODO: speciall case for n=1 (synhronization, and costly tables not needed)
 *
 * Note: dg should never return null (null is used as marker in buffer), and
 *  all valid entries are assumed to be notnull (for verification).
 */
class RememberN(T) {
private:
	T[] buffer; // cyclic buffer used for storing history
				// in case of overflow it is reallocated, size is doubled,
				// data copied into begging region.
				// Note: buffer is never shrinken.
	int i; // where is begging of circular buffer
	int j; // where is end of circular buffer

public:
	invariant() {
		assert(starting_position_in_buffer >= 0);
		assert(min_rank > -1);
		assert(min_position >= 0);

		// it is stupid to remember anything when nobody is interested in it
		assert(positions.length > 0);

		assert(min_rank < positions.length);

		if (starting_position_in_buffer > 0) {
			assert(min_rank >= 0);
		}

		assert(0 < i && i < buffer.length);
		assert(0 < j && j < buffer.length);
		assert(i != j);

		// todo: unused part of buffer contains null values or T.init
		// this ensures that old data was correctly discarded

		foreach (p; positions) {
			assert(p >= 0);
			assert(p >= min_position);
			assert(p <= starting_position_in_buffer + buffer.lentgh); // weak bound
		}

		assert(ranks_ordered_by_possion.length == positions.length);

		assert(rank_ordered_by_possition[rank_ordered_by_position.length-1]
			<= starting_position_in_buffer + buffer.lentgh); // weak bound

		foreach (i, r; ranks_ordered_by_position[0 .. $-1]) {
			assert(r >= 0);
			assert(r < positions.length);
			assert(positions[r] <= positions[ranks_ordered_by_position[i+1]]);
		}

		if (min_rank >= 0) {
			assert(ranks_ordered_by_position[0] == min_rank);
			assert(positions[min_rank] == min_position);
		}

		assert(min_position == starting_position_in_buffer);

		assert(position_of_rank_in_ranks_ordered_by_position.length == positions.length);

		foreach (r, por; position_of_rank_in_ranks_ordered_by_position) {
			assert(ranks_ordered_by_position[por] == r);
		}

		// todo: position_of_rank_in_ranks_ordered_by_position is permutation of [0, ..., n-1]
	}

private:
	int starting_position_in_buffer; /// indicates what are the oldest data in the begining of the circular buffer

	int[] positions; /// position of each itterator
					/// position[rank] is a value which will be return by getNext(rank)
					// position[rank] can hypothycally, point beyond the buffer,
					//on getNext all missing values will be generated
	int min_position; /// minimal value from positions
	int min_rank; /// rank of iterator which have minimal possiton

	/// indexes to position[],
	/// so position[ranks_ordered_by_position[i]] is not decressing.
	int[] ranks_ordered_by_position;
	/// indexes to ranks_ordered_by_position,
	/// so ranks_ordered_by_position[position_of_rank_in_ranks_ordered_by_position[i]] == i, for all i
	int[] position_of_rank_in_ranks_ordered_by_position;

public:
	/// construct n new generators based on this
	this(int n_) {
		positions = new int[n_];
		positions[] = 0;
		min_position = 0;
		min_rank = 0;
		ranks_ordered_by_position = new int[n_];
		foreach (i; 0 .. n_) {
			ranks_ordered_by_position[i] = i;
		}
		position_of_rank_in_ranks_ordered_by_position = new int[n_];
		foreach (i; 0 .. n_) {
			position_of_rank_in_ranks_ordered_by_position[ranks_ordered_by_position[i]] = i;
		}
		i = 0;
		j = 1;
		starting_position_in_buffer = 0;
		buffer = new int[2];
	}

	/// last valid position
	private int ending_position_in_buffer() {
		if (i < j) {
			return starting_position_in_buffer + (j-i);
		} else {
			return starting_position_in_buffer + bl - (j-i);
		}
	}

	/// returns size of buffer
	private int bl()
	out(ret) {
		assert(bl > 0);
	}
	body {
		return buffer.length;
	}

	/// number of duplicated generators
	private int n()
	out (res) {
		assert(res > 0);
	}
	body {
		return positions.length;
	}

	/// internal function for retriving position t
	private T get(int t)
	out(res) {
		assert(res !is T.init);
	}
	body {
		assert(!(t < starting_position_in_buffer));
		assert(!(t > ending_position_in_buffer()));
		if (i < j) {
			auto tt = t-starting_position_in_buffer+i;
			assert(tt >= i);
			if (tt < bl) {
				return buffer[tt];
			} else {
				auto tt2 = tt-bl;
				assert(tt2 < bl);
				assert(tt2 >= 0);
				assert(tt2 < i);
				return buffer[tt2];
			}
		}
	}

	/** retrivies value at possition t from buffer
	 *
	 * Note: it must be in buffer
	 */
	private T opIndex(int t)
	out(res) {
		assert(res !is T.init);
	}
	body {
		if (t < starting_position_in_buffer) {
			throw new Exception("requested position is already forgoten");
		}
		if (t > ending_position_in_buffer()) {
			throw new Exception("requested position is not yet in buffer");
		}
		return get(t);
	}

	/** like opIndex but if t is grater than last valid possition,
	 * we take new values from stream until it will be valid
	 *
	 * Note: it is valid to pass value which exceeds last valid position,
	 *       greatly (not just by 1)
	 */
	private T opIndexAuto(int t, T delegate() takeNewDG)
	out(res) {
		assert(res !is T.init);
	}
	body {
		if (t < starting_position_in_buffer) {
			throw new Exception("requested position is already forgoten");
		}
		if (t <= ending_position_in_buffer()) {
			return get(t);
		}

		if (i != j) { // buffer have space
			;
		} else { // buffer full
			// buffer extension needed
			int oldn = n();
			buffer.length = buffer.length * 2;
			//if (!(i < j)) { // in fact, we are in else because i == j!
				buffer[j .. j+i] = buffer[0 .. i];
				buffer[0 .. i] = T.init;
				j = j+i;
			//}
			int newn = n();
			assert(oldn == newn);
		}

		assert(i != j);

		j += 1;
		if (j > bl()) {
			j -= bl();
			//j = 0; // equivalent
		}
		assert(0 < j && j <= bl); // j points one position after last element, so there we can write

		// needed to aquire new element.
		T new_elem = takeNewDG(t);

		if(j == bl) {
			buffer[0] = new_elem;
		} else {
			buffer[j] = new_elem;
		}

		return new_elem;
	}

	/** Returns value of buffer.
	 *
	 * Increments position of rank by one.
	 */
	public T getNext(int rank, T delegate(int) takeNewDG)
	in {
		assert(rank >= 0);
		assert(rank < n);
	}
	out (res) {
		assert(res !is T.init);
	}
	body {
		auto k = positions[rank];
		auto v = opIndexAuto(k, takeNewDG);
		positions[rank]++;
		auto new_k = positions[rank];
		auto p = position_of_rank_in_ranks_ordered_by_position[rank];
		if (p < n) { // if we ware newest, so we are after this getNext
			// rearangments needed
			for (int p2 = p+1; p2 < n; p2++) { // actually this loop, will not execture if (p<n)
				/// Todo: if (position[ranks_ordered_by_position[
				assert(0);
			}
		}
		return v;
	}

	/** returns position of iterator rank */
	public int myPosition(rank)
	in {
		assert(rank >= 0);
		assert(rank < n);
	}
	out(res) {
		assert(res >= 0 && res < ending_position_in_buffer());
	}
	body {
		assert(rank >= 0);
		assert(rank < n);
		return positions[rank];
	}
}

/** Duplicator class, uses RememberN as implementation */
class Duplicator(T) {
	RememberN!(T) hist;
	Generator!(T) org;
	Generator!(T)[] duplicated;

	/** Construct n new generators from generator org_ */
	this(Generator!(T) org_, int n) {
		org = org_;
		hist = new RememberN!(T)(n);
		foreach (i, ref d; duplicated) {
			d = new Duplicated!(T)(this, i);
		}
	}

	/** Take generator of rank i */
	Duplicated!(T) opIndex(int i) {
		return duplicated[i];
	}

private:
	/** Run getNext on generator of rank 'rank'
	 *
	 * Used in Duplicated!(T) class
	 */
	T getNext(int rank) {
		hist.getNext(rank, (int new_position) { return org.getNext(); } ); // lazy delegate
	}
}

/** Duplicated generator, this class is returned by opIndex(int) method in Duplicator. */
private class Duplicated(T) : Generator!(T) {
	const int myrank;
	Duplicator!(T) duplicator;

private:
	/// constructor
	this(Duplicator!(T) duplicator_, int myrank_) {
		myrank = myrank_;
		duplicator = duplicator_;
	}

protected:
	/// implementation of iteration
	void iter() {
		while (true) {
			yield(duplicator.getNext(myrank));
		}
	}
}


/** Buffering generator,
 *
 * It first reads n new values to buffer, then it yield n values from this buffer, and then repeats.
 */
class BufferingGenerator(T) : Generator!(T) {
private:
	T[] buffer;
	Generator!(T) G;

	invariant() {
		assert(i <= buffer.length);
	}

public:
	/// construct buffer of size.
	this(Generator!(T) G_, int size) {
		G = G_;
		buffer = new T[size];
	}

protected:
	/// implemntation
	void iter() {
		while (true) {
			foreach (ref b; buffer) {
				b = G.getNext();
			}
			foreach (b; buffer) {
				yield(b);
			}
		}
	}
}
