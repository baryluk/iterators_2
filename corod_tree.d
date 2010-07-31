/** This modules shows simple example of recursive iterator.
 *
 * Recursion of trees is good example of how easly it can be impleented.
 */
module corod_tree;

import corod;


/// Wrapper for tree
class Tree {
	TreeNode root; // root of tree

public
	/// create tree
	this(TreeNode root_) {
		root = root_;
	}

	// returns iterator
	// -1 - preorder
	// 0 - inorder
	// 1 - postorder
	TreeIterator getIterator(int p) {
		return root.getIterator(p);
	}


private:
	/// called from TreeIterator to start iteration
	static void iter(TreeIterator ti, int p, TreeNode s) {
		switch (p) {
			case 0: s.inorder_iter(ti); break;
			case -1: s.preorder_iter(ti); break;
			case 1: s.postorder_iter(ti); break;
		}
	}
}

class TreeNode {
	int v; /// some value in tree

private:
	TreeNode left;
	TreeNode right;

public:
	/// create node with two subtrees
	this(int v_, TreeNode left_, TreeNode right_) {
		v = v_;
		left = left_;
		right = right_;
	}

	/// create leaf
	this(int v_) {
		v = v_;
	}

private:
	/// iterate in in-order
	void inorder_iter(TreeIterator ti) {
		if (left !is null) { left.inorder_iter(ti); }
		ti.yield(v);
		if (right !is null) { right.inorder_iter(ti); }
	}

	/// iterate in pre-order
	void preorder_iter(TreeIterator ti) {
		ti.yield(v);
		if (left !is null) { left.preorder_iter(ti); }
		if (right !is null) { right.preorder_iter(ti); }
	}

	/// iterate in post-order
	void postorder_iter(TreeIterator ti) {
		if (left !is null) { left.postorder_iter(ti); }
		if (right !is null) { right.postorder_iter(ti); }
		ti.yield(v);
	}

	TreeIterator getIterator(int p) {
		return this.new TreeIterator(p);
	}

public:
	/// iterator as nested class
	class TreeIterator : FiberGenerator!(int) {
	private:
		int p;
	public:
		this(int p_) {
			p = p_;
		}

	protected:
		void iter() {
			Tree.iter(this, p, this.outer);
		}
	}
}



import std.stdio;

void main(string[] args) {
	auto tree = new Tree(
	new TreeNode(
		15,
		new TreeNode(
			33,
			new TreeNode(
				14,
				new TreeNode(7),
				new TreeNode(11)
			),
			new TreeNode(
				6,
				new TreeNode(74),
				new TreeNode(71)
			)
		),
		new TreeNode(
			56,
			new TreeNode(24, null, null),
			new TreeNode(
				3,
				null,
				new TreeNode(55)
			)
		)
	)
	);

	foreach (v; tree.getIterator(-1)) {
		writefln("v = %d", v);
	}
}


/// todo: permuttaions

