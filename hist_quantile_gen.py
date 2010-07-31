#!/usr/bin/python

# utility to compute b, and k for usage in 98sigmod-quantiles.pdf


import sys

N = int(float(sys.argv[1])) # 1e5
epsilon = float(sys.argv[2]) # 0.01


#def binom(n,k):
#	return reduce(lambda a,b: a*(n-b)/(b+1),xrange(k),1)

def binom(n, k):
	if 0 <= k <= n:
		p = 1
		for t in xrange(min(k, n - k)):
			p = (p * (n - t)) // (t + 1)
		return p
	else:
		return 0

def f1(b, h):
	#print (h-2)*binom(b+h-2, h-1) - binom(b+h-3, h-3) + binom(b+h-3, h-2), int(2*epsilon*N)
	return ((h-2)*binom(b+h-2, h-1) - binom(b+h-3, h-3) + binom(b+h-3, h-2) <= int(2*epsilon*N))

def f2(b, h, k):
	return k*binom(b+h-2, h-1) >= N

min_bk = 100000000000000000
min_bki = []

for b in xrange(2, 40):
	h = 3
	max_h = 0
	while f1(b, h):
		#print b, h
		max_h = h
		h += 1
	h = max_h
	if max_h > 0:
		k = int(N / binom(b+h-2, max_h-1))+1
		if f2(b, h, k) == False:
			print "error"
		bk = b*k
		#print "b=%d h=%d k=%d, kd=%d" % (b,max_h,k,bk)
		if bk < min_bk:
			min_bk = bk
			min_bki = [(b,k)]
		elif bk == min_bk:
			min_bki.append((b,k))

print min_bki, (min_bki[0][0]*min_bki[0][1])*4/1024

