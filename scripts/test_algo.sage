## Parameters

#D,q = 16, 97
#D,q = 16, 47
#D,q = 24, 73
#D,q = 24, 71

#D,q = 512,  12289
#D,q = 512,   5119
#D,q = 512,   3583
#D,q = 768,  18433
#D,q = 768,   6911
#D,q = 768,   3329
#D,q = 1024, 12289
D,q = 1024, 5119


## Derived constants and preliminary tests

assert q < 2**15

assert all(f in [2,3] for f in dict(factor(D)))
M = 3*D if D % 3 == 0 else 2*D
assert D == euler_phi(M)

phi_M = sage.rings.polynomial.cyclotomic.cyclotomic_coeffs(M, sparse=True)
assert phi_M == ({D: 1, D//2: -1, 0: 1} if D % 3 == 0 else {D: 1, 0: 1})
assert max(phi_M.keys()) == D
assert phi_M[D] == 1
del phi_M[D]

prim_pow = [n for n in range(M) if gcd(n,M) == 1]

tau = 2*pi
def root_fft(k, m):
	return CDF(exp(I*tau*k/m))


assert q in Primes()
try:
	ext = next(ext for ext in [1,2] if (q^ext-1) % M == 0)
except StopIteration:
	print("q value is invalid (M does not divide q^k-1 for k=1,2)")
	exit(1)
if ext == 2 and q % (6 if D % 3 == 0 else 4) == 1:
	print("q value is inconvenient (ext=2 and q % 4/6 = 1)")
	exit(1)

fq = FiniteField(q)
fqx.<alpha> = FiniteField((q,ext))
g_ntt = fqx.zeta(M)
def root_ntt(k, m):
	return g_ntt^(M//m * k)

print(f"D={D}, q={q} -> M={M}, rem={phi_M}, ext={ext}, g_ntt={g_ntt}")


## Utility functions

def eval(p, x):
	assert len(p) == D
	x_k = 1
	s = 0
	for k in range(D):
		s += p[k] * x_k
		x_k *= x
	return s

def poly_mod_mul(p1, p2):
	assert len(p1) == len(p2) == D
	p3 = [0]*(2*D-1)
	for i1, x1 in enumerate(p1):
		for i2, x2 in enumerate(p2):
			p3[i1+i2] += x1*x2
	# modulo phi_M
	for i1 in range(2*D-2, D-1, -1):
		c1 = p3[i1]
		for i2, c2 in phi_M.items():
			p3[i1 % D + i2] -= c1 * c2
	return p3[:D]

def tbt_mul(v1, v2):
	return [z1 * z3 for z1, z3 in zip(v1, v2)]


## FFT tests

def cpx_eq(x, y):
	return abs(x-y) < 1e-9
def cpx_vec_eq(v1, v2):
	assert len(v1) == len(v2)
	return all(cpx_eq(x,y) for x,y in zip(v1, v2))

def print_vecs(*vecs):
	for v in vecs:
		print([
			round(z,3) if isinstance(z, float) else
			z.n(digits=3) if isinstance(z, ComplexDoubleElement) or isinstance(z, RealDoubleElement) else
			z
			for z in v
		])

import random
def rand_real_vec():
	return [random.random() for _ in range(D)]

def cpx_dft(p):
	assert len(p) == D
	return [eval(p, root_fft(k, M)) for k in prim_pow]

def test_fft(fft, ifft, fmt=lambda x: x, mul=tbt_mul):
	# test fft and compare with dft
	p1 = rand_real_vec()
	#p1 = list(range(D))
	v1 = fft(p1)
	vdft = fmt(cpx_dft(p1))
	assert cpx_vec_eq(v1, vdft)
	#print_vecs(v1)

	print(fft.__name__, "passes tests")
	if ifft is None: return
	
	# test ifft and round trip
	prt = ifft(v1)
	assert cpx_vec_eq(p1, prt)
	
	# multiplication via fft
	p2 = rand_real_vec()
	v2 = fft(p2)
	v3 = mul(v1, v2)
	p3 = ifft(v3)
	
	# compare with naive multiplication
	p4 = poly_mod_mul(p1, p2)
	assert cpx_vec_eq(p3, p4)
	
	print(ifft.__name__, "passes tests")

# evaluates polynomial p of degree < d = euler_phi(m)
# at (w_i) the roots of phi_m (the d primitive m-th roots of unity)
# recursive, complex to complex
def fft1(p, d=D, m=M):
	assert len(p) == d == euler_phi(m)
	if m == 2: # (d = 1)
		return p
	elif m == 6: # (d = 2)
		# with m = 6, the 2 roots are not opposites...
		# they are conjugates though.
		# * p(w_0) = p0 + w_0 p1
		# * p(w_1) = p0 + w_1 p1 = p0 + conj(w_0 p1)
		a, b = p[0], root_fft(1,6) * p[1]
		return [a + b, a + b.conjugate()]
	elif d % 3 == 0:
		p0, p1, p2 = p[::3], p[1::3], p[2::3]
		# p(x) = p0(x³) + x p1(x³) + x² p2(x³)
		v0 = fft1(p0, d//3, m//3)
		v1 = fft1(p1, d//3, m//3)
		v2 = fft1(p2, d//3, m//3)
		# pj have degree d/3 → pj(w_i³) = vj[i]
		# p(w_(d/2+i)) = p(-w_i)
		v = [0]*d
		j = CDF(exp(I*tau/3))
		j2 = j*j
		for i in range(d//3):
			w = root_fft(prim_pow[i], m)
			a, b, c = v0[i], w * v1[i], (w*w) * v2[i]
			v[i] = a + b + c
			v[i+d//3] = a + j*b + j2*c
			v[i+2*d//3] = a + j2*b + j*c
		return v
	elif d % 2 == 0:
		p0, p1 = p[::2], p[1::2]
		# p(x) = p0(x²) + x p1(x²)
		v0 = fft1(p0, d//2, m//2)
		v1 = fft1(p1, d//2, m//2)
		# pj have degree d/2 → pj(w_i²) = vj[i]
		# p(w_(d/2+i)) = p(-w_i)
		v = [0]*d
		for i in range(d//2):
			w = root_fft(prim_pow[i], m)
			a, b = v0[i], w * v1[i]
			v[i] = a + b
			v[i+d//2] = a - b
		return v
	else:
		raise NotImplementedError()

def ifft1(v, d=D, m=M):
	assert len(v) == d == euler_phi(m)
	if m == 2: # (d = 1)
		return v
	elif m == 6: # (d = 2)
		a, b = v[0], v[1]
		# a = p0 + w0 p1, b = p0 + w1 p1
		# a + b = 2p0 + (w0+w1)p1 = 2p0 + p1
		# a - b = (w0-w1)p1
		dr = root_fft(1,6)-root_fft(5,6)
		p1 = (a - b) / dr
		p0 = (a + b - p1) / 2
		return [p0, p1]
	elif d % 3 == 0:
		v0, v1, v2 = [0]*(d//3), [0]*(d//3), [0]*(d//3)
		j = CDF(exp(I*tau/3))
		j2 = j*j
		for i in range(d//3):
			w = root_fft(prim_pow[i], m)
			a, b, c = v[i], v[i+d//3], v[i+2*d//3]
			v0[i] = (a + b + c) / 3
			v1[i] = (a + j2*b + j*c) / (3*w)
			v2[i] = (a + j*b + j2*c) / (3*w*w)
		p0 = ifft1(v0, d//3, m//3)
		p1 = ifft1(v1, d//3, m//3)
		p2 = ifft1(v2, d//3, m//3)
		return [c for t in zip(p0,p1,p2) for c in t]
	elif d % 2 == 0:
		v0, v1 = [0]*(d//2), [0]*(d//2)
		for i in range(d//2):
			w = root_fft(prim_pow[i], m)
			a, b = v[i], v[i+d//2]
			v0[i] = (a + b)/2
			v1[i] = (a - b)/2 / w
		p0 = ifft1(v0, d//2, m//2)
		p1 = ifft1(v1, d//2, m//2)
		return [c for t in zip(p0,p1) for c in t]
	else:
		raise NotImplementedError()

# recursive, real to complex
# because p has real coefficients, p(conj(x)) = conj(p(x))
# so we can restrict the evaluation to (w'_i = w_2i, 0≤i<d/2)
# because we have: w_(d-1-i) = conj(w_i)
# we can recover the evaluations on the rest of the roots with:
# p(w_(2i+1)) = p(conj(w_(d-2i-2))) = conj(p(w'_(d/2-i-1)))
# note: interestingly, this doesn't quite work for the cyclic convolution
# note: -w'_i = w'_(i+d/4)
def fft2(p, d=D, m=M):
	assert len(p) == d == euler_phi(m)
	#print(f"d={d}, m={m}")
	if m == 4: # (d = 2)
		# evaluate in z=i
		# could be a no-op depending on layout!
		return [p[0] + I*p[1]]
	elif m == 6: # (d = 2)
		# definitely not a no-op
		return [p[0] + root_fft(1,6) * p[1]]
	elif d % 3 == 0:
		p0, p1, p2 = p[::3], p[1::3], p[2::3]
		v0 = fft2(p0, d//3, m//3)
		v1 = fft2(p1, d//3, m//3)
		v2 = fft2(p2, d//3, m//3)
		v = [0]*(d//2)
		j = root_fft(1,3)
		j2 = j*j
		for i in range(d//6):
			w = root_fft(prim_pow[2*i], m)
			a, b, c = v0[i], w * v1[i], (w*w) * v2[i]
			v[i] = a + b + c
			v[i+d//6] = a + j*b + j2*c
			v[i+2*d//6] = a + j2*b + j*c
		return v
	elif d % 2 == 0:
		p0, p1 = p[::2], p[1::2]
		v0 = fft2(p0, d//2, m//2)
		v1 = fft2(p1, d//2, m//2)
		v = [0]*(d//2)
		for i in range(d//4):
			w = root_fft(prim_pow[2*i], m)
			a, b = v0[i], w * v1[i]
			v[i] = a + b
			v[i+d//4] = a - b
		return v
	else:
		raise NotImplementedError()

def ifft2(v, d=D, m=M):
	assert len(v)*2 == d == euler_phi(m)
	if m == 4: # (d = 2)
		return [v[0].real(), v[0].imag()]
	elif m == 6: # (d = 2)
		w = root_fft(1,6)
		p1 = v[0].imag() / w.imag()
		p0 = v[0].real() - p1*w.real()
		return [p0, p1]
	elif d % 3 == 0:
		v0, v1, v2 = [0]*(d//6), [0]*(d//6), [0]*(d//6)
		j = root_fft(1,3)
		j2 = j*j
		for i in range(d//6):
			w = root_fft(prim_pow[2*i], m)
			a, b, c = v[i], v[i+d//6], v[i+2*d//6]
			v0[i] = (a + b + c) / 3
			v1[i] = (a + j2*b + j*c) / (3*w)
			v2[i] = (a + j*b + j2*c) / (3*w*w)
		p0 = ifft2(v0, d//3, m//3)
		p1 = ifft2(v1, d//3, m//3)
		p2 = ifft2(v2, d//3, m//3)
		return [c for t in zip(p0,p1,p2) for c in t]
	elif d % 2 == 0:
		v0, v1 = [0]*(d//4), [0]*(d//4)
		for i in range(d//4):
			w = root_fft(prim_pow[2*i], m)
			a, b = v[i], v[i+d//4]
			v0[i] = (a + b)/2
			v1[i] = (a - b)/2 / w
		p0 = ifft2(v0, d//2, m//2)
		p1 = ifft2(v1, d//2, m//2)
		return [c for t in zip(p0,p1) for c in t]
	else:
		raise NotImplementedError()

def digit_rev(x, n):
	digits = []
	while n % 3 == 0:
		digits.append((x % 3, 3))
		x //= 3
		n //= 3
	while n % 2 == 0:
		digits.append((x % 2, 2))
		x //= 2
		n //= 2
	assert x == 0 and n == 1
	for d, r in digits:
		x = x*r + d
	return x

def fft3(p):
	assert len(p) == D
	v = p.copy()
	d = 1
	m = M//D # 2 or 3
	io = D//2 # offset between real and imaginary parts
	while m < M:
		if D % (2*d) == 0:
			d, m = 2*d, 2*m
		elif D % (3*d) == 0:
			d, m = 3*d, 3*m
		else:
			raise NotImplementedError()

		s = D//d # step
		for o in range(s): # offset
			if m == 4: # (d = 2)
				# evaluate in z=i
				# z = v[o] + I*v[o+s]
				# v[o], v[o+io] = z.real(), z.imag()
				# this is a no-op! (s=io)
				pass

			elif m == 6: # (d = 2)
				z = v[o] + root_fft(1,6)*v[o+s]
				v[o], v[o+io] = z.real(), z.imag()

			elif d % 3 == 0:
				j = root_fft(1,3)
				j2 = j*j
				for i in range(d//6):
					w = root_fft(prim_pow[digit_rev(i, d//3)], m)
					i0 = o + 3*s*i
					v0 = v[i0]     + I*v[i0+io]
					v1 = v[i0+s]   + I*v[i0+s+io]
					v2 = v[i0+2*s] + I*v[i0+2*s+io]

					a, b, c = v0, w*v1, w*w*v2
					va = a + b + c
					vb = a + j*b + j2*c
					vc = a + j2*b + j*c

					v[i0],     v[i0+io]     = va.real(), va.imag()
					v[i0+s],   v[i0+s+io]   = vb.real(), vb.imag()
					v[i0+2*s], v[i0+2*s+io] = vc.real(), vc.imag()

			elif d % 2 == 0:
				for i in range(d//4):
					w = root_fft(prim_pow[digit_rev(i, d//2)], m)
					i0 = o + 2*s*i
					v0 = v[i0]   + I*v[i0+io]
					v1 = v[i0+s] + I*v[i0+s+io]

					a, b = v0, w * v1
					va = a + b
					vb = a - b

					v[i0],   v[i0+io]   = va.real(), va.imag()
					v[i0+s], v[i0+s+io] = vb.real(), vb.imag()

			else:
				raise NotImplementedError()

	return v

def ifft3(v):
	assert len(v) == D
	v = v.copy()
	d = D
	m = M
	io = D//2 # offset between real and imaginary parts
	while d > 1:
		s = D//d # step
		for o in range(s): # offset
			if m == 4: # (d = 2)
				# still a no-op!
				pass

			elif m == 6: # (d = 2)
				z = v[o] + I*v[o+io]
				w = root_fft(1,6)
				v[o+s] = z.imag() / w.imag()
				v[o] = z.real() - v[o+s]*w.real()

			elif d % 3 == 0:
				j = root_fft(1,3)
				j2 = j*j
				for i in range(d//6):
					w = root_fft(prim_pow[digit_rev(i, d//3)], m)
					i0 = o + 3*s*i
					va = v[i0]     + I*v[i0+io]
					vb = v[i0+s]   + I*v[i0+s+io]
					vc = v[i0+2*s] + I*v[i0+2*s+io]

					v0 = (va + vb + vc)
					v1 = (va + j2*vb + j*vc) / w
					v2 = (va + j*vb + j2*vc) / (w*w)

					v[i0],     v[i0+io]     = v0.real(), v0.imag()
					v[i0+s],   v[i0+s+io]   = v1.real(), v1.imag()
					v[i0+2*s], v[i0+2*s+io] = v2.real(), v2.imag()

			elif d % 2 == 0:
				for i in range(d//4):
					w = root_fft(prim_pow[digit_rev(i, d//2)], m)
					i0 = o + 2*s*i
					va = v[i0]   + I*v[i0+io]
					vb = v[i0+s] + I*v[i0+s+io]

					v0 = (va + vb)
					v1 = (va - vb) / w

					v[i0],   v[i0+io]   = v0.real(), v0.imag()
					v[i0+s], v[i0+s+io] = v1.real(), v1.imag()

			else:
				raise NotImplementedError()
		
		if d % 3 == 0:
			d, m = d//3, m//3
		elif d % 2 == 0:
			d, m = d//2, m//2
		else:
			raise NotImplementedError()
	
	for i in range(D):
		v[i] /= D//2

	return v

def split_complex(l):
	return [z.real() for z in l] + [z.imag() for z in l]
def mul_complex(v1, v2):
	io = D//2
	v3 = [0]*D
	for i in range(io):
		x = v1[i] + I*v1[i+io]
		y = v2[i] + I*v2[i+io]
		z = x * y
		v3[i], v3[i+io] = z.real(), z.imag()
	return v3


## NTT tests

def rand_int_vec():
	return [fq.random_element() for _ in range(D)]

def int_dft(p):
	assert len(p) == D
	return [eval(p, root_ntt(k, M)) for k in prim_pow]

def test_ntt(ntt, intt, fmt=lambda x: x, mul=tbt_mul):
	# test ntt and compare with cpx_dft
	p1 = rand_int_vec()
	v1 = ntt(p1)
	vdft = fmt(int_dft(p1))
	#print_vecs(v1, vdft)
	assert v1 == vdft

	print(ntt.__name__, "passes tests")
	if intt is None: return
	
	# test intt and round trip
	prt = intt(v1)
	#print_vecs(p1, prt)
	assert p1 == prt

	# multiplication via fft
	p2 = rand_int_vec()
	v2 = ntt(p2)
	v3 = mul(v1, v2)
	p3 = intt(v3)
	
	# compare with naive multiplication
	p4 = poly_mod_mul(p1, p2)
	assert p3 == p4
	
	print(intt.__name__, "passes tests")

def ntt2(p, d=D, m=M):
	assert len(p) == d == euler_phi(m)
	if d == 1:
		return [p[0]]
	elif d == 2 and ext == 2:
		assert root_ntt(prim_pow[0],m).conjugate() == root_ntt(prim_pow[1],m)
		return [p[0] + root_ntt(1,m)*p[1]]
	elif m == 6:
		return [
			p[0] + root_ntt(1,6)*p[1],
			p[0] + root_ntt(5,6)*p[1],
		]
	elif d % 3 == 0:
		p0, p1, p2 = p[::3], p[1::3], p[2::3]
		v0 = ntt2(p0, d//3, m//3)
		v1 = ntt2(p1, d//3, m//3)
		v2 = ntt2(p2, d//3, m//3)
		hd = d//ext//3
		v = [0]*(3*hd)
		j = root_ntt(1,3)
		j2 = j*j
		for i in range(hd):
			w = root_ntt(prim_pow[ext*i], m)
			a, b, c = v0[i], w * v1[i], (w*w) * v2[i]
			v[i] = a + b + c
			v[i+hd] = a + j*b + j2*c
			v[i+2*hd] = a + j2*b + j*c
		return v
	elif d % 2 == 0:
		p0, p1 = p[::2], p[1::2]
		v0 = ntt2(p0, d//2, m//2)
		v1 = ntt2(p1, d//2, m//2)
		hd = d//ext//2
		v = [0]*(2*hd)
		for i in range(hd):
			w = root_ntt(prim_pow[ext*i], m)
			a, b = v0[i], w * v1[i]
			v[i] = a + b
			v[i+hd] = a - b
		return v
	else:
		raise NotImplementedError()

def intt2(v, d=D, m=M):
	assert len(v)*ext == d == euler_phi(m)
	if d == 1:
		return [int(v[0])]
	elif d == 2:
		if ext == 2:
			v0, v1 = v[0], v[0].conjugate()
		else:
			v0, v1 = v[0], v[1]
		w1, w2 = root_ntt(prim_pow[0],m), root_ntt(prim_pow[1],m)
		p1 = (v0-v1) / (w1 - w2)
		p0 = v[0] - w1 * p1
		return [int(p0), int(p1)]
	elif d % 3 == 0:
		hd = d//ext//3
		v0, v1, v2 = [0]*hd, [0]*hd, [0]*hd
		j = root_ntt(1,3)
		j2 = j*j
		for i in range(hd):
			w = root_ntt(prim_pow[ext*i], m)
			a, b, c = v[i], v[i+hd], v[i+2*hd]
			v0[i] = (a + b + c) / 3
			v1[i] = (a + j2*b + j*c) / (3*w)
			v2[i] = (a + j*b + j2*c) / (3*w*w)
		p0 = intt2(v0, d//3, m//3)
		p1 = intt2(v1, d//3, m//3)
		p2 = intt2(v2, d//3, m//3)
		return [c for t in zip(p0,p1,p2) for c in t]
	elif d % 2 == 0:
		hd = d//ext//2
		v0, v1 = [0]*hd, [0]*hd
		for i in range(hd):
			w = root_ntt(prim_pow[ext*i], m)
			a, b = v[i], v[i+hd]
			v0[i] = (a + b)/2
			v1[i] = (a - b)/2 / w
		p0 = intt2(v0, d//2, m//2)
		p1 = intt2(v1, d//2, m//2)
		return [c for t in zip(p0,p1) for c in t]
	else:
		raise NotImplementedError()

def fqk_to_coefs(x):
	coefs = x.polynomial().list()
	return coefs + [fq(0)]*(2-len(coefs))

def ntt3(p):
	assert len(p) == D
	v = p.copy()
	d = 1
	m = M//D # 2 or 3
	rd = D

	io = D//2 # offset between integer and "imaginary" part, if ext=2
	def read_fqk(o):
		if ext == 2:
			return v[o] + alpha*v[o+D//2]
		else:
			return v[o]
	def write_fqk(o, x):
		if ext == 2:
			v[o], v[o+D//2] = fqk_to_coefs(x)
		else:
			v[o] = x

	while m < M:
		if D % (2*d) == 0:
			d, m, rd = 2*d, 2*m, rd//2
		elif D % (3*d) == 0:
			d, m, rd = 3*d, 3*m, rd//3
		else:
			raise NotImplementedError()

		s = D//d # step
		for o in range(s): # offset
			if d == 2 and ext == 2:
				va = v[o] + root_ntt(1,m)*v[o+s]
				write_fqk(o, va)

			elif m == 6:
				v0, v1 = v[o], v[o+s]
				va = v0 + v1 * root_ntt(1,6)
				vb = v0 + v1 * root_ntt(5,6)
				write_fqk(o, va)
				write_fqk(o+s, vb)

			elif d % 3 == 0:
				j = root_ntt(1,3)
				j2 = j*j
				for i in range(d//ext//3):
					assert digit_rev(i, d//3) == digit_rev(i*rd*3, D//ext)*ext
					w = root_ntt(prim_pow[digit_rev(i, d//3)], m)
					i0 = o + 3*s*i
					v0 = read_fqk(i0)
					v1 = read_fqk(i0+s)
					v2 = read_fqk(i0+2*s)
					a, b, c = v0, w*v1, w*w*v2
					va = a + b + c
					vb = a + j*b + j2*c
					vc = a + j2*b + j*c
					write_fqk(i0, va)
					write_fqk(i0+s, vb)
					write_fqk(i0+2*s, vc)
				
			elif d % 2 == 0:
				for i in range(d//ext//2):
					assert digit_rev(i, d//2) == digit_rev(i*rd*2, D//ext)*ext
					w = root_ntt(prim_pow[digit_rev(i, d//2)], m)
					i0 = o + 2*s*i
					v0 = read_fqk(i0)
					v1 = read_fqk(i0+s)
					a, b = v0, w * v1
					va = a + b
					vb = a - b
					write_fqk(i0, va)
					write_fqk(i0+s, vb)

			else:
				raise NotImplementedError()
			
	return v

def intt3(v):
	assert len(v) == D
	v = v.copy()
	d = D
	m = M
	io = D//2 # offset between integer and "imaginary" part, if ext=2
	def read_fqk(o):
		if ext == 2:
			return v[o] + alpha*v[o+D//2]
		else:
			return v[o]
	def write_fqk(o, x):
		if ext == 2:
			v[o], v[o+D//2] = fqk_to_coefs(x)
		else:
			v[o] = x

	while d > 1:
		s = D//d # step
		for o in range(s): # offset
			if d == 2:
				va = read_fqk(o)
				if ext == 2:
					vb = va.conjugate()
				else:
					vb = read_fqk(o+s)
				w1, w2 = root_ntt(1,m), root_ntt(prim_pow[1],m)
				v1 = (va-vb) / (w1-w2)
				v0 = va - w1*v1
				v[o], v[o+s] = v0, v1
				
			elif d % 3 == 0:
				j = root_ntt(1,3)
				j2 = j*j
				for i in range(d//ext//3):
					w = root_ntt(prim_pow[digit_rev(i, d//3)], m)
					i0 = o + 3*s*i
					va = read_fqk(i0)
					vb = read_fqk(i0+s)
					vc = read_fqk(i0+2*s)
					v0 = (va +    vb +    vc)
					v1 = (va + j2*vb + j *vc) / w
					v2 = (va + j *vb + j2*vc) / (w*w)
					write_fqk(i0,     v0)
					write_fqk(i0+s,   v1)
					write_fqk(i0+2*s, v2)
				
			elif d % 2 == 0:
				for i in range(d//ext//2):
					w = root_ntt(prim_pow[digit_rev(i, d//2)], m)
					i0 = o + 2*s*i
					va = read_fqk(i0)
					vb = read_fqk(i0+s)
					v0 = (va + vb)
					v1 = (va - vb) / w
					write_fqk(i0,   v0)
					write_fqk(i0+s, v1)

			else:
				raise NotImplementedError()
		
		if d % 3 == 0:
			d, m = d//3, m//3
		elif d % 2 == 0:
			d, m = d//2, m//2
		else:
			raise NotImplementedError()

	for i in range(D):
		v[i] = fq(v[i] / (D//2))

	return v


def split_fqk(l):
	if ext == 1:
		return l
	return [y for plane in zip(*(fqk_to_coefs(x) for x in l)) for y in plane]
def mul_fqk(v1, v2):
	io = D//2
	return split_fqk([
		(v1[i]+alpha*v1[i+io]) * (v2[i]+alpha*v2[i+io])
		for i in range(io)
	])


test_fft(fft1, ifft1)
test_fft(fft2, ifft2, fmt=lambda l: [l[2*i] for i in range(D//2)])
test_fft(fft3, ifft3,
	fmt=lambda l: split_complex([l[digit_rev(i, D)] for i in range(D//2)]),
	mul=mul_complex)

test_ntt(ntt2, intt2, fmt=lambda l: [l[ext*i] for i in range(D//ext)])
test_ntt(ntt3, intt3,
	fmt=lambda l: split_fqk([l[digit_rev(i, D)] for i in range(D//ext)]),
	mul=tbt_mul if ext == 1 else mul_fqk)
