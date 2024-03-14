load("scripts/test_algo.sage")

polyq.<x> = CyclotomicField(M)
polyz = polyq.ring_of_integers()
polyzq = polyz.quo(polyz.ideal(q), 'n')

def gen_small_poly():
	return polyz([randrange(-1, 2) for _ in range(D)])

def gen_small_inv_poly():
	f = gen_small_poly()
	while True:
		try:
			1/polyzq(f)
		except ZeroDivisionError:
			f = gen_small_poly()
		else:
			break
	return f

# in the following:
# - L is the m-th cyclotomic field
# - K is the m/2-th cyclotomic field
# - Lx is the quadratic extension of K isomorphic to L

# F: member of Lx
# returns member of L isomorphic to F
def from_ext(F, L, d):
	Fl = F.list()
	return L([x for pair in zip(*Fl) for x in pair])

# f: member of L
# returns product of conjugates + norm of f in L/K
def norm(f, K, Lx, m, ext):
	fl = f.list()
	if ext == 2:
		f0, f1 = K(fl[0::2]), K(fl[1::2])
	elif ext == 3:
		f0, f1, f2 = K(fl[0::3]), K(fl[1::3]), K(fl[2::3])
	x = Lx.gen()
	if m == 6:
		fx, fc = f0 + f1*x, f0 + f1 - f1*x
	elif ext == 2:
		fx, fc = f0 + f1*x, f0 - f1*x
	elif ext == 3:
		j = K.zeta(3)
		x2, j2 = x*x, j*j
		fx = f0 + f1*x + f2*x2
		fc = (f0 + f1*j*x + f2*j2*x2)*(f0 + f1*j2*x + f2*j*x2)
	fn = (fx * fc).list()
	for i in range(1, ext):
		assert fn[i] == 0
	return fc, fn[0]

def adjoint(f, d, m, L):
	fv = fft2(f.list(), d, m)
	fv = [z.conjugate() for z in fv]
	fal = ifft2(fv, d, m)
	return L(fal)

def reduce(f, g, F, G, d, m, L):
	fa, ga = adjoint(f, d, m, L), adjoint(g, d, m, L)
	while True:
		k = (F*fa+G*ga) / (f*fa + g*ga)
		k = L([round(c) for c in k.list()])
		if k == 0: break
		F, G = F - k*f, G - k*g
	return F, G

class BadGcd(Exception):
	pass

def tower_solver(f, g, d, m):
	print(f"d = {d}, m = {m}")
	if d == 1:
		f, g = int(f), int(g)
		g,u,v = xgcd(f, -g)
		if q % g != 0:
			raise BadGcd()
		k = q // g
		return v*k, u*k
	
	else:
		if m == 6:
			ext = 2
			d2, m2 = d//2, m//3
		else:
			ext = 3 if d % 3 == 0 else 2
			d2, m2 = d//ext, m//ext
		L.<x> = CyclotomicField(m)
		f, g = L(f), L(g)
		
		K.<y> = CyclotomicField(m2)
		t = polygen(K, 't')
		if m == 6:
			Lx.<x> = K.extension(t^ext - t - y)
		else:
			Lx.<x> = K.extension(t^ext - y) # isomorphic to CyclotomicField(m)
		fc, fn = norm(f, K, Lx, m, ext)
		gc, gn = norm(g, K, Lx, m, ext)

		F2, G2 = tower_solver(fn, gn, d2, m2)
		
		F = from_ext(gc * Lx(F2), L, d)
		G = from_ext(fc * Lx(G2), L, d)
		
		return reduce(f, g, F, G, d, m, L)

def gen_ntru():
	while True:
		f = gen_small_inv_poly()
		g = gen_small_inv_poly()

		try:
			F, G = tower_solver(f, g, D, M)
		except BadGcd:
			print("q not divisible by gcd(N(f), N(g)) -> abort")
			continue
		
		r = f*G - g*F
		if r != q:
			print("(f,g,F,G) not a valid solution -> abort")
			continue
		print("success!")

		print("f =", f)
		print("g =", g)
		print("F =", F)
		print("G =", G)

		f, g = polyzq(f), polyzq(g)
		h = polyzq(f / g)
		print("h =", h)
		assert g*h == f
		return f, g, F, G, h

gen_ntru()