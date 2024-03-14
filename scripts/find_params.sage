d = 16

m = 3*d if d % 3 == 0 else 2*d
print(f"d = {d}, m = {m}:")
for k in [1, 2]:
  remaining = 6
  for q in range(3, 32768, 2):
    if q not in Primes() or (q^k-1) % m != 0:
      continue
    if k == 2 and q % (6 if d % 3 == 0 else 4) == 1:
      # inconvenient
      continue
    print(f"k = {k}, q = {q}")
    remaining -= 1
    if remaining == 0: break
