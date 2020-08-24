# Copyright Nathan Lipiarski 2020
# This program based on the algorithm in https://www-jstor-org.offcampus.lib.washington.edu/stable/2153497

using Primes

setprecision(BigFloat, 2^16)

N = 20000

xs = zeros(BigInt, N)
ys = zeros(BigInt, N)

# First convergent is 3 // 1
xs[1] = 3
ys[1] = 1

# Second convergent is 22 // 7
xs[2] = 22
ys[2] = 7

b = 1 / 100
n = 2

while n < N
	@assert (xs[n-1] * ys[n] - xs[n] * ys[n-1] == (-1) ^ (n + 1)) "This check shouldn't fail"
	t_n = xs[n] // ys[n]
	alpha_prime = abs(cos(t_n)) / (ys[n]^2 * abs(sin(t_n))) - ys[n - 1] / ys[n]
	@assert (alpha_prime > 0) "Need to increase floating point accuracy"
	B = max(b * ys[n]^2, ys[n] + 1)

	while ys[n] < B && n < N
		global n = n + 1
		a_n = floor(alpha_prime)
		xs[n] = a_n * xs[n - 1] + xs[n - 2]
		ys[n] = a_n * ys[n - 1] + ys[n - 2]
		alpha_prime = 1 / (alpha_prime - a_n)
	end
end

for i in 1:N
	# We need approximations of the form 2a / 2b + 1
	if xs[i] % 2 == 1 || ys[i] % 2 == 0
		continue
	end

	a = BigInt(xs[i] // 2)
	if isprime(a)
		println("We found a prime!")
		println(a)
		println("\n")
	end
end
