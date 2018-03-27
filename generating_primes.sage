def IntegerToSequence(n,basis):
    L = []
    m = n
    while m > 0 :
        L.append(m%basis)
        m = m // basis
    return L
    
def construct_prime_with_bits(n, size_word) :
    #construct a prime number p such that p = 2^n +/- 2^(n-32) +/- ...
    Res = []
    Ncoeff = n // size_word
    if n % size_word != 0 :
        Ncoeff = Ncoeff + 1
    for i in range(1,3^Ncoeff) :
        p = 2^n
        L = IntegerToSequence(i,3)
        for j in range(len(L)) :
            p = p + 2^(size_word * j) * (L[j]-1)
        if p.is_prime() and not(p in Res) :
            Res.append(p)
        f = L[len(L) - 1] # f vaut 1 ou 2
        if f == 2 :
            p2 = p - 2*2^((len(L)-1) * size_word)
        else :
            p2 = p - 2^((len(L)-1) * size_word)
        if p2.is_prime() and not(p2 in Res) :
            Res.append(p2)
    return Res
    
tt = cputime()

u = cputime()
P64 =  construct_prime_with_bits(256,64)
print cputime(u)
u = cputime()
P63 =  construct_prime_with_bits(252,63)
print cputime(u)
u = cputime()
P51 =  construct_prime_with_bits(255,51)
print cputime(u)
u = cputime()
P52 =  construct_prime_with_bits(260,52)
print cputime(u)
u = cputime()
P53 =  construct_prime_with_bits(265,53)
print cputime(u)
u = cputime()
P32 =  construct_prime_with_bits(256,32)
print cputime(u)
u = cputime()
P28 =  construct_prime_with_bits(252,28)
print cputime(u)
u = cputime()
P26 =  construct_prime_with_bits(260,26)
print cputime(u)

print 'primes with bits built in ', cputime(tt), 'seconds.'

primes = P64 + P63 + P51 + P52 + P53 + P32 + P28 + P26

def construct_prime_epsilon(n,window_e):
    Res = []
    for e in range(window_e) :
       if (2^n + e).is_prime() :
           Res.append([2^n + e,+e])
       if (2^n - e).is_prime() :
           Res.append([2^n - e,-e])
    return Res
    
L = []
for k in range(0,9) :
	tt = cputime()
	L = L + construct_prime_epsilon(256 + k ,2^12) + construct_prime_epsilon(256 - k ,2^12)
	print cputime(tt), 'for k =', k

for [p,e] in L :
    primes.append(p)

primes.append(6^98 -7)
primes.append(8^91 +5)
primes.append(7^98 -2)
primes.append(9^87 +4)
primes.append(8^95 -9)
primes.append(9^99 +4)
primes.append(2^226 - 5)
primes.append(2^228+3)
primes.append(2^233-3)
primes.append(2^235-15)
primes.append(2^243-9)
primes.append(2^266-3)
primes.append(2^273+5)
primes.append(2^285-9)
primes.append(2^291-19)
primes.append(2^292+13)
primes.append(2^295+9)
primes.append(2^301+27)
primes.append(2^308+27)
primes.append(2^310+15)
primes.append(2^317+9)
primes.append(2^319+9)
primes.append(2^320+27)
primes.append(2^321-9)
primes.append(2^327+9)
primes.append(2^328+15)
primes.append(2^336-3)
primes.append(2^341+5)
    
#
#
# the list 'primes' contains the primes
#
#

f=open("primes.sage", "w")

l = len(primes)
i = 1

f.write('primes = [')
f.write(str(primes[0]))

while i < l :
	f.write(',')
	f.write(str(primes[i]))
	i = i + 1
	
f.write(']')
f.close()
