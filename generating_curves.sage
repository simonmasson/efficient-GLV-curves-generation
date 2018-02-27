load('Hase_CM_curves.sage')

from sage.rings.factorint import factor_trial_division
    
def write_possible_curves(tables_s_Delta, file, p) :
    proof.arithmetic(false)
    Fp = GF(p)
    for d in [2,3,5,7] :
        for t in tables_s_Delta[d] :
            s = t[1][0]
            D = t[1][1]
            Delta = t[0]
            if D == 'D' :
                D = -1
            if legendre_symbol(D,p) == -1 :
                [E,A,B] = DefineReductionQCurve(p,  d, s, D)
                if E != 'singular curve' and ((A not in GF(p)) and (B not in GF(p))) :
                    P = E.random_element()
                    while P == 0 :
                        P = E.random_element()
                    #
                    # Q1 and Q2 are used to know
                    # the supersingularity and to
                    # check the traces of the curves
                    #
                    Q1 = p*P
                    if Q1[0] != P[0] :  #(p+1)*P != 0 and (p-1)*P != 0 :
                        K.<a> = NumberField(x^2 - Delta)
                        #we factorize in the maximal order (implicitly)
                        i1 = K.factor(p)[0][0]
                        # Maybe it is the frobenius of one of the twists
                        frob1 = (i1^2).gens_reduced()[0]
                        traces = []
                        traces.append(frob1.trace())
                        frob2 = -frob1
                        traces.append(frob2.trace())
                        [cond,D0] = separate_square_factors(Delta)
                        if D0 == -3 :
                            #the curve is isogenous to a j=0 curve, it has six twists
                            junity = (-1 + a/cond)/2
                            frob3 = frob1 * junity
                            frob4 = frob2 * junity
                            frob5 = frob3 * junity
                            frob6 = frob4 * junity
                            traces.append(frob3.trace())
                            traces.append(frob4.trace())
                            traces.append(frob5.trace())
                            traces.append(frob6.trace())
                        if D0 == -1 :
                            #the curve is isogenous to a j=1728 curve, it has four twists
                            iunity = a/cond
                            frob3 = frob1 * iunity
                            frob4 = frob2 * iunity
                            traces.append(frob3.trace())
                            traces.append(frob4.trace())
                        twist_orders = []
                        detect_pb = true
                        Q2 = p*Q1 + P # Q2 = (p^2 + 1)*P
                        for t in traces :
                            #we check if N*P == 0
                            if int(t) * P == Q2 :
                                order = int(p^2 + 1 - t)
                                detect_pb = not(detect_pb)
                            else :
                                twist_orders.append(int(p^2+1 - t))
                        if detect_pb :
                            print 'there is a problem with the trace computation'
                        trial = factor_trial_division(order,2^9)
                        boolisprime = trial[len(trial) - 1][0].is_prime()
                        if boolisprime :
                            print 'p =', p
                            print 'subgroup of size ', ceil(log(trial[len(trial) - 1][0])/log(2))
                            f.write('[')
                            f.write(str(p))
                            f.write(',')
                            f.write(str(A))
                            f.write(',')
                            f.write(str(B))
                            f.write(',')
                            f.write(str(order))
                            f.write(',')
                            f.write(str(trial[len(trial) - 1][0]))
                            f.write(',[')
                            f.write(str(twist_orders[0]))
                            for i in range(1,len(twist_orders)-1) :
                                f.write(',')
                                f.write(str(twist_orders[i]))
                            f.write('],')
                            f.write(str(E.base_field().gen()))
                            f.write('^2+')
                            f.write(str(-(E.base_field().gen())^2))
                            f.write(',')
                            f.write(str(d))
                            f.write(',')
                            f.write(str(s))
                            f.write(',')
                            f.write(str(D))
                            f.write('],\n')

#####################################
# if primes.sage is computed :
load('primes.sage')

# else :
# load('generating_primes.sage')
#####################################

f=open("GLV4_curves.txt", "w")
f.write('curves = [   ')
i = 1
length = len(primes)
for p in primes :
	tt = cputime()
	write_possible_curves(tables_s_Delta, f, p)
	print 'computed in ', cputime(tt), 's'
	print i,'/', length, '\n'
	i = i+1
	
f.write('   ]')
f.close()

print ' Warining ! Remove the last comma in the file GLV4_curves.txt !'
