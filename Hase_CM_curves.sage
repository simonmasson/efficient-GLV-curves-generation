load('CM_Q-curves.sage')

def get_Qfamily_equation_from_j(d,j0) :
    [jnum,jden] = get_hasegawa_j_inv(d)
    #We transform j0 from number field element to a polynomial ring element
    L.<r> = PolynomialRing(jnum.parent())
    if j0 in QQ :
        newj0 = j0[0]
    else :
        newj0 = j0[0] + r * j0[1]
    return L(newj0*jden - jnum)
    
def compute_Hase_CM_curves(d, T1, T2) :
    table_thm = []
    for t in T1:
        Q = get_Qfamily_equation_from_j(d,t[1][0])
        # Q is in QQ[sqrtD][s][r] but degree 0 (no r), and finally Q is in QQ[s*sqrtD]
        Q = Q[0]
        GCD = gcd(Q[0],Q[1])
        if GCD != 1:
            #GCD is in QQ[s][D]
            if GCD.degree() >1 :
                print 'Problem !', GCD
            if GCD.degree() == 1 :
                if GCD[1].degree()!= 2 or GCD[0].degree()>=1 :
                    print 'ALERT : GCD is not in the form a * s^2*D + b : ', GCD
                else :
                    table_thm.append([t[0],separate_square_factors(-GCD[0][0] / GCD[1][2])])
            else :
                if GCD.degree() == 0 :
                    for s in GCD[0].roots():
                        if s[0] == 0 :
                            table_thm.append([t[0], [s[0],0]])
                        else :
                            table_thm.append([t[0], [s[0],'D']])
    for t in T2:
        Q = get_Qfamily_equation_from_j(d,t[1][0])

        # Normally, no root with r and sqrtD algebraicly independant
        GCD = gcd(gcd(Q[0][0],Q[0][1]),gcd(Q[1][0],Q[1][1]))
        if GCD != 1 :
            print 'there is a root to find !'
            print t[0] ,';',GCD
        # Searching a root with r = sqrtD
        # Looking for rational roots of Q with sqrtD := r
        r = t[1][0].parent().gen()
        R = (Q[0][0](D = r^2) + r*Q[0][1](D= r^2)) + r * (Q[1][0](D = r^2) + r*Q[1][1](D= r^2))
        for root in R.roots():
            if root[0] in QQ :
                table_thm.append([t[0],[abs(root[0]), ZZ(r^2)]])
    return table_thm

tables_s_Delta = []
tables_s_Delta.append(0)
tables_s_Delta.append(0)
tables_s_Delta.append(compute_Hase_CM_curves(2,T1, T2))
tables_s_Delta.append(compute_Hase_CM_curves(3,T1, T2))
tables_s_Delta.append(0)
tables_s_Delta.append(compute_Hase_CM_curves(5,T1, T2))
tables_s_Delta.append(0)
tables_s_Delta.append(compute_Hase_CM_curves(7,T1, T2))
