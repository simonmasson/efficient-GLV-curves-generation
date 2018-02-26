def get_hasegawa_j_inv(d) :
    Q0.<s> = PolynomialRing(QQ)
    Q1.<D> = PolynomialRing(Q0)
    Q2.<sqrtD> = PolynomialRing(Q1)
    if d == 5 :
        Q3.<sqrtD> = Q2.quo(sqrtD^2 +1)
    else :
        Q3.<sqrtD> = Q2.quo(sqrtD^2 - D)
    
    if d == 2 :
        C = 9*(1+s*sqrtD)
        sC = 9*(1-s*sqrtD)
        jnum = -12^3*(C-24)^3
        jden = (C^2*sC)
        
    if d == 3 :
        C = 2*(1+s*sqrtD)
        sC = 2*(1-s*sqrtD)
        jnum = 2^8*3^3*(2*C+1)^3
        jden = (C*sC^3)
        
    if d == 5 :
        jnum = -64*(3*(6*s^2+6*s-1)-20*(s^2-s)*sqrtD)^3
        jden = (1+s^2)*(1+s*sqrtD)^4
        
    if d == 7 :
        C = 7*(27+s^2*sqrtD^2)
        sC = C
        jnum = (27+s^2*sqrtD^2)*(85+96*s*sqrtD+ 15*s^2*sqrtD^2)^3
        jden = (1-s^2*sqrtD^2)*(1-s*sqrtD)^6
        
    return [Q3(jnum),Q3(jden)]
    
def get_hasegawa_reduction_coefficients(p, d, s, Delta) :
    proof.arithmetic(false)
    Fp = GF(p)
    K.<x> = PolynomialRing(Fp)
    if p % 4 == 3 :
        #we write sqrt(Delta) with sqrt(-1)
        #Creation of FF_{p^2}
        L.<X> = GF(p^2, modulus=x^2 + 1)
        #Expression of sqrt(Delta) with sqrt(-1)   (noted Dbar)
        S.<y> = PolynomialRing(L)
        Dbar = S(y^2 - Delta).roots()[0][0]
    if p % 4 == 1 :
        if p == 2^255 - 19:
            #Creation of FF_{p^2}
            L.<X> = GF(p^2, name='X', modulus=x^2+2)
            #Expression of sqrt(Delta) with L.gen()   (noted Dbar)
            S.<y> = PolynomialRing(L)
            Dbar = S(y^2 - Delta).roots()[0][0]
        else :
            #we find the element to write FF_{p^2}
            #Creation of FF_{p^2}
            L.<X> = GF(p^2, name='X')
            #Expression of sqrt(Delta) with L.gen()   (noted Dbar)
            S.<y> = PolynomialRing(L)
            Dbar = S(y^2 - Delta).roots()[0][0]
    #tt = cputime()
    s = Fp(QQ(s))
    #s = L(GF(p)(s.numerator())) * L(GF(p)(s.denominator()))^(-1)
    #print 'time=',cputime(tt)
    if d == 2 :
        C = 9*(1+s*Dbar)
        A = 2*(C-24)
        B = -8*(C-16)
    
    if d == 3 :
        C = 2*(1+s*Dbar)
        A = -3*(2*C+1)
        B = C^2 + 10 * C - 2
    
    if d == 5 :
        A = -27*s*(11*s - 2) * (3*(6*s^2+6*s-1) - 20 * s * (s-1) * Dbar)
        B = 54*s^2*(11*s-2)^2 * ((13*s^2+59*s - 9) - 2 * (s-1) * (20*s+9) * Dbar)
    
    if d == 7 :
        C = 7*(27+s^2*Dbar^2)
        A = -3*C*(85+96*s*Dbar+15*s^2*Dbar^2)
        B = 14*C*(9*(3*s^4*Dbar^4 + 130 *s^2*Dbar^2 + 171) + 16*(9*s^2*Dbar^2 + 163)*s*Dbar)

    return [A,B]
    
def get_Qfamily_equation_from_j(d,j0) :
    [jnum,jden] = get_hasegawa_j_inv(d)
    #We transform j0 from number field element to a polynomial ring element
    L.<r> = PolynomialRing(jnum.parent())
    if j0 in QQ :
        newj0 = j0[0]
    else :
        newj0 = j0[0] + r * j0[1]
    return L(newj0*jden - jnum)
    
def compute_table_thm(d, list_Disc_h1, list_Disc_h2) :
    T1 = construct_CM_j_roots(list_Disc_h1)
    T2 = construct_CM_j_roots(list_Disc_h2)
    table_thm = []
    for t in T1:
        Q = get_Qfamily_equation_from_j(d,t[1][0])
        #Q is in QQ[sqrtD][s][r] but degree 0 (no r), and finally Q is in QQ[s*sqrtD]
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

        #Recherche d'une solution avec r et sqrtD algébriquement indépendants (aucune)
        GCD = gcd(gcd(Q[0][0],Q[0][1]),gcd(Q[1][0],Q[1][1]))
        if GCD != 1 :
            print 'there is a root to find !'
            print t[0] ,';',GCD
        #Recherche d'une solution avec r et sqrtD algébriquement dépendants, thm => r = sqrtD
        #On remplace sqrtD par r dans Q puis on cherche les racines rationnelles
        r = t[1][0].parent().gen()
        R = (Q[0][0](D = r^2) + r*Q[0][1](D= r^2)) + r * (Q[1][0](D = r^2) + r*Q[1][1](D= r^2))
        for root in R.roots():
            if root[0] in QQ :
                table_thm.append([t[0],[abs(root[0]), ZZ(r^2)]])
    return table_thm
    
tt = cputime()

tables_s_Delta = []
tables_s_Delta.append(0)
tables_s_Delta.append(0)
tables_s_Delta.append(compute_table_thm(2,list_Disc_h1, list_Disc_h2))
tables_s_Delta.append(compute_table_thm(3,list_Disc_h1, list_Disc_h2))
tables_s_Delta.append(0)
tables_s_Delta.append(compute_table_thm(5,list_Disc_h1, list_Disc_h2))
tables_s_Delta.append(0)
tables_s_Delta.append(compute_table_thm(7,list_Disc_h1, list_Disc_h2))

print 'thm_6_Smith.sage loaded in ', cputime(tt), ' seconds.'
