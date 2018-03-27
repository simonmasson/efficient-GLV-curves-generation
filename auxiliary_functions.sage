KK.<x> = QQ[]
   
def separate_square_factors(r) :
    p = r.numerator()
    q = r.denominator()
    s = 1
    if p*q>0 :
        D = 1
    else :
        D = -1
    for l in p.factor() :
        if l[1] >1:
            s = s * l[0]^(l[1]//2)
            D = D * l[0]^(l[1]%2)
        else :
            D = D * l[0]
    for m in q.factor() :
        if m[1]>1:
            s = s/(m[0]^(m[1]//2))
            D = D*(m[0]^(m[1]%2))
            s = s/(m[0]^(m[1]%2))
        else :
            s = s / m[0]
            D = D * m[0]
    return [abs(s),D]
    
def roots_d2(P) :
    if P.degree() >2 :
        print "Polynomial of degree ", P.degree(), ", too difficult to find roots"
    if P.degree() == -1 :
        print "Zero polynomial, every element is a root"
    if P.degree() == 0 :
        print "Non-zero constant polynomial, no root"
    if P.degree() == 1 :
        return [-P[0]/P[1]]
    if P.degree() == 2 :
        a = P[2]
        b = P[1]
        Delta = P.discriminant()
        [s,D] = separate_square_factors(Delta)
        K = NumberField(x^2 - D, 'r')
        return [(-b+s*K.gen())/(2*a),(-b-s*K.gen())/(2*a)]


#
# j-invariants and coefficients of Hasegawa curves
#

def get_hasegawa_j_inv(d) :
    Q0.<s> = PolynomialRing(QQ)
    Q1.<D> = PolynomialRing(Q0)
    Q2.<sqrtD> = PolynomialRing(Q1)
    if d == 5 :
        Q3.<sqrtD> = Q2.quo(sqrtD^2 +1)
    else :
        Q3.<sqrtD> = Q2.quo(sqrtD^2 - D)
    
    # formulas given by Benjamin Smith.
    # can be replaced by 1728 * 4*A^3/(4*A^3 + 27*B^2)
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
    
def DefineReductionQCurve(p, d,s,Delta) :
    if Delta == 'D':#d == 5 :
        Delta = -1
    [A,B] = get_hasegawa_reduction_coefficients(p, d, s, Delta)
    if 4*A^3 + 27 *B^2 != 0 :
        return [EllipticCurve(A.parent(), [A,B]), A, B]
    return ['singular curve', A, B]
