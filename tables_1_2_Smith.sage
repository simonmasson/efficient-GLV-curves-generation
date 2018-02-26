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
        
def construct_Discs_D0_f(list_D0, list_f) :
    list_Discs = []
    for i in range(len(list_D0)):
        list_Discs.append(-list_D0[i]*list_f[i]^2)
    return list_Discs
    
def construct_CM_j_roots(list_Discs) :
    T = []
    for D in list_Discs :
        T.append([D,roots_d2(KK(sage.schemes.elliptic_curves.cm.hilbert_class_polynomial(D,algorithm=None)))])
    return T

def print_table_CM_j_roots(table) :
    for l in table :
        D = l[0]
        [f,D0] = separate_square_factors(D)
        if len(l[1]) == 1:
            print D0,'*', f,'^2  ;  ', [(0 if j == 0 else j.factor()) for j in l[1]]
        if len(l[1]) > 1 :
            GCD = gcd([coeff for coeff in l[1][0]])
            j = l[1][0]
            print D0,'*',f,'^2  ;  ', GCD , '(',j[0]/GCD,'+/-', j[1]/GCD,'sqrt(', j.parent().gen()^2 ,'))'
            
list_D0_h1 = [3,3,3,4,4,7,7,8,11,19,43,67,163]
list_f_h1 = [1,2,3,1,2,1,2,1,1,1,1,1,1]

list_D0_h2 = [3,3,3,4,4,4,7,8,8,11,15,15,20,24,35,40,51,52,88,91,115,123,148,187,232,235,267,403,427]
list_f_h2 = [4,5,7,3,4,5,4,2,3,3,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

tt = cputime()

list_Disc_h1 = construct_Discs_D0_f(list_D0_h1, list_f_h1)
list_Disc_h2 = construct_Discs_D0_f(list_D0_h2, list_f_h2)

T1 = construct_CM_j_roots(list_Disc_h1)
T2 = construct_CM_j_roots(list_Disc_h2)

print 'tables 1 and 2 from Smith loaded in ', cputime(tt), ' seconds.'
