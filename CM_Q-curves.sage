load('auxiliary_functions.sage')

# Discriminant of CM curves with Hilbert polynomial
# of degree 1 or 2

def construct_Discs_D0_f(list_D0, list_f) :
    list_Discs = []
    for i in range(len(list_D0)):
        list_Discs.append(-list_D0[i]*list_f[i]^2)
    return list_Discs

list_D0_h1 = [3,3,3,4,4,7,7,8,11,19,43,67,163]
list_f_h1 = [1,2,3,1,2,1,2,1,1,1,1,1,1]

list_D0_h2 = [3,3,3,4,4,4,7,8,8,11,15,15,20,24,35,40,51,52,88,91,115,123,148,187,232,235,267,403,427]
list_f_h2 = [4,5,7,3,4,5,4,2,3,3,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

list_Disc_h1 = construct_Discs_D0_f(list_D0_h1, list_f_h1)
list_Disc_h2 = construct_Discs_D0_f(list_D0_h2, list_f_h2)

######################################################

# j-invariants of CM curves with Hilbert polynomial
# of degree 1 or 2

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

T1 = construct_CM_j_roots(list_Disc_h1)
T2 = construct_CM_j_roots(list_Disc_h2)
