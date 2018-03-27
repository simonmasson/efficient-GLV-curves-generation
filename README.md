# Efficient GLV curves generation
Generation of efficient four-dimensional GLV curves with high security (256-bit)

## How to use ?
`sage : load('generating_curves.sage')`<br>
NB : if you want to use other primes, just change the line 102 of `generating_curves.sage` with `primes = [my, new, primes]`.

## Files description
 - <b>auxiliary_functions.sage</b><br> 
 `separate_square_factors(r)` returns s and D such that r = s^2*D with s rational and D square-free.<br>
 `roots_d2(P)` returns the roots of a degree 2 polynomial.<br>
 `get_hasegawa_j_inv(d)` returns the j-invariant of a Hasegawa Q-curve of degree d.<br>
 `get_hasegawa_reduction_coefficients(p, d, s, D)` returns A and B in GF(p^2) defining the Hasegawa
 curve y^2 = x^3 + A*x + B parametrized by s, d and Delta.

 - <b>CM_Q-curves.sage</b><br>
 `contruct_Discs_D0_f` return the list of possible discriminants for a CM Q-curve with deg(H_D) = 1 or 2.<br>
 `construct_CM_j_roots` returns the possible j-invariants for a CM Q-curve with deg(H_D) = 1 or 2.
 
 - <b>Hase_CM_curves.sage</b><br>
  `get_Qfamily_equation_from_j(d, j0)` returns the equation satisfied by (s, Delta) to get j(d, s, Delta) = j0.<br>
  `compute_Hase_CM_curves.sage` returns the list of (d, s, Delta) for which a CM Q-curve of degree d arises.
 
 - <b>generating_primes.sage</b><br> 
 Computes a list of 256-bit primes with efficient finite field arithmetic, and stores it in a file `primes.sage`.
 
 - <b>generating_curves.sage</b><br>
 `get_possible_traces(p, Delta)` returns the list of the twisted curves traces for a given quadratic finite field F_{p^2} and a given endomorphisms algebra <b>Q</b>(√Δ).<br>
 `get_curve_order_and_twist_orders(P, Q1, traces)` (Q1 = [p]P) return [a, L] where a is the order of the curve, and L is the list of the twisted curves orders.

write_curve(cpt, file, p, A, B, order, N, twist_orders, E, d, s, D) :
 Computes a list of curves with efficient arithmetic, for which we can apply the four-dimensional GLV method. The list is stored in
 `GLV4_curves.txt`. 

- <b>primes.sage</b><br>
File created by `generating_primes.sage`.

- <b>GLV4_curves.txt</b><br>
File created by `generating_curves.sage`. Each element of the list represents a curve, given in the form : `[p, A, B, order, N, [twists], Delta, d, s, D]` where the curve is given by E_{d, s, D} : y^2 = x^3 + A*x + B with A,B in F_{p^2} = F_p(sqrt(Delta)), #E = order, N|#E is prime, and [twists] are the orders of the twists.
