# Efficient GLV curves generation
Generation of efficient four-dimensional GLV curves with high security (256-bit)

Included :<br>
 - <b>tables_1_2_Smith.sage</b><br> 
 Computes the possible j-invariants for a fixed discriminant D_0f^2 for which a CM Hasegawa QQ-curve arises.
 
 - <b>thm_6_Smith.sage</b><br>
 Computes (s,\Delta) for which E_{d,s,\Delta} is a degree d QQ-curve for which the four-dimensional GLV method is possible.
 
 - <b>generating_primes.sage</b><br> 
 Computes a list of primes with efficient finite field arithmetic, and stores it in a file `primes.sage`.
 
 - <b>generating_curves.sage</b><br>
 Computes a list of curves with efficient arithmetic, for which we can apply the four-dimensional GLV method. The list is stored in
 `GLV4_curves.txt`.

- <b>primes.sage</b><br>
File created by `generating_primes.sage`.

- <b>GLV4_curves.txt</b><br>
File created by `generating_curves.sage`.
