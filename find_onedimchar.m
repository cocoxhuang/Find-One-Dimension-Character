function find_onedimchar(C,ell : primes_bound := 3000);
    Z := Integers();
    L := LSeries(C);
    P_ell<T> := PolynomialRing(GF(ell)); /*GF(ell) returns the finite field of ell*/
    J := Jacobian(C);
    cond := ell*Conductor(C);
    print "Conductor of curve is: ", Conductor(C);
    bad_primes := [p[1]: p in Factorization(Conductor(C))]; /* Factorization(cond) returns <prime, multiplicity>*/
    print "Primes dividing conductor are: ", bad_primes;
    G := DirichletGroup(cond,GF(ell)); /* returns the group of Dirichlet characters modulo cond with values in GF(ell)*/
    
    /* compute the the possible generators of psi in DirichletGroup of conductor cond in the finite field of ell*/
    gens_G := Generators(G);
    n := #gens_G;
    /* exps_G is a list of [order of generator 1, ..., order of generator n] */
    exps_G := [Order(G.i) : i in [1..n]]; /* G_i returns the i-th generator */

    /* Recall G_i is the i-th generator os G */
    /* Computes the the possible elements: 
    G_1^{x[1]} * G_2^{x[2]} ... * G_n{x[n]} as generators of characters psi, psi^{-1} chi*/
    /* X = {x | x= <x[1], x[2], ... x[n]>} */
    X := Set(CartesianProduct([[0..e-1] : e in exps_G]));
    p := 3;    
    /* Repeat untill we only have two generators: generator for psi and psi^{-1} chi */
    while #X gt 2 do
        if cond mod p ne 0 then
            Jp := BaseExtend(J,GF(p));
            invcharpol := EulerFactor(Jp);   
            charpol := P_ell ! Reverse(Coefficients(invcharpol));
            eigvals_rhoell_frobp := [r[1] : r in Roots(charpol)]; /* roots(charpol) returns <root, multiplicity> */
            X := [x : x in X | Evaluate(&*[(G.i)^(x[i]) : i in [1..n]],p) in eigvals_rhoell_frobp];
        end if;
        p := NextPrime(p);
        if p gt primes_bound then
            printf "Exceeded the bound for primes";
            break;
        end if;
    end while;
    X_chars := [&*[(G.i)^x[i] : i in [1..n]] : x in X];
    return X_chars;
end function;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//     Example: Finding the 1-dim subquotient reps of \ell-torsion representation of a genus 2 Jacobian    //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

P<x> := PolynomialRing(Rationals());
C := HyperellipticCurveOfGenus(2,[(1)*x^5 + (-31)*x^4 + (-110)*x^3 + (21)*x^2 + (-1)*x^1 + (0),(1)*x^1 + (0)]);
onedim_subreps_of5tors := find_onedimchar(C,5);

print "------------------------";
print "(Order might be inaccurate) <Conductor(psi), Order(psi)> and <Conductor(psi^{-1} chi), Order(psi^{-1} chi)> are: ", [<Conductor(chi), Order(chi)> : chi in onedim_subreps_of5tors];
print "(Order might be inaccurate) In the format of [psi(1), psi(2) ,... ], value lists of psi, psi^{-1} chi are: ", [ValueList(AssociatedPrimitiveCharacter(chi)) : chi in onedim_subreps_of5tors];
