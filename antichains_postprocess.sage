import random
import sys

var('λ d n')

def test_same(p1, p2, odd_n):
    for λ0 in range(1, 7):
        for j in range(1, 5):
            λ0 = random.randint(1, 6)
            n0 = 2*j + odd_n
            if p1.subs(λ == λ0, n == n0) != p2.subs(λ == λ0, n == n0):
                print(f"{p1} and {p2} differ at λ={λ0}, n={n0}",
                      "p1", p1.subs(λ == λ0, n == n0),
                      "p2", p2.subs(λ == λ0, n == n0),
                      file=sys.stderr)
                return False
    return True


def clean(p):
    denom = lcm(x.denominator() for x in
                p.simplify_rational().polynomial(QQ['λ','n']).coefficients())
    if denom < 0:
        denom = -denom

    equiv = (denom * p).simplify_rational()
    equiv = equiv.polynomial(ring=PolynomialRing(QQ['λ'], 'n'))
    return equiv, denom


# this function returns lower_layer + upper_layer, but cleaned up (simplified,
# grouped by d, multiplied by a common expression
def simplify_cluster_odd(lower_layer, upper_layer, j):
    assert(bool((n+3) * binomial(n, (n-3)/2) == binomial(n, (n+1)/2) * (n-1)))
    obv = (lower_layer + upper_layer).subs(
        d == (n-1)/2,
        binomial(n, (n-3)/2) == binomial(n, (n+1)/2) * (n-1)/(n+3))

    poly1 = lower_layer * (1+λ)^(j*(d+2)) * (2*d) / binomial(2*d+1, d-1) / λ^j
    poly2 = upper_layer * (1+λ)^(j*(d+2)) * (2*d+4) / binomial(2*d+1, d+1) / λ^j
    equiv, denom = clean((poly1 + poly2).subs(d == (n-1)/2))
    multiplier = binomial(n, (n+1)/2) * λ^j / ((n+3) * (1+λ)^(j*((n+3)/2)) * denom)

    # equiv * multiplier and obv should be the same. just to be sure, test it
    # with a few random cases.
    assert(test_same(equiv * multiplier, obv, 1))
    return equiv, multiplier


def simplify_cluster_even(lower_layer, upper_layer, j):
    assert(bool(lower_layer == upper_layer))

    poly = 2 * lower_layer*(1+λ)^(j*(d+1)) / λ^j / binomial(2*d, d-1)
    equiv, denom = clean(poly.subs(d == n/2))
    multiplier = binomial(n, n/2 - 1) * λ^j / ((1+λ)^(j*(n/2 + 1)) * denom)

    assert(test_same(equiv * multiplier,
                     (lower_layer + upper_layer).subs(d == n/2),
                     0))
    return equiv, multiplier


while True:
    try:
        j, odd = map(int, input().split())
        odd_string = "odd" if odd else "even"
        lower, upper = map(lambda x: sage_eval(x, locals={'λ':λ, 'd': d, 'n': n}),
                           [input(), input()])

        if odd:
            equiv, multiplier = simplify_cluster_odd(lower, upper, j)
        else:
            equiv, multiplier = simplify_cluster_even(lower, upper, j)
        print(f"\nL{j} for {odd_string} n: {multiplier} * ({equiv})")
        print(f"\nL{j} for {odd_string} n at λ = 1:",
              (multiplier*equiv).subs(λ==1).simplify())
    except EOFError:
        break
