"""Tests for the number theory, RSA, group theory and attack modules.

Run with:  pytest -q
"""

import math
import random

import pytest

from rsa_group import attacks, group_theory as gt, number_theory as nt, rsa


# ------------------------------------------------------------ number theory

@pytest.mark.parametrize("seed", range(100))
def test_gcd_and_bezout(seed):
    rng = random.Random(seed)
    a, b = rng.randint(-10**12, 10**12), rng.randint(1, 10**12)
    g, x, y = nt.extended_gcd(a, b)
    assert abs(g) == math.gcd(a, b) == nt.gcd(a, b)
    assert a * x + b * y == g


@pytest.mark.parametrize("seed", range(100))
def test_mod_inverse_and_mod_pow_match_builtin(seed):
    rng = random.Random(seed)
    n = rng.randint(2, 10**18)
    a = rng.randint(1, n - 1)
    if math.gcd(a, n) == 1:
        assert nt.mod_inverse(a, n) == pow(a, -1, n)
    else:
        with pytest.raises(ValueError):
            nt.mod_inverse(a, n)
    e = rng.randint(0, 10**30)
    assert nt.mod_pow(a, e, n) == pow(a, e, n)


def test_miller_rabin_small_numbers():
    sieve = [p for p in range(2, 5000) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    assert [k for k in range(5000) if nt.is_probable_prime(k)] == sieve


def test_miller_rabin_rejects_carmichael_numbers():
    # Fermat's test is fooled by these; Miller–Rabin is not
    for n in (561, 1105, 1729, 2465, 2821, 6601, 8911, 41041, 825265, 321197185):
        assert not nt.is_probable_prime(n)


def test_miller_rabin_known_large_prime():
    assert nt.is_probable_prime(2**127 - 1)          # Mersenne prime
    assert not nt.is_probable_prime(2**128 + 1)       # 2^128 + 1 = 59649589127497217 · 5704689200685129054721


@pytest.mark.parametrize("bits", [8, 16, 64, 256])
def test_random_prime_has_exact_bit_length(bits):
    p = nt.random_prime(bits)
    assert p.bit_length() == bits and nt.is_probable_prime(p)


# ------------------------------------------------------------ RSA

@pytest.fixture(scope="module")
def keypair():
    return rsa.generate_keypair(1024)


def test_key_properties(keypair):
    pub, priv = keypair
    assert pub.bits == 1024
    assert priv.p * priv.q == pub.n
    phi = (priv.p - 1) * (priv.q - 1)
    assert pub.e * priv.d % phi == 1


def test_encrypt_decrypt_roundtrip(keypair):
    pub, priv = keypair
    rng = random.Random(0)
    for _ in range(20):
        m = rng.randrange(pub.n)
        c = rsa.encrypt(m, pub)
        assert rsa.decrypt(c, priv) == m
        assert rsa.decrypt_crt(c, priv) == m


def test_text_roundtrip(keypair):
    pub, priv = keypair
    msg = "Euler's theorem: m^(ed) ≡ m (mod n) 🔐"
    assert rsa.decrypt_text(rsa.encrypt_text(msg, pub), priv) == msg


def test_message_too_long(keypair):
    pub, _ = keypair
    with pytest.raises(ValueError):
        rsa.encrypt_text("x" * 200, pub)


def test_signatures(keypair):
    pub, priv = keypair
    sig = rsa.sign(b"transfer 100 to Alice", priv)
    assert rsa.verify(b"transfer 100 to Alice", sig, pub)
    assert not rsa.verify(b"transfer 900 to Alice", sig, pub)
    assert not rsa.verify(b"transfer 100 to Alice", sig + 1, pub)


def test_textbook_example():
    # Classic worked example: p = 61, q = 53, e = 17 → n = 3233, d = 2753
    pub, priv = rsa.keypair_from_primes(61, 53, e=17)
    assert (pub.n, priv.d) == (3233, 2753)
    assert rsa.encrypt(65, pub) == 2790
    assert rsa.decrypt(2790, priv) == 65


def test_textbook_rsa_is_multiplicatively_malleable(keypair):
    # Enc(a)·Enc(b) = Enc(a·b): why real RSA needs padding
    pub, priv = keypair
    a, b = 12345, 67890
    product = rsa.encrypt(a, pub) * rsa.encrypt(b, pub) % pub.n
    assert rsa.decrypt(product, priv) == a * b


# ------------------------------------------------------------ group theory

@pytest.mark.parametrize("n", [15, 21, 33, 35, 55, 77, 91, 143, 221, 3233])
def test_euler_and_lagrange(n):
    assert len(gt.units(n)) == nt.euler_phi(n)
    assert gt.verify_euler_theorem(n)
    assert gt.verify_lagrange(n)


@pytest.mark.parametrize("n", [15, 21, 33, 35, 55, 77, 91, 143])
def test_carmichael_lambda_divides_phi(n):
    assert nt.euler_phi(n) % gt.carmichael_lambda(n) == 0


def test_cyclic_groups():
    assert gt.is_cyclic(7) and gt.is_cyclic(9) and gt.is_cyclic(18)
    assert not gt.is_cyclic(15) and not gt.is_cyclic(8)   # (ℤ/15)* ≅ C2 × C4


def test_rsa_map_is_a_permutation_even_on_non_units():
    # Covers messages sharing a factor with n, which Euler's theorem alone does not
    assert gt.rsa_is_permutation(3233, 17)
    assert not gt.rsa_is_permutation(3233, 2)              # gcd(2, φ(n)) ≠ 1


# ------------------------------------------------------------ attacks

def test_trial_division():
    assert attacks.trial_division(3233) == 53
    assert attacks.trial_division(101) is None


@pytest.mark.parametrize("bits", [32, 48, 64])
def test_pollard_rho_breaks_small_keys(bits):
    pub, priv = rsa.generate_keypair(bits)
    recovered = attacks.break_rsa(pub)
    assert {recovered.p, recovered.q} == {priv.p, priv.q}
    m = 42
    assert rsa.decrypt(rsa.encrypt(m, pub), recovered) == m


@pytest.mark.parametrize("n", [3233, 5959, 10403, 1022117])
def test_fermat_factor_small(n):
    p = attacks.fermat_factor(n)
    assert 1 < p < n and n % p == 0


def test_fermat_breaks_2048_bit_key_with_close_primes():
    # p and q differ only in their low bits, so a = (p + q)/2 is found at once
    p = nt.random_prime(1024)
    q = p + 2
    while not nt.is_probable_prime(q):
        q += 2
    pub, priv = rsa.keypair_from_primes(p, q)
    recovered = attacks.break_rsa(pub, factor=attacks.fermat_factor)
    assert {recovered.p, recovered.q} == {p, q}
    assert rsa.decrypt(rsa.encrypt(42, pub), recovered) == 42


def test_fermat_gives_up_on_well_separated_primes():
    with pytest.raises(RuntimeError):
        attacks.fermat_factor(1_000_003 * 999_999_937, max_steps=100)


@pytest.mark.parametrize("n", [15, 21, 35, 77, 91, 143, 221, 3233])
def test_shor_reduction_factors(n):
    (p, q), attempts = attacks.shor_classical(n, seed=1)
    assert p * q == n and 1 < p < n
    assert attempts[-1].factors == (p, q)


def test_shor_step_on_the_standard_example():
    # a = 7, n = 15: order 4, 7² = 49 ≡ 4, gcd(3, 15) = 3 and gcd(5, 15) = 5
    step = attacks.shor_step(15, 7)
    assert step.order == 4 and step.outcome == "success"
    assert sorted(step.factors) == [3, 5]
