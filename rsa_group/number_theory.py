"""Number-theory primitives that RSA is built on.

Everything here is implemented from first principles (no pow(x, -1, n),
no math.gcd for the core algorithms) so each step of the report can be
traced in code.
"""

from __future__ import annotations

import secrets

# Small primes used to reject most composites before Miller–Rabin
SMALL_PRIMES = [p for p in range(3, 1000) if all(p % d for d in range(2, int(p ** 0.5) + 1))]


def gcd(a: int, b: int) -> int:
    """Greatest common divisor by the Euclidean algorithm: gcd(a, b) = gcd(b, a mod b)."""
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def extended_gcd(a: int, b: int) -> tuple[int, int, int]:
    """Return (g, x, y) with a·x + b·y = g = gcd(a, b)  (Bézout's identity).

    Iterative version: keeps the invariants
        old_r = a·old_s + b·old_t   and   r = a·s + b·t
    while running the Euclidean algorithm on (old_r, r).
    """
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_r, old_s, old_t


def mod_inverse(a: int, n: int) -> int:
    """Multiplicative inverse of a in (ℤ/nℤ)*, i.e. the x with a·x ≡ 1 (mod n)."""
    g, x, _ = extended_gcd(a % n, n)
    if g != 1:
        raise ValueError(f"{a} has no inverse modulo {n} (gcd = {g})")
    return x % n


def mod_pow(base: int, exponent: int, modulus: int) -> int:
    """base^exponent mod modulus by right-to-left square-and-multiply.

    Uses O(log exponent) multiplications, which is what makes RSA practical:
    a 2048-bit exponent needs about 3,000 modular multiplications, not 2^2048.
    """
    if modulus == 1:
        return 0
    if exponent < 0:
        raise ValueError("negative exponent; use mod_inverse first")
    result = 1
    base %= modulus
    while exponent:
        if exponent & 1:
            result = result * base % modulus
        base = base * base % modulus
        exponent >>= 1
    return result


def is_probable_prime(n: int, rounds: int = 40) -> bool:
    """Miller–Rabin primality test.

    Write n − 1 = 2^s · d with d odd. For a random witness a, n passes if
    a^d ≡ 1 or a^(2^j · d) ≡ −1 (mod n) for some 0 ≤ j < s. A composite n
    passes one round with probability at most 1/4, so 40 rounds give an
    error probability below 2^-80.
    """
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0:
        return False
    for p in SMALL_PRIMES:
        if n == p:
            return True
        if n % p == 0:
            return False

    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1

    for _ in range(rounds):
        a = 2 + secrets.randbelow(n - 3)          # a ∈ [2, n − 2]
        x = mod_pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False                          # a is a witness: n is composite
    return True


def random_prime(bits: int) -> int:
    """Random prime of exactly `bits` bits.

    The top two bits are set so that the product of two such primes has
    exactly 2·bits bits; the low bit is set so the candidate is odd.
    """
    if bits < 3:
        raise ValueError("need at least 3 bits")
    while True:
        candidate = secrets.randbits(bits) | (0b11 << (bits - 2)) | 1
        if is_probable_prime(candidate):
            return candidate


def euler_phi(n: int) -> int:
    """Euler's totient by trial factorisation — for small n only (demos and tests)."""
    result, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            result -= result // p
        p += 1
    if m > 1:
        result -= result // m
    return result
