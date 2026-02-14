"""Why RSA needs large keys: classical attacks and Shor's reduction.

* trial_division / pollard_rho — classical factoring. Both are exponential in
  the bit length, which is why small moduli fall instantly and 2048-bit ones
  do not.
* fermat_factor — breaks keys of any size when p and q are too close together,
  which is why key generation must pick the two primes independently.
* break_rsa — once n is factored, φ(n) and the private exponent d follow
  immediately (report §2.5).
* shor_classical — the classical half of Shor's algorithm (report §3.2–3.3):
  factoring reduces to finding the order r of a random a modulo n. A quantum
  computer finds r in polynomial time with the Quantum Fourier Transform;
  here r is found by brute force, which is exponential, so only small n work.
  Everything else — choosing a, the gcd checks and the post-processing — is
  exactly what Shor's algorithm does classically.
"""

from __future__ import annotations

import random
from dataclasses import dataclass

from .group_theory import element_order
from .number_theory import gcd, is_probable_prime, mod_inverse, mod_pow
from .rsa import PrivateKey, PublicKey


# ------------------------------------------------------------ classical factoring

def trial_division(n: int) -> int | None:
    """Smallest prime factor of n, or None if n is prime. O(√n) divisions."""
    if n % 2 == 0:
        return 2
    f = 3
    while f * f <= n:
        if n % f == 0:
            return f
        f += 2
    return None


def pollard_rho(n: int, seed: int = 0) -> int:
    """A non-trivial factor of composite n by Pollard's rho (Brent's cycle detection).

    Iterates x ↦ x² + c mod n. Modulo the smallest prime factor p the sequence
    cycles after about √p steps (birthday paradox), and gcd(|x − y|, n)
    exposes p. Expected cost is O(n^(1/4)) — far better than trial division,
    but still exponential in the bit length.
    """
    if n % 2 == 0:
        return 2
    if is_probable_prime(n):
        raise ValueError(f"{n} is prime")
    rng = random.Random(seed)
    while True:
        y, c, m = rng.randrange(1, n), rng.randrange(1, n), 128
        g = r = q = 1
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = gcd(q, n)
                k += m
            r *= 2
        if g == n:                      # overshot: back up one step at a time
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = gcd(abs(x - ys), n)
        if g != n:
            return g


def _isqrt(n: int) -> int:
    """Floor of √n by Newton's method on integers."""
    if n < 2:
        return n
    x = 1 << ((n.bit_length() + 1) // 2)           # start above √n
    while True:
        y = (x + n // x) // 2
        if y >= x:
            return x
        x = y


def fermat_factor(n: int, max_steps: int = 1_000_000) -> int:
    """A factor of odd composite n by Fermat's method: find n = a² − b² = (a − b)(a + b).

    Starts at a = ⌈√n⌉ and steps a upward until a² − n is a perfect square.
    If n = pq, the answer is a = (p + q)/2, reached after about
    (√p − √q)² / (2√n) steps — a single step when p and q share their top
    half of bits, however large n is.
    """
    if n % 2 == 0:
        return 2
    a = _isqrt(n)
    if a * a == n:
        return a
    a += 1
    for _ in range(max_steps):
        b2 = a * a - n
        b = _isqrt(b2)
        if b * b == b2:
            return a - b
        a += 1
    raise RuntimeError(f"p and q are not close enough for Fermat's method on {n}")


def break_rsa(public: PublicKey, factor=pollard_rho) -> PrivateKey:
    """Recover the private key from the public key by factoring n."""
    p = factor(public.n)
    q = public.n // p
    phi = (p - 1) * (q - 1)
    d = mod_inverse(public.e, phi)
    return PrivateKey(n=public.n, d=d, p=p, q=q, d_p=d % (p - 1),
                      d_q=d % (q - 1), q_inv=mod_inverse(q, p))


# ------------------------------------------------------------ Shor's reduction

@dataclass
class ShorAttempt:
    a: int
    outcome: str                    # "lucky gcd", "odd order", "trivial root", "success"
    order: int | None = None
    factors: tuple[int, int] | None = None


def shor_step(n: int, a: int) -> ShorAttempt:
    """One run of Shor's classical post-processing for a chosen base a (report §3.3)."""
    g = gcd(a, n)
    if g > 1:                                        # step 2: lucky classical factor
        return ShorAttempt(a, "lucky gcd", factors=(g, n // g))

    r = element_order(a, n)                          # step 4: the quantum part (brute force here)
    if r % 2 == 1:
        return ShorAttempt(a, "odd order", order=r)

    half = mod_pow(a, r // 2, n)
    if half == n - 1:                                # a^(r/2) ≡ −1: only trivial factors
        return ShorAttempt(a, "trivial root", order=r)

    p = gcd(half - 1, n)                             # step 5: classical post-processing
    return ShorAttempt(a, "success", order=r, factors=(p, n // p))


def shor_classical(n: int, seed: int | None = None, max_tries: int = 50) -> tuple[tuple[int, int], list[ShorAttempt]]:
    """Factor n by repeating shor_step with random bases until one succeeds.

    Returns the factors and the log of every attempt. For n = pq with odd
    distinct primes, each random a succeeds with probability at least 1/2.
    """
    if n % 2 == 0:
        return (2, n // 2), []
    rng = random.Random(seed)
    attempts = []
    for _ in range(max_tries):
        attempt = shor_step(n, rng.randrange(2, n - 1))
        attempts.append(attempt)
        if attempt.factors:
            return attempt.factors, attempts
    raise RuntimeError(f"no factor found for {n} in {max_tries} tries")
