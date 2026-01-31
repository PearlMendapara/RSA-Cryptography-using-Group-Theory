"""Group-theoretic view of RSA.

RSA lives in the multiplicative group of units
    G = (ℤ/nℤ)* = { a ∈ {1, …, n−1} : gcd(a, n) = 1 },   |G| = φ(n).

Lagrange's theorem says the order of every element divides |G|, so
a^φ(n) = 1 for every unit a — this is Euler's theorem. Because e·d = 1 + kφ(n),
    (m^e)^d = m · (m^φ(n))^k = m   for every unit m,
which is exactly why decryption undoes encryption.

These functions work by enumeration, so they are meant for small n
(demonstrations and tests), not real key sizes.
"""

from __future__ import annotations

from .number_theory import euler_phi, gcd, mod_pow


def units(n: int) -> list[int]:
    """Elements of (ℤ/nℤ)*."""
    return [a for a in range(1, n) if gcd(a, n) == 1]


def element_order(a: int, n: int) -> int:
    """Smallest r > 0 with a^r ≡ 1 (mod n). Requires gcd(a, n) = 1."""
    if gcd(a, n) != 1:
        raise ValueError(f"{a} is not a unit modulo {n}")
    r, x = 1, a % n
    while x != 1:
        x = x * a % n
        r += 1
    return r


def is_cyclic(n: int) -> bool:
    """True if (ℤ/nℤ)* has a generator, i.e. an element of order φ(n)."""
    phi = euler_phi(n)
    return any(element_order(a, n) == phi for a in units(n))


def carmichael_lambda(n: int) -> int:
    """Exponent of the group: the lcm of all element orders.

    λ(n) divides φ(n), and d = e⁻¹ mod λ(n) also works as an RSA private
    exponent (used by PKCS #1).
    """
    from math import lcm
    result = 1
    for a in units(n):
        result = lcm(result, element_order(a, n))
    return result


def verify_euler_theorem(n: int) -> bool:
    """Check a^φ(n) ≡ 1 (mod n) for every unit a."""
    phi = euler_phi(n)
    return all(mod_pow(a, phi, n) == 1 for a in units(n))


def verify_lagrange(n: int) -> bool:
    """Check that every element order divides the group order φ(n)."""
    phi = euler_phi(n)
    return all(phi % element_order(a, n) == 0 for a in units(n))


def rsa_is_permutation(n: int, e: int) -> bool:
    """Check that m ↦ m^e mod n is a bijection on {0, …, n−1}.

    Holds for squarefree n whenever gcd(e, λ(n)) = 1, including messages that
    are not units (a CRT argument covers those), so every message decrypts.
    """
    return len({mod_pow(m, e, n) for m in range(n)}) == n
