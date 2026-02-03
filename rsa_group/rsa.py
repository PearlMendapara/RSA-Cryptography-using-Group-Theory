"""Textbook RSA: key generation, encryption, decryption and signatures.

Follows Section 2 of the report step by step.

⚠ Educational implementation. Real systems must use padding (OAEP for
encryption, PSS for signatures) and constant-time arithmetic. Textbook RSA is
deterministic and malleable: it leaks when the same message is sent twice and
lets an attacker multiply ciphertexts. Use a vetted library (e.g.
`cryptography`) for anything real.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass

from .number_theory import gcd, mod_inverse, mod_pow, random_prime

DEFAULT_E = 65537          # 2^16 + 1: prime, and only two set bits → fast encryption


@dataclass(frozen=True)
class PublicKey:
    n: int
    e: int

    @property
    def bits(self) -> int:
        return self.n.bit_length()


@dataclass(frozen=True)
class PrivateKey:
    n: int
    d: int
    # CRT parameters, kept so decryption can run ~3–4× faster (see decrypt_crt)
    p: int
    q: int
    d_p: int
    d_q: int
    q_inv: int


def generate_keypair(bits: int = 2048, e: int = DEFAULT_E) -> tuple[PublicKey, PrivateKey]:
    """Report §2.2, steps 1–5.

    1. Choose distinct primes p, q of bits/2 bits each
    2. n = p·q
    3. φ(n) = (p − 1)(q − 1)
    4. Use public exponent e with gcd(e, φ(n)) = 1
    5. d = e⁻¹ mod φ(n)  (extended Euclidean algorithm)
    """
    if bits < 16:
        raise ValueError("modulus must be at least 16 bits")
    half = bits // 2
    while True:
        p = random_prime(half)
        q = random_prime(bits - half)
        if p == q:
            continue
        phi = (p - 1) * (q - 1)
        if gcd(e, phi) == 1:
            break

    n = p * q
    d = mod_inverse(e, phi)
    private = PrivateKey(n=n, d=d, p=p, q=q,
                         d_p=d % (p - 1), d_q=d % (q - 1), q_inv=mod_inverse(q, p))
    return PublicKey(n=n, e=e), private


def keypair_from_primes(p: int, q: int, e: int = DEFAULT_E) -> tuple[PublicKey, PrivateKey]:
    """Build a key pair from chosen primes — handy for small worked examples."""
    phi = (p - 1) * (q - 1)
    if gcd(e, phi) != 1:
        raise ValueError(f"e = {e} is not coprime to φ(n) = {phi}")
    n, d = p * q, mod_inverse(e, phi)
    private = PrivateKey(n=n, d=d, p=p, q=q,
                         d_p=d % (p - 1), d_q=d % (q - 1), q_inv=mod_inverse(q, p))
    return PublicKey(n=n, e=e), private


# ---------------------------------------------------------------- integers

def encrypt(m: int, key: PublicKey) -> int:
    """c = m^e mod n  (report §2.3)."""
    if not 0 <= m < key.n:
        raise ValueError("message must satisfy 0 ≤ m < n")
    return mod_pow(m, key.e, key.n)


def decrypt(c: int, key: PrivateKey) -> int:
    """m = c^d mod n  (report §2.3)."""
    return mod_pow(c, key.d, key.n)


def decrypt_crt(c: int, key: PrivateKey) -> int:
    """Decryption via the Chinese Remainder Theorem.

    Compute m_p = c^(d mod p−1) mod p and m_q = c^(d mod q−1) mod q, then
    recombine with Garner's formula:  m = m_q + q · (q⁻¹ (m_p − m_q) mod p).
    Two exponentiations with half-size moduli and exponents are roughly
    3–4× cheaper than one full-size exponentiation.
    """
    m_p = mod_pow(c, key.d_p, key.p)
    m_q = mod_pow(c, key.d_q, key.q)
    h = key.q_inv * (m_p - m_q) % key.p
    return m_q + h * key.q


# ---------------------------------------------------------------- text

def _to_int(data: bytes) -> int:
    return int.from_bytes(data, "big")


def _to_bytes(x: int) -> bytes:
    return x.to_bytes((x.bit_length() + 7) // 8, "big")


def encrypt_text(message: str, key: PublicKey) -> int:
    """Encode UTF-8 text as one integer and encrypt it (message must be shorter than n)."""
    m = _to_int(message.encode("utf-8"))
    if m >= key.n:
        raise ValueError(f"message too long for a {key.bits}-bit key; "
                         "real systems encrypt a symmetric key instead (hybrid encryption)")
    return encrypt(m, key)


def decrypt_text(c: int, key: PrivateKey) -> str:
    return _to_bytes(decrypt_crt(c, key)).decode("utf-8")


# ---------------------------------------------------------------- signatures

def _hash_to_int(message: bytes, n: int) -> int:
    return _to_int(hashlib.sha256(message).digest()) % n


def sign(message: bytes, key: PrivateKey) -> int:
    """σ = H(m)^d mod n with H = SHA-256  (report §2.4)."""
    return mod_pow(_hash_to_int(message, key.n), key.d, key.n)


def verify(message: bytes, signature: int, key: PublicKey) -> bool:
    """Accept iff σ^e mod n = H(m)  (report §2.4)."""
    return mod_pow(signature, key.e, key.n) == _hash_to_int(message, key.n)
