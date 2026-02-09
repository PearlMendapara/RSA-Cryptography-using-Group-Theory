"""Walk through the report's RSA section in code.

    python demo.py             # 2048-bit key
    python demo.py --bits 512  # faster
"""

from __future__ import annotations

import argparse
import time

from rsa_group import attacks, group_theory as gt, number_theory as nt, rsa


def header(title: str) -> None:
    print("\n" + "=" * 64 + f"\n{title}\n" + "=" * 64)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bits", type=int, default=2048, help="RSA modulus size")
    args = parser.parse_args()

    header("1. Worked example with small primes (p = 61, q = 53, e = 17)")
    pub, priv = rsa.keypair_from_primes(61, 53, e=17)
    phi = (61 - 1) * (53 - 1)
    print(f"n = {pub.n}, φ(n) = {phi}, d = e⁻¹ mod φ(n) = {priv.d}")
    print(f"check: e·d mod φ(n) = {17 * priv.d % phi}")
    m = 65
    c = rsa.encrypt(m, pub)
    print(f"encrypt m = {m}: c = {m}^17 mod {pub.n} = {c}")
    print(f"decrypt:        m = {c}^{priv.d} mod {pub.n} = {rsa.decrypt(c, priv)}")

    header("2. The group (ℤ/3233ℤ)* behind it")
    n = pub.n
    print(f"|(ℤ/{n}ℤ)*| = φ({n}) = {len(gt.units(n))}")
    print(f"Carmichael λ({n}) = {gt.carmichael_lambda(n)}  (every element order divides it)")
    print(f"order of 2 = {gt.element_order(2, n)}, order of 65 = {gt.element_order(65, n)}")
    print(f"Euler's theorem holds for all units: {gt.verify_euler_theorem(n)}")
    print(f"m ↦ m^17 permutes all of ℤ/{n}ℤ (units and non-units): {gt.rsa_is_permutation(n, 17)}")

    header(f"3. Real-size key ({args.bits}-bit modulus)")
    t0 = time.perf_counter()
    pub, priv = rsa.generate_keypair(args.bits)
    print(f"key generated in {time.perf_counter() - t0:.2f} s   (e = {pub.e})")
    print(f"n = {str(pub.n)[:40]}…  ({pub.bits} bits)")

    msg = "Harvest now, decrypt later."
    c = rsa.encrypt_text(msg, pub)
    print(f"\nplaintext : {msg}")
    print(f"ciphertext: {str(c)[:40]}…")
    print(f"decrypted : {rsa.decrypt_text(c, priv)}")

    t0 = time.perf_counter(); [rsa.decrypt(c, priv) for _ in range(20)]; plain = time.perf_counter() - t0
    t0 = time.perf_counter(); [rsa.decrypt_crt(c, priv) for _ in range(20)]; crt = time.perf_counter() - t0
    print(f"\ndecryption, 20 runs: plain {plain:.3f} s vs CRT {crt:.3f} s  ({plain / crt:.1f}× faster)")

    sig = rsa.sign(msg.encode(), priv)
    print(f"\nsignature valid: {rsa.verify(msg.encode(), sig, pub)}")
    print(f"tampered message accepted: {rsa.verify(b'Harvest now, decrypt never.', sig, pub)}")

    header("4. Why key size matters: factoring small moduli")
    for bits in (32, 48, 64):
        weak_pub, _ = rsa.generate_keypair(bits)
        t0 = time.perf_counter()
        recovered = attacks.break_rsa(weak_pub)
        print(f"{bits}-bit n = {weak_pub.n}: factored as {recovered.p} × {recovered.q} "
              f"in {time.perf_counter() - t0:.4f} s → private key recovered")

    header("5. Shor's reduction: factoring via order finding")
    n = 3233
    print(f"n = {n}. Try bases a and look at the order r of a modulo n:")
    for a in (2, 13, 3):
        step = attacks.shor_step(n, a)
        half = pow(a, step.order // 2, n) if step.order % 2 == 0 else None
        if step.outcome == "trivial root":
            note = f"a^(r/2) = {half} ≡ −1 (mod n) → only trivial factors, retry"
        elif step.outcome == "odd order":
            note = "r is odd → a^(r/2) undefined, retry"
        else:
            note = (f"a^(r/2) = {half}; gcd({half} − 1, n) = {step.factors[0]}, "
                    f"gcd({half} + 1, n) = {step.factors[1]}")
        print(f"  a = {a:>2}: r = {step.order:>3}  {note}")
    (p, q), _ = attacks.shor_classical(n, seed=3)
    print(f"→ {n} = {p} × {q}")
    print("A quantum computer finds r in polynomial time; here it is brute-forced,")
    print("which is exponential — exactly the gap Shor's algorithm closes.")


if __name__ == "__main__":
    main()
