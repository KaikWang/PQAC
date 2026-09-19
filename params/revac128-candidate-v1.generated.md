# Generated parameter-profile audit

- Profile: revac128-candidate-v1
- Status: candidate
- Profile ID (SHA3-256): ef5a80fcf061b3a7c0307d8d046239bac641143e7cef3e6b07c468bf12f2f14a
- params.h SHA3-256: a0e8ede1050d291a94cf9b8f2ec816268b300a2786d71a6d267f0a81f91b397e
- Computational target: 128 bits
- Statistical target: 128 bits
- Balanced gadget centered coverage: 265720
- Centered modulus requirement: 202884

## Active implementation parameters

| Component | Frozen candidate values |
|---|---|
| Ring | n=256, d=4, q=405769 (19 bits) |
| Accumulator gadget | balanced base 3, dimension 12, 32-bit handles |
| Update signs | fixed +1 (PKC theorem-compatible) |
| Truncated sampler | b=14, k=5, l=2 |
| Issue proof | n=64, rank=20, q=212704515569 (38 bits), m1=104, m2=58, l=7 |
| Show proof | n=64, rank=23, q=216093359578546177 (58 bits), m1=183, m2=74, l=4 |
| LaBRADOR | q_p=274877906837, grouped right, zk=True |
| Synchronization | public DELTA4, private DELTA6, BML=20 |

## Deterministic structural sizes

- Accumulator value: 2432 bytes (2.38 KiB).
- Compact accumulator witness (without right): 12192 bytes (11.91 KiB).
- Full accumulator witness: 946080 bytes (923.91 KiB).
- Accumulator Ba matrix: 3735552 bytes (3648.00 KiB) in q-bit packed form.
- These are algebraic packed sizes; codec framing, signatures, compression and allocator memory are excluded.

## Release blockers and warnings

- the active 19-bit accumulator fails the PKC'26 beta_1 < q condition for a 500-deletion window
- the screened Issue m2=64 and Show m2=86 performance profile is not compiled; the active 58-bit Show modulus is not the selected approximately 75-bit modulus
- the separate theorem-screened q_acc arithmetic and cross-ring relation are not implemented
- the balanced odd-base and truncated-sampler lemmas are not yet composed with the complete accumulator reduction
- the same-handle commitment binding and hiding parameters are not included in a concrete security reduction
- the Labrador grouped-right norm schedule is an implementation policy rather than a proved accumulator bound
- source parameter scripts still use LOG2_EPS=-40; the composed statistical-distance budget needs rederivation
- compiled Show modulus fits the single-limb backend, while the current ordinary-rejection parameter model selects a larger modulus
- accumulator witness beta remains the permissive prototype bound

## Rule

This records a reproducible candidate, not a 128-bit security claim, until the generator succeeds with --require-release.
