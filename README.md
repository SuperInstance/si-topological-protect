# si-topological-protect

> **Proof of Concept:** Topological protection of agent consensus — homology classes determine which agreements are structurally immune to small perturbations.

## The Insight

In condensed matter physics, **topological protection** means certain quantum states can't be disrupted by local perturbations — they're protected by the topology of the system (topological insulators, anyonic braiding).

For agent consensus, we model agreements as simplices in a **simplicial complex**:
- **0-simplices** = individual agents
- **1-simplices** = pairwise agreements
- **2-simplices** = three-way agreements
- **k-simplices** = (k+1)-party agreements

The **homology** of this complex tells us:
| Betti Number | Meaning | Protection |
|-------------|---------|-----------|
| β₀ = components | Connected subgroups | Protected connectivity |
| β₁ = loops | Cyclic agreements | Loop-breaking requires ≥2 defections |
| β₂ = voids | Multi-party fills | Void-breaking requires ≥3 defections |

## What This Proves

1. **Triangles (2-simplices) have β₁=0** — the loop is "filled" by the face
2. **Empty triangles have β₁=1** — the loop is unprotected
3. **Euler characteristic**: χ = V - E + F encodes the topology
4. **Cycles need ≥2 defections to break** — topologically robust
5. **Chains need only 1** — topologically fragile

## Usage

```rust
use si_topological_protect::*;

// Build consensus complex for 5 agents
let mut ct = ConsensusTopology::new(5);
ct.add_agreement(vec![0, 1, 2]); // three-way agreement
ct.add_agreement(vec![2, 3, 4]); // another three-way
ct.add_agreement(vec![0, 4]);    // bridge agreement

// Topology
let betti = ct.betti();           // Betti numbers
let euler = ct.euler();           // Euler characteristic
let report = ct.protection_report(); // protection for each agreement

// Check specific agreement
let level = ct.complex.is_protected(&[0, 1, 2]);
let robust = ct.complex.robustness(&[0, 1, 2]);
```

## Modules

- `Simplex` — a k-simplex (subset of agents forming an agreement)
- `AgreementComplex` — simplicial complex with Betti numbers, Euler characteristic
- `ConsensusTopology` — fleet-level analysis with protection reports
- `ProtectionLevel` — None/Weak/Moderate/Strong classification
- `robustness()` — minimum defections to break an agreement
- `euler_characteristic()` — topological invariant

## Connection to Conservation Law

Topological protection IS conservation at the structural level:
- A protected agreement (non-trivial H₁) conserves its structure against perturbations
- Conservation of γ + η = C is protected by the topology of the budget allocation graph
- Breaking a conservation law = crossing a topological barrier = changing Betti numbers
- Fleet stability = high Betti numbers = more protected = more conserved

## Mathematical Background

### Simplicial Homology
For simplicial complex K with chain groups Cₖ:
Hₖ(K) = ker(∂ₖ) / im(∂ₖ₊₁)

where ∂ₖ is the boundary operator.

### Betti Numbers
βₖ = dim(Hₖ) = dim(ker(∂ₖ)) - dim(im(∂ₖ₊₁))

### Euler Characteristic
χ(K) = Σₖ(-1)ᵏ · |{k-simplices}|

For connected planar graphs: χ = V - E + F.

### Topological Protection
An agreement is protected if it lives in a non-trivial homology class.
Perturbations that don't change the homology can't break it.

## Tests: 16

Covers: simplex boundaries, Betti numbers for triangles/disconnected/complete/filled graphs, Euler characteristic, protection levels, robustness, consensus topology, protection reports, dimension counting.

## License

MIT
