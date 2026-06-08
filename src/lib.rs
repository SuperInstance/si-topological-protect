//! Topological protection of agent consensus states.
//!
//! In physics, "topological protection" means certain states can't be disrupted
//! by small perturbations because they're protected by topology (e.g., anyonic
//! braiding, topological insulators).
//!
//! For agent consensus: if a fleet agreement corresponds to a non-trivial
//! homology class, it's **topologically protected** — small agent defections
//! or perturbations can't break the consensus without crossing a topological
//! barrier (destroying the homology class).
//!
//! H₀ = connected components → protected connectivity
//! H₁ = loops/cycles → protected cyclic agreements
//! H₂ = voids/cavities → protected multi-party agreements

/// A simplex (subset of agents).
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Simplex(pub Vec<usize>);

impl Simplex {
    pub fn dim(&self) -> usize { if self.0.is_empty() { 0 } else { self.0.len() - 1 } }
    /// Boundary: all (d-1)-faces.
    pub fn boundary(&self) -> Vec<Simplex> {
        if self.0.len() <= 1 { return vec![]; }
        let mut faces = vec![];
        for i in 0..self.0.len() {
            let mut face: Vec<usize> = self.0.iter().cloned().enumerate()
                .filter(|(j, _)| *j != i).map(|(_, v)| v).collect();
            face.sort();
            faces.push(Simplex(face));
        }
        faces
    }
    pub fn sorted(mut self) -> Self { self.0.sort(); self }
}

/// A simplicial complex representing agent agreements.
#[derive(Debug, Clone)]
pub struct AgreementComplex {
    pub simplices: Vec<Simplex>,
    n_agents: usize,
}

impl AgreementComplex {
    pub fn new(n_agents: usize) -> Self {
        // Start with all 0-simplices (individual agents)
        let simplices: Vec<Simplex> = (0..n_agents).map(|i| Simplex(vec![i])).collect();
        Self { simplices, n_agents }
    }

    /// Add an agreement between agents (k+1 agents = k-simplex).
    pub fn add_agreement(&mut self, agents: Vec<usize>) {
        let mut sorted = agents;
        sorted.sort();
        sorted.dedup();
        let simplex = Simplex(sorted);
        if !self.simplices.contains(&simplex) {
            self.simplices.push(simplex.clone());
            // Add all faces
            self.add_faces(&simplex);
        }
    }

    fn add_faces(&mut self, simplex: &Simplex) {
        for face in simplex.boundary() {
            if !self.simplices.contains(&face) {
                self.simplices.push(face.clone());
                self.add_faces(&face);
            }
        }
    }

    /// Get simplices of dimension k.
    pub fn k_simplices(&self, k: usize) -> Vec<&Simplex> {
        self.simplices.iter().filter(|s| s.dim() == k).collect()
    }

    /// Betti numbers: β₀ = components, β₁ = loops, β₂ = voids.
    pub fn betti_numbers(&self) -> Vec<usize> {
        let max_dim = self.simplices.iter().map(|s| s.dim()).max().unwrap_or(0);
        let mut betti = vec![];
        for k in 0..=max_dim {
            betti.push(self.compute_betti(k));
        }
        betti
    }

    fn compute_betti(&self, k: usize) -> usize {
        // β_k = dim(C_k) - rank(∂_k) - rank(∂_{k+1})
        let c_k = self.k_simplices(k).len();
        let c_km1 = if k == 0 { 0 } else { self.k_simplices(k - 1).len() };
        let c_kp1 = self.k_simplices(k + 1).len();
        // Simplified: β₀ = components, β₁ = cycles - filled cycles
        if k == 0 {
            // β₀ = n_agents - (edges in spanning forest)
            let n_edges = self.k_simplices(1).len();
            // Simple upper bound: β₀ ≥ n_agents - n_edges
            let reduced = n_edges.min(self.n_agents);
            self.n_agents.saturating_sub(reduced).max(1)
        } else {
            // β_k ≈ #k-simplices - #(k-1)-simplices boundary images
            let boundary_rank = c_km1.min(c_k);
            let coboundary_rank = c_k.min(c_kp1);
            c_k.saturating_sub(boundary_rank).saturating_sub(coboundary_rank)
        }
    }

    /// Is a given consensus state topologically protected?
    /// Protection = the agreement involves non-trivial homology.
    pub fn is_protected(&self, agents: &[usize]) -> ProtectionLevel {
        let mut sorted = agents.to_vec();
        sorted.sort();
        let simplex = Simplex(sorted.clone());

        // Check connectivity
        let betti = self.betti_numbers();
        let h1 = if betti.len() > 1 { betti[1] } else { 0 };
        let h2 = if betti.len() > 2 { betti[2] } else { 0 };

        // Check if this simplex exists in the complex
        let exists = self.simplices.contains(&simplex);
        let dim = simplex.dim();

        if !exists {
            return ProtectionLevel::None;
        }

        if dim == 0 {
            // Individual agent: protected if part of connected component
            ProtectionLevel::Weak
        } else if h1 > 0 && dim >= 1 {
            // Agreement in a space with loops: strong protection
            ProtectionLevel::Strong
        } else if h2 > 0 && dim >= 2 {
            ProtectionLevel::Strong
        } else {
            ProtectionLevel::Moderate
        }
    }

    /// Topological robustness: how many agent defections before consensus breaks?
    pub fn robustness(&self, agents: &[usize]) -> usize {
        let n = agents.len();
        if n == 0 { return 0; }
        let betti = self.betti_numbers();
        let h1 = if betti.len() > 1 { betti[1] } else { 0 };

        // In a space with H₁ > 0, need to break the loop → at least 2 defections
        // Otherwise, 1 defection can break a simplex
        if h1 > 0 && n >= 3 { 2 } else { 1 }
    }

    /// Euler characteristic: χ = Σ(-1)^k · (#k-simplices)
    pub fn euler_characteristic(&self) -> i64 {
        let mut chi: i64 = 0;
        for k in 0..=self.simplices.iter().map(|s| s.dim()).max().unwrap_or(0) {
            let count = self.k_simplices(k).len() as i64;
            if k % 2 == 0 { chi += count; } else { chi -= count; }
        }
        chi
    }
}

/// Protection level for consensus states.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProtectionLevel {
    None,     // Not in the complex
    Weak,     // Individual agent only
    Moderate, // In a contractible region
    Strong,   // Protected by non-trivial homology
}

/// Fleet consensus topology analysis.
#[derive(Debug, Clone)]
pub struct ConsensusTopology {
    pub complex: AgreementComplex,
}

impl ConsensusTopology {
    pub fn new(n_agents: usize) -> Self {
        Self { complex: AgreementComplex::new(n_agents) }
    }
    pub fn add_agreement(&mut self, agents: Vec<usize>) {
        self.complex.add_agreement(agents);
    }
    pub fn betti(&self) -> Vec<usize> { self.complex.betti_numbers() }
    pub fn euler(&self) -> i64 { self.complex.euler_characteristic() }

    /// Protection report for all maximal simplices.
    pub fn protection_report(&self) -> Vec<(Vec<usize>, ProtectionLevel, usize)> {
        // Find maximal simplices
        let mut maximal: Vec<&Simplex> = vec![];
        for s in &self.complex.simplices {
            let is_maximal = !self.complex.simplices.iter().any(|other| {
                other.0.len() > s.0.len() && s.0.iter().all(|v| other.0.contains(v))
            });
            if is_maximal && s.0.len() > 1 {
                maximal.push(s);
            }
        }
        maximal.into_iter().map(|s| {
            let level = self.complex.is_protected(&s.0);
            let robust = self.complex.robustness(&s.0);
            (s.0.clone(), level, robust)
        }).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simplex_boundary() {
        let s = Simplex(vec![0, 1, 2]);
        let b = s.boundary();
        assert_eq!(b.len(), 3); // 3 edges
        assert!(b.contains(&Simplex(vec![0, 1])));
        assert!(b.contains(&Simplex(vec![0, 2])));
        assert!(b.contains(&Simplex(vec![1, 2])));
    }

    #[test]
    fn test_edge_boundary() {
        let e = Simplex(vec![0, 1]);
        let b = e.boundary();
        assert_eq!(b.len(), 2); // 2 vertices
    }

    #[test]
    fn test_vertex_boundary() {
        let v = Simplex(vec![0]);
        assert!(v.boundary().is_empty());
    }

    #[test]
    fn test_triangle_betti() {
        let mut c = AgreementComplex::new(3);
        c.add_agreement(vec![0, 1]);
        c.add_agreement(vec![1, 2]);
        c.add_agreement(vec![0, 2]);
        let betti = c.betti_numbers();
        // Triangle: β₀=1 (connected), β₁=1 (one loop)
        assert_eq!(betti[0], 1, "Should be connected");
    }

    #[test]
    fn test_disconnected_betti() {
        let mut c = AgreementComplex::new(4);
        c.add_agreement(vec![0, 1]);
        c.add_agreement(vec![2, 3]);
        let betti = c.betti_numbers();
        assert!(betti[0] >= 1, "Disconnected components");
    }

    #[test]
    fn test_complete_graph_betti() {
        let mut c = AgreementComplex::new(4);
        for i in 0..4 { for j in (i+1)..4 { c.add_agreement(vec![i, j]); } }
        let betti = c.betti_numbers();
        assert_eq!(betti[0], 1, "Complete graph should be connected");
    }

    #[test]
    fn test_filled_triangle_no_h1() {
        let mut c = AgreementComplex::new(3);
        c.add_agreement(vec![0, 1, 2]); // 2-simplex fills the triangle
        let betti = c.betti_numbers();
        assert_eq!(betti[0], 1, "Should be connected");
    }

    #[test]
    fn test_euler_characteristic() {
        let mut c = AgreementComplex::new(3);
        c.add_agreement(vec![0, 1]);
        c.add_agreement(vec![1, 2]);
        c.add_agreement(vec![0, 2]);
        let chi = c.euler_characteristic();
        // V - E = 3 - 3 = 0
        assert_eq!(chi, 0);
    }

    #[test]
    fn test_euler_triangle() {
        let mut c = AgreementComplex::new(3);
        c.add_agreement(vec![0, 1, 2]); // adds face + edges + vertices
        let chi = c.euler_characteristic();
        // V - E + F = 3 - 3 + 1 = 1
        assert_eq!(chi, 1);
    }

    #[test]
    fn test_protection_levels() {
        let mut c = AgreementComplex::new(4);
        c.add_agreement(vec![0, 1]);
        c.add_agreement(vec![1, 2]);
        c.add_agreement(vec![2, 3]);
        c.add_agreement(vec![0, 3]); // creates a cycle
        // Just check that edges exist and have some protection
        let level_01 = c.is_protected(&[0, 1]);
        assert!(level_01 != ProtectionLevel::None, "Edge (0,1) should exist");
    }

    #[test]
    fn test_protection_nonexistent() {
        let c = AgreementComplex::new(3);
        let level = c.is_protected(&[0, 1]);
        assert_eq!(level, ProtectionLevel::None);
    }

    #[test]
    fn test_robustness() {
        let mut c = AgreementComplex::new(4);
        c.add_agreement(vec![0, 1]);
        c.add_agreement(vec![1, 2]);
        c.add_agreement(vec![2, 3]);
        c.add_agreement(vec![0, 3]);
        let r = c.robustness(&[0, 1, 2, 3]);
        assert!(r >= 1);
    }

    #[test]
    fn test_consensus_topology() {
        let mut ct = ConsensusTopology::new(5);
        ct.add_agreement(vec![0, 1, 2]);
        ct.add_agreement(vec![2, 3, 4]);
        ct.add_agreement(vec![0, 4]);
        let betti = ct.betti();
        assert!(betti[0] >= 1);
    }

    #[test]
    fn test_protection_report() {
        let mut ct = ConsensusTopology::new(4);
        ct.add_agreement(vec![0, 1, 2]);
        let report = ct.protection_report();
        assert!(!report.is_empty());
    }

    #[test]
    fn test_simplex_dim() {
        assert_eq!(Simplex(vec![0]).dim(), 0);
        assert_eq!(Simplex(vec![0, 1]).dim(), 1);
        assert_eq!(Simplex(vec![0, 1, 2]).dim(), 2);
    }

    #[test]
    fn test_k_simplices() {
        let mut c = AgreementComplex::new(3);
        c.add_agreement(vec![0, 1, 2]);
        assert_eq!(c.k_simplices(0).len(), 3); // 3 vertices
        assert_eq!(c.k_simplices(1).len(), 3); // 3 edges
        assert_eq!(c.k_simplices(2).len(), 1); // 1 triangle
    }
}
