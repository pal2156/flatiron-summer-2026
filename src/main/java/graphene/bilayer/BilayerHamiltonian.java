package graphene.bilayer;

import graphene.lattice.HermitianEigenSolver;
import graphene.lattice.HoneycombLattice;
import graphene.lattice.Vector2;

/**
 * The 4x4 Bloch Hamiltonian for AB- (Bernal-) stacked bilayer graphene.
 *
 * <p>Physics reasoning: stack a second graphene layer on top of the first,
 * shifted so that one sublattice of the top layer (A2) sits directly above
 * one sublattice of the bottom layer (A1); this is the "Bernal" or "AB"
 * stacking realized in natural graphite. The other two sublattices, B1
 * (bottom) and B2 (top), sit above/below the empty centers of the other
 * layer's hexagons, with no atom directly across from them. This geometric
 * distinction between the "dimer" sites (A1, A2) and the "non-dimer" sites
 * (B1, B2) is the entire story of bilayer graphene's low-energy physics.
 *
 * <p>With the four sites ordered (A1, B1, A2, B2), the minimal two-parameter
 * Hamiltonian is
 * <pre>
 *     H(k) = [ -V         -g0 f(k)   g1        0        ]
 *            [ -g0 f(k)*  -V         0         0        ]
 *            [ g1         0          +V        -g0 f(k) ]
 *            [ 0          0          -g0 f(k)*  +V       ]
 * </pre>
 * where f(k) is exactly the nearest-neighbor structure factor from
 * {@code graphene.monolayer.MonolayerHamiltonian} (each layer is, on its
 * own, an ordinary honeycomb lattice with intralayer hopping g0), g1 is the
 * direct A1-A2 "dimer" hopping (the dominant interlayer coupling), and V is
 * half of an applied interlayer bias: an external field pushes charge onto
 * one layer at the expense of the other, raising one layer's on-site energy
 * by V and lowering the other's by V.
 *
 * <p>Two weaker interlayer hoppings, g3 (B1-B2, responsible for "trigonal
 * warping" of the low-energy bands) and g4 (A1-B2 and A2-B1), are left out
 * entirely: not set to zero as an approximation, but genuinely absent from
 * this model, so that every term present can be tied to one clearly
 * identified physical process (Castro Neto et al., Rev. Mod. Phys. 81, 109
 * (2009), Sec. II.C).
 *
 * <p>At k = K exactly (where f(K) = 0), H block-diagonalizes into the dimer
 * pair (A1, A2), with eigenvalues +/- sqrt(V^2 + g1^2), and the non-dimer
 * pair (B1, B2), with eigenvalues +/- V. At V = 0 this reproduces the
 * textbook result: two bands touching exactly at K (from B1/B2) and two
 * bands split off by +/- g1 (from A1/A2, the dimer states). Away from k = K,
 * virtual hopping through the gapped dimer states mixes back into the B1/B2
 * subspace and gives it a small, quadratic-in-k dispersion: the touching
 * bands are quadratic ("massive"), unlike the linear Dirac cone of a single
 * layer, which is exactly the qualitative difference the constructor-level
 * Javadoc above and the generated figures are meant to illustrate.
 */
public final class BilayerHamiltonian {

    /** Commonly quoted interlayer ("dimer") hopping energy for Bernal bilayer graphene, in eV. */
    public static final double DEFAULT_GAMMA_1_EV = 0.39;

    private final double gamma0;
    private final double gamma1;
    private final double interlayerBiasV;

    /** Unbiased bilayer (V = 0). */
    public BilayerHamiltonian(double gamma0Ev, double gamma1Ev) {
        this(gamma0Ev, gamma1Ev, 0.0);
    }

    /**
     * @param gamma0Ev intralayer nearest-neighbor hopping (same role as t in the monolayer model), eV
     * @param gamma1Ev interlayer dimer (A1-A2) hopping, eV
     * @param interlayerBiasVEv half the applied interlayer potential difference, eV; layer 1 sits at
     *     -interlayerBiasV, layer 2 at +interlayerBiasV
     */
    public BilayerHamiltonian(double gamma0Ev, double gamma1Ev, double interlayerBiasVEv) {
        this.gamma0 = gamma0Ev;
        this.gamma1 = gamma1Ev;
        this.interlayerBiasV = interlayerBiasVEv;
    }

    public double gamma0() {
        return gamma0;
    }

    public double gamma1() {
        return gamma1;
    }

    public double interlayerBiasV() {
        return interlayerBiasV;
    }

    /** The same nearest-neighbor structure factor as monolayer graphene: f(k) = sum_i exp(i k . delta_i). */
    private double[] structureFactor(Vector2 k) {
        double re = 0.0;
        double im = 0.0;
        for (Vector2 delta : HoneycombLattice.NEAREST_NEIGHBOR_VECTORS) {
            double phase = k.dot(delta);
            re += Math.cos(phase);
            im += Math.sin(phase);
        }
        return new double[] {re, im};
    }

    /** Diagonalizes H(k) and returns its four band energies in ascending order. */
    public double[] bandEnergies(Vector2 k) {
        double[] f = structureFactor(k);

        double[][] real = {
            {-interlayerBiasV, -gamma0 * f[0], gamma1, 0.0},
            {-gamma0 * f[0], -interlayerBiasV, 0.0, 0.0},
            {gamma1, 0.0, interlayerBiasV, -gamma0 * f[0]},
            {0.0, 0.0, -gamma0 * f[0], interlayerBiasV}
        };
        double[][] imag = {
            {0.0, -gamma0 * f[1], 0.0, 0.0},
            {gamma0 * f[1], 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, -gamma0 * f[1]},
            {0.0, 0.0, gamma0 * f[1], 0.0}
        };
        return HermitianEigenSolver.eigenvalues(real, imag);
    }

    /**
     * The direct gap between the two low-energy ("non-dimer") bands at wave vector k:
     * the smallest positive band energy minus the largest negative one, i.e. E[2] - E[1]
     * out of the four ascending eigenvalues. Zero when V = 0, since those two bands touch.
     */
    public double directGap(Vector2 k) {
        double[] e = bandEnergies(k);
        return e[2] - e[1];
    }

    /**
     * The true minimum direct gap in a neighborhood of K, found by a radial/angular
     * search rather than simply evaluating at k = K. This matters physically: under a
     * nonzero bias the low-energy bands develop a shallow "Mexican hat" shape, and their
     * closest approach is at a small but nonzero |k - K|, not exactly at K (McCann, Phys.
     * Rev. B 74, 161403(R) (2006); Castro Neto et al., Rev. Mod. Phys. 81, 109 (2009),
     * Sec. II.C). This method reports that true minimum, which is the gap an experiment
     * (or a transport measurement) would actually see.
     */
    public double minimumGapNearDiracPoint(Vector2 diracPoint, double searchRadius, int radialSamples, int angularSamples) {
        double minGap = directGap(diracPoint);
        for (int i = 1; i <= radialSamples; i++) {
            double r = searchRadius * i / radialSamples;
            for (int j = 0; j < angularSamples; j++) {
                double theta = 2.0 * Math.PI * j / angularSamples;
                Vector2 k = diracPoint.plus(new Vector2(r * Math.cos(theta), r * Math.sin(theta)));
                minGap = Math.min(minGap, directGap(k));
            }
        }
        return minGap;
    }
}
