package graphene.monolayer;

import graphene.lattice.HermitianEigenSolver;
import graphene.lattice.HoneycombLattice;
import graphene.lattice.Vector2;

/**
 * The 2x2 Bloch Hamiltonian for pi-electrons hopping on the honeycomb lattice.
 *
 * <p>Physics reasoning: graphene's electrical conductivity comes from a single
 * electron per carbon atom, left over in the p_z orbital after three of the
 * four valence electrons form the in-plane sigma bonds that hold the lattice
 * together. These p_z orbitals overlap with their neighbors to form the
 * "pi band," and it is this single band of delocalized electrons that we
 * treat here in the tight-binding (linear combination of atomic orbitals)
 * approximation.
 *
 * <p>Because the unit cell has two atoms (A and B), a Bloch state is a
 * 2-component spinor (c_A(k), c_B(k)), one amplitude per sublattice. Hopping
 * an electron from an A site to a neighboring B site with amplitude -t
 * generates an off-diagonal Hamiltonian term; summing over the three nearest
 * neighbors {@code delta_1, delta_2, delta_3} and Fourier transforming gives
 * <pre>
 *     H(k) = [ e0(k)        -t f(k)  ]
 *            [ -t f(k)*      e0(k)   ],      f(k) = sum_i exp(i k . delta_i)
 * </pre>
 * The complex structure factor f(k) carries all of the lattice geometry.
 * Diagonalizing this 2x2 Hermitian matrix gives two bands,
 * E_pm(k) = e0(k) +/- t|f(k)|: the filled "pi" (bonding, valence) band and
 * the empty "pi*" (antibonding, conduction) band. They touch with zero gap
 * exactly where f(k) = 0, which happens only at the two inequivalent
 * Brillouin-zone corners K and K': the Dirac points.
 *
 * <p>Optionally, hopping to the six next-nearest neighbors (same sublattice,
 * amplitude -t') is included. Because next-nearest-neighbor hopping never
 * changes sublattice, it only ever shifts the diagonal e0(k); it does not
 * touch the K/K' band-touching condition (which is set purely by f(k) = 0),
 * but it does break the +E/-E electron-hole symmetry of the bands, since
 * e0(k) is generally different at +k and at points related by inversion. A
 * nonzero t' is the leading correction that makes the conduction and valence
 * bands asymmetric, as observed in real graphene (Castro Neto et al., Rev.
 * Mod. Phys. 81, 109 (2009), Sec. II.A).
 */
public final class MonolayerHamiltonian {

    /** Commonly quoted nearest-neighbor hopping energy for graphene, in eV. */
    public static final double DEFAULT_HOPPING_T_EV = 2.7;

    private final double t;
    private final double tPrime;

    /** Pure nearest-neighbor model (electron-hole symmetric bands). */
    public MonolayerHamiltonian(double hoppingTEv) {
        this(hoppingTEv, 0.0);
    }

    /** Nearest- plus next-nearest-neighbor model (breaks electron-hole symmetry). */
    public MonolayerHamiltonian(double hoppingTEv, double hoppingTPrimeEv) {
        this.t = hoppingTEv;
        this.tPrime = hoppingTPrimeEv;
    }

    public double hoppingT() {
        return t;
    }

    public double hoppingTPrime() {
        return tPrime;
    }

    /**
     * The complex nearest-neighbor structure factor f(k) = sum_i exp(i k . delta_i),
     * returned as {re, im}. Its magnitude, scaled by t, is exactly the energy gap
     * 2t|f(k)| between the two bands at wave vector k; f(k) = 0 defines the Dirac
     * points.
     */
    public double[] structureFactor(Vector2 k) {
        double re = 0.0;
        double im = 0.0;
        for (Vector2 delta : HoneycombLattice.NEAREST_NEIGHBOR_VECTORS) {
            double phase = k.dot(delta);
            re += Math.cos(phase);
            im += Math.sin(phase);
        }
        return new double[] {re, im};
    }

    /** The next-nearest-neighbor structure factor g(k), real-valued by inversion symmetry of the six NNN vectors. */
    private double nnnStructureFactor(Vector2 k) {
        double sum = 0.0;
        for (Vector2 delta : HoneycombLattice.NEXT_NEAREST_NEIGHBOR_VECTORS) {
            sum += Math.cos(k.dot(delta));
        }
        return sum;
    }

    /**
     * Diagonalizes H(k) and returns the two band energies {E_-(k), E_+(k)} in
     * ascending order (valence band first, conduction band second).
     */
    public double[] bandEnergies(Vector2 k) {
        double[] f = structureFactor(k);
        double diagonal = -tPrime * nnnStructureFactor(k);

        double[][] real = {
            {diagonal, -t * f[0]},
            {-t * f[0], diagonal}
        };
        double[][] imag = {
            {0.0, -t * f[1]},
            {t * f[1], 0.0}
        };
        return HermitianEigenSolver.eigenvalues(real, imag);
    }
}
