package graphene.lattice;

/**
 * Geometry of the graphene honeycomb lattice.
 *
 * <p>Graphene's carbon atoms sit on a honeycomb lattice, which is <em>not</em> itself
 * a Bravais lattice: there is no single set of lattice vectors that reaches every
 * atom from a given atom. Instead, the honeycomb is two interpenetrating triangular
 * Bravais sublattices, conventionally labeled A and B. Every A atom has three nearest
 * neighbors, all of which are B atoms, and vice versa; this A/B bipartite structure
 * is what ultimately produces the two-band (2x2 Hamiltonian) physics of monolayer
 * graphene and the two independent Dirac cones at the K and K' points.
 *
 * <p>Basis: one A atom and one B atom per unit cell, with the A atom taken as the
 * origin and the B atom displaced by {@link #DELTA_1} (one of the three
 * nearest-neighbor vectors).
 *
 * <p>Sign/axis convention follows Castro Neto, Guinea, Peres, Novoselov &amp; Geim,
 * "The electronic properties of graphene," Rev. Mod. Phys. 81, 109 (2009), Sec. II.A
 * and Fig. 2, which we also use as the numerical ground truth for the bands computed
 * elsewhere in this project.
 */
public final class HoneycombLattice {

    private HoneycombLattice() {
        // Namespace of physical constants and geometry only; not instantiable.
    }

    /** Carbon-carbon bond length, in angstroms. This is the fundamental length scale of the lattice. */
    public static final double BOND_LENGTH_ANGSTROM = 1.42;

    /**
     * Real-space Bravais lattice vectors of the triangular sublattice (i.e. the
     * translations that map the A sublattice onto itself, or equivalently the B
     * sublattice onto itself). Their magnitude, sqrt(3)*a, is the graphene lattice
     * constant, not the bond length a itself.
     */
    public static final Vector2 A1 =
            new Vector2(1.5 * BOND_LENGTH_ANGSTROM, 0.5 * Math.sqrt(3) * BOND_LENGTH_ANGSTROM);
    public static final Vector2 A2 =
            new Vector2(1.5 * BOND_LENGTH_ANGSTROM, -0.5 * Math.sqrt(3) * BOND_LENGTH_ANGSTROM);

    /**
     * The three vectors from an A atom to its three nearest-neighbor B atoms. These
     * are exactly the hopping vectors that enter the nearest-neighbor tight-binding
     * sum: an electron on an A site can hop to a B site at A_position + delta_i for
     * i = 1, 2, 3.
     */
    public static final Vector2 DELTA_1 =
            new Vector2(0.5 * BOND_LENGTH_ANGSTROM, 0.5 * Math.sqrt(3) * BOND_LENGTH_ANGSTROM);
    public static final Vector2 DELTA_2 =
            new Vector2(0.5 * BOND_LENGTH_ANGSTROM, -0.5 * Math.sqrt(3) * BOND_LENGTH_ANGSTROM);
    public static final Vector2 DELTA_3 = new Vector2(-BOND_LENGTH_ANGSTROM, 0.0);

    public static final Vector2[] NEAREST_NEIGHBOR_VECTORS = {DELTA_1, DELTA_2, DELTA_3};

    /**
     * The six next-nearest-neighbor vectors, which connect an atom to the six other
     * atoms of the <em>same</em> sublattice. These are just +/-A1, +/-A2, and
     * +/-(A1 - A2): next-nearest-neighbor hopping never changes sublattice, which is
     * why it appears on the diagonal of the Bloch Hamiltonian rather than the
     * off-diagonal, and why it is responsible for breaking electron-hole symmetry
     * (see {@link graphene.monolayer.MonolayerHamiltonian}).
     */
    public static final Vector2[] NEXT_NEAREST_NEIGHBOR_VECTORS = {
        A1, A1.scale(-1), A2, A2.scale(-1), A1.minus(A2), A2.minus(A1)
    };

    /**
     * Reciprocal lattice vectors, defined by a_i . b_j = 2*pi*delta_ij. These set the
     * periodicity of k-space: any physical quantity that depends on k (band energies,
     * wavefunctions up to a phase) must be invariant under k -> k + m*b1 + n*b2 for
     * integers m, n.
     */
    public static final Vector2 B1 = new Vector2(
            2.0 * Math.PI / (3.0 * BOND_LENGTH_ANGSTROM),
            2.0 * Math.PI / (Math.sqrt(3) * BOND_LENGTH_ANGSTROM));
    public static final Vector2 B2 = new Vector2(
            2.0 * Math.PI / (3.0 * BOND_LENGTH_ANGSTROM),
            -2.0 * Math.PI / (Math.sqrt(3) * BOND_LENGTH_ANGSTROM));

    /** Center of the hexagonal Brillouin zone. */
    public static final Vector2 GAMMA = Vector2.ZERO;

    /**
     * The K Dirac point, a corner of the hexagonal Brillouin zone. This is where the
     * valence (pi) and conduction (pi*) bands touch with zero gap and linear
     * (Dirac-cone) dispersion; see {@link graphene.dirac.DiracCone}. Its analytic
     * position, K = (2*pi/3a, 2*pi/(3*sqrt(3)*a)), follows directly from requiring
     * that the nearest-neighbor structure factor R(k) (see MonolayerHamiltonian)
     * vanish.
     */
    public static final Vector2 K =
            new Vector2(2.0 * Math.PI / (3.0 * BOND_LENGTH_ANGSTROM),
                    2.0 * Math.PI / (3.0 * Math.sqrt(3) * BOND_LENGTH_ANGSTROM));

    /**
     * The K' Dirac point, the other inequivalent corner of the Brillouin zone. K and
     * K' are related by time-reversal (kx, ky) -> (kx, -ky) and are not connected by
     * any reciprocal lattice vector, so they host physically distinct (but
     * degenerate in energy) Dirac cones.
     */
    public static final Vector2 K_PRIME =
            new Vector2(2.0 * Math.PI / (3.0 * BOND_LENGTH_ANGSTROM),
                    -2.0 * Math.PI / (3.0 * Math.sqrt(3) * BOND_LENGTH_ANGSTROM));

    /**
     * The M point: the midpoint of the Brillouin zone edge joining K and K'. Because
     * K and K' share the same kx coordinate in this convention, M falls exactly on
     * the kx axis. Gamma-M-K-Gamma is the conventional high-symmetry path used to
     * present graphene's band structure in the literature (e.g. Castro Neto et al.,
     * Fig. 3).
     */
    public static final Vector2 M = new Vector2(K.x(), 0.0);

    /**
     * All six corners of the hexagonal Brillouin zone, obtained by 120-degree
     * rotations of K and K' (each is a 3-fold-symmetric point of the lattice),
     * returned in counterclockwise angular order so callers can use them
     * directly as a polygon (e.g. to draw the zone outline or test whether a
     * point lies inside it).
     */
    public static Vector2[] brillouinZoneCorners() {
        Vector2[] corners = new Vector2[6];
        for (int i = 0; i < 3; i++) {
            double angle = i * 2.0 * Math.PI / 3.0;
            corners[2 * i] = rotate(K, angle);
            corners[2 * i + 1] = rotate(K_PRIME, angle);
        }
        java.util.Arrays.sort(corners, java.util.Comparator.comparingDouble(v -> Math.atan2(v.y(), v.x())));
        return corners;
    }

    /** The three Brillouin zone corners equivalent to K by the lattice's 3-fold rotation symmetry. */
    public static Vector2[] kPoints() {
        return new Vector2[] {K, rotate(K, 2.0 * Math.PI / 3.0), rotate(K, 4.0 * Math.PI / 3.0)};
    }

    /** The three Brillouin zone corners equivalent to K' by the lattice's 3-fold rotation symmetry. */
    public static Vector2[] kPrimePoints() {
        return new Vector2[] {K_PRIME, rotate(K_PRIME, 2.0 * Math.PI / 3.0), rotate(K_PRIME, 4.0 * Math.PI / 3.0)};
    }

    private static Vector2 rotate(Vector2 v, double angleRadians) {
        double cos = Math.cos(angleRadians);
        double sin = Math.sin(angleRadians);
        return new Vector2(cos * v.x() - sin * v.y(), sin * v.x() + cos * v.y());
    }
}
