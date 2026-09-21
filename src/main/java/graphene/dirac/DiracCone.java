package graphene.dirac;

import graphene.lattice.HoneycombLattice;
import graphene.lattice.Vector2;
import graphene.monolayer.MonolayerHamiltonian;

/**
 * The linear ("Dirac") approximation to graphene's pi bands near a Dirac point.
 *
 * <p>Physics reasoning: expanding the nearest-neighbor structure factor f(k)
 * (see {@link MonolayerHamiltonian}) to first order in q = k - K around a
 * Dirac point K shows that |f(K+q)| grows linearly in |q|, with a slope set
 * purely by the lattice geometry and the hopping energy t: there is no k^2
 * term, because the three nearest-neighbor phases exp(i k . delta_i)
 * interfere destructively at K in just the right way to cancel it. Keeping
 * only this linear term turns the tight-binding Hamiltonian into a
 * two-dimensional massless Dirac equation,
 * <pre>
 *     E(q) = +/- hbar * vF * |q|,      vF = 3 t a / (2 hbar),
 * </pre>
 * where the "speed of light" of the Dirac equation is replaced by the Fermi
 * velocity vF. This is the origin of graphene's much-discussed "relativistic"
 * low-energy electrons: they are not actually relativistic, but the honeycomb
 * lattice's geometry alone forces them to obey the same linear,
 * direction-independent dispersion relation as massless relativistic
 * particles (Wallace, Phys. Rev. 71, 622 (1947); Castro Neto, Guinea, Peres,
 * Novoselov &amp; Geim, Rev. Mod. Phys. 81, 109 (2009), Eq. 7).
 *
 * <p>Because it is only the leading term of a Taylor expansion, the linear
 * formula is accurate close to K and increasingly wrong further away, where
 * the hexagonal symmetry of the lattice, not the circular symmetry implied
 * by |q|, reasserts itself as "trigonal warping" of the equal-energy
 * contours (visible in {@code figures/monolayer_dirac_cone_zoom.png}).
 */
public final class DiracCone {

    /** Reduced Planck constant, in eV*s. */
    public static final double HBAR_EV_S = 6.582119569e-16;

    private static final double ANGSTROM_IN_METERS = 1e-10;

    private final MonolayerHamiltonian hamiltonian;
    private final Vector2 diracPoint;
    private final double fittedSlopeEvAngstrom;

    public DiracCone(MonolayerHamiltonian hamiltonian) {
        this(hamiltonian, HoneycombLattice.K);
    }

    public DiracCone(MonolayerHamiltonian hamiltonian, Vector2 diracPoint) {
        this.hamiltonian = hamiltonian;
        this.diracPoint = diracPoint;
        this.fittedSlopeEvAngstrom = fitSlope();
    }

    public Vector2 diracPoint() {
        return diracPoint;
    }

    /**
     * The textbook analytic Fermi velocity, vF = 3*t*a/(2*hbar), obtained by
     * Taylor-expanding f(k) to linear order around K by hand (see class
     * Javadoc). Independent of any numerics: this is a closed-form formula
     * in terms of the two inputs to the tight-binding model, t and a.
     */
    public double analyticFermiVelocity() {
        double bondLengthMeters = HoneycombLattice.BOND_LENGTH_ANGSTROM * ANGSTROM_IN_METERS;
        return 3.0 * hamiltonian.hoppingT() * bondLengthMeters / (2.0 * HBAR_EV_S);
    }

    /**
     * The Fermi velocity extracted directly from the full tight-binding
     * dispersion, by radially averaging E_+(K+q) over many directions at
     * several small radii |q| and fitting a line E = slope*|q| through the
     * origin (least squares). Deliberately computed independently of
     * {@link #analyticFermiVelocity()}, so the two serve as a check on each
     * other; see {@link DiracConeAnalysis} for the numeric comparison.
     */
    public double fitFermiVelocity() {
        return fittedSlopeEvAngstrom * ANGSTROM_IN_METERS / HBAR_EV_S;
    }

    /** The fitted linear-cone slope, in eV per angstrom^-1 (i.e. hbar*vF expressed in this project's units). */
    public double fittedSlopeEvAngstrom() {
        return fittedSlopeEvAngstrom;
    }

    /** The linear-cone energy prediction E(q) = fittedSlope * |q|, using the numerically fitted Fermi velocity. */
    public double linearApproximation(double qMagnitudeInverseAngstrom) {
        return fittedSlopeEvAngstrom * qMagnitudeInverseAngstrom;
    }

    /**
     * Least-squares slope (eV per angstrom^-1) of the conduction band E_+(K+q)
     * versus |q|, fit through the origin and averaged over direction so that
     * an accidental anisotropy at one particular angle cannot bias the fit.
     * The sampled radii (0.004-0.03 A^-1) are small compared to |K| ~ 1.7
     * A^-1, deep inside the regime where the linear approximation should
     * hold to high accuracy.
     */
    private double fitSlope() {
        int radialSamples = 5;
        int angularSamples = 24;
        double minRadius = 0.004;
        double maxRadius = 0.03;

        double numerator = 0.0;
        double denominator = 0.0;
        for (int i = 0; i < radialSamples; i++) {
            double r = radialSamples == 1 ? minRadius : minRadius + (maxRadius - minRadius) * i / (radialSamples - 1);
            for (int j = 0; j < angularSamples; j++) {
                double theta = 2.0 * Math.PI * j / angularSamples;
                Vector2 q = new Vector2(r * Math.cos(theta), r * Math.sin(theta));
                double conductionEnergy = hamiltonian.bandEnergies(diracPoint.plus(q))[1];
                numerator += r * conductionEnergy;
                denominator += r * r;
            }
        }
        return numerator / denominator;
    }
}
