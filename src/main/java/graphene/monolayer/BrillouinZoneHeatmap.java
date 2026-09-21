package graphene.monolayer;

import graphene.lattice.HoneycombLattice;
import graphene.lattice.Vector2;
import java.awt.image.BufferedImage;

/**
 * Rasterizes a single graphene band, E_-(k) or E_+(k), as a colored k-space
 * image. This is pure rendering machinery (no new physics beyond what
 * {@link MonolayerHamiltonian} already computes) shared by the full-zone
 * band-structure figure and the zoomed-in Dirac-cone figure, which are just
 * the same computation over a small versus a large k-space window.
 */
final class BrillouinZoneHeatmap {

    private BrillouinZoneHeatmap() {
        // Static utility.
    }

    /**
     * Renders {@code bandIndex} (0 = valence, 1 = conduction) of {@code hamiltonian}
     * over a square k-space window of half-width {@code halfWidth} centered at
     * {@code center}, mapped through the shared diverging colormap on [vmin, vmax].
     *
     * @param maskOutsideHexagon if true, pixels outside the first Brillouin zone are
     *     left transparent, so only the physically distinct k-points are shown
     */
    static BufferedImage renderPanel(
            MonolayerHamiltonian hamiltonian,
            int bandIndex,
            Vector2 center,
            double halfWidth,
            int pixels,
            double vmin,
            double vmax,
            boolean maskOutsideHexagon) {

        double[][] energy = computeEnergyGrid(hamiltonian, bandIndex, center, halfWidth, pixels);
        BufferedImage image = new BufferedImage(pixels, pixels, BufferedImage.TYPE_INT_ARGB);
        Vector2[] hexagon = HoneycombLattice.brillouinZoneCorners();

        for (int row = 0; row < pixels; row++) {
            double ky = kyAtRow(row, center, halfWidth, pixels);
            for (int col = 0; col < pixels; col++) {
                double kx = kxAtCol(col, center, halfWidth, pixels);

                if (maskOutsideHexagon && !insideHexagon(new Vector2(kx, ky), hexagon)) {
                    image.setRGB(col, row, 0); // fully transparent
                    continue;
                }

                int argb = Colormap.sample(energy[row][col], vmin, vmax).getRGB() | 0xFF000000;
                image.setRGB(col, row, argb);
            }
        }
        return image;
    }

    /**
     * Computes {@code bandIndex}'s energy on the same pixel grid used by
     * {@link #renderPanel}, without colorizing it. Exposed so that
     * {@link ContourTracer} can draw iso-energy contour lines exactly aligned
     * with the heatmap pixels underneath them.
     */
    static double[][] computeEnergyGrid(
            MonolayerHamiltonian hamiltonian, int bandIndex, Vector2 center, double halfWidth, int pixels) {
        double[][] grid = new double[pixels][pixels];
        for (int row = 0; row < pixels; row++) {
            double ky = kyAtRow(row, center, halfWidth, pixels);
            for (int col = 0; col < pixels; col++) {
                double kx = kxAtCol(col, center, halfWidth, pixels);
                grid[row][col] = hamiltonian.bandEnergies(new Vector2(kx, ky))[bandIndex];
            }
        }
        return grid;
    }

    private static double kxAtCol(int col, Vector2 center, double halfWidth, int pixels) {
        return center.x() - halfWidth + 2.0 * halfWidth * col / (pixels - 1);
    }

    private static double kyAtRow(int row, Vector2 center, double halfWidth, int pixels) {
        // Image row 0 is the top of the picture, but +ky points up in physics, so flip.
        return center.y() + halfWidth * (1.0 - 2.0 * row / (pixels - 1));
    }

    /** x pixel coordinate (within one panel) that a given k-point maps to under {@link #renderPanel}. */
    static double pixelX(Vector2 k, Vector2 center, double halfWidth, int pixels) {
        return (k.x() - (center.x() - halfWidth)) / (2.0 * halfWidth) * (pixels - 1);
    }

    /** y pixel coordinate (within one panel) that a given k-point maps to under {@link #renderPanel}. */
    static double pixelY(Vector2 k, Vector2 center, double halfWidth, int pixels) {
        return (center.y() + halfWidth - k.y()) / (2.0 * halfWidth) * (pixels - 1);
    }

    /** Standard convex-polygon point test: {@code point} is inside iff it is on the same side of every edge. */
    private static boolean insideHexagon(Vector2 point, Vector2[] hexagonCcw) {
        int n = hexagonCcw.length;
        for (int i = 0; i < n; i++) {
            Vector2 a = hexagonCcw[i];
            Vector2 b = hexagonCcw[(i + 1) % n];
            double edgeX = b.x() - a.x();
            double edgeY = b.y() - a.y();
            double cross = edgeX * (point.y() - a.y()) - edgeY * (point.x() - a.x());
            if (cross < 0) {
                return false;
            }
        }
        return true;
    }
}
