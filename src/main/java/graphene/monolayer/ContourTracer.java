package graphene.monolayer;

import java.awt.geom.Line2D;
import java.util.ArrayList;
import java.util.List;

/**
 * Traces iso-energy contour lines through a rectangular grid of band energies
 * (marching squares). This is what turns a plain color gradient into a
 * recognizable "cone": concentric, evenly-spaced, circular contour rings are
 * the visual signature of a linear |k| dispersion, in exactly the same way
 * contour lines on a topographic map reveal a conical hill.
 */
final class ContourTracer {

    private ContourTracer() {
        // Static utility.
    }

    /** Traces every segment of the {@code level} contour of {@code grid}, in pixel (col, row) coordinates. */
    static List<Line2D.Double> trace(double[][] grid, double level) {
        List<Line2D.Double> segments = new ArrayList<>();
        int rows = grid.length;
        int cols = grid[0].length;

        for (int row = 0; row < rows - 1; row++) {
            for (int col = 0; col < cols - 1; col++) {
                double topLeft = grid[row][col];
                double topRight = grid[row][col + 1];
                double bottomLeft = grid[row + 1][col];
                double bottomRight = grid[row + 1][col + 1];

                List<double[]> crossings = new ArrayList<>(4);
                addIfCrossing(crossings, col, row, col + 1, row, topLeft, topRight, level); // N edge
                addIfCrossing(crossings, col + 1, row, col + 1, row + 1, topRight, bottomRight, level); // E edge
                addIfCrossing(crossings, col, row + 1, col + 1, row + 1, bottomLeft, bottomRight, level); // S edge
                addIfCrossing(crossings, col, row, col, row + 1, topLeft, bottomLeft, level); // W edge

                if (crossings.size() == 2) {
                    segments.add(toSegment(crossings.get(0), crossings.get(1)));
                } else if (crossings.size() == 4) {
                    // Ambiguous saddle cell: use the cell-average value to pick a consistent pairing.
                    double average = (topLeft + topRight + bottomLeft + bottomRight) / 4.0;
                    if (average > level) {
                        segments.add(toSegment(crossings.get(0), crossings.get(3))); // N-W
                        segments.add(toSegment(crossings.get(1), crossings.get(2))); // E-S
                    } else {
                        segments.add(toSegment(crossings.get(0), crossings.get(1))); // N-E
                        segments.add(toSegment(crossings.get(2), crossings.get(3))); // S-W
                    }
                }
            }
        }
        return segments;
    }

    private static void addIfCrossing(
            List<double[]> crossings, int x1, int y1, int x2, int y2, double v1, double v2, double level) {
        if ((v1 - level >= 0) == (v2 - level >= 0)) {
            return; // both endpoints on the same side: this edge is not crossed
        }
        double t = (level - v1) / (v2 - v1);
        crossings.add(new double[] {x1 + t * (x2 - x1), y1 + t * (y2 - y1)});
    }

    private static Line2D.Double toSegment(double[] a, double[] b) {
        return new Line2D.Double(a[0], a[1], b[0], b[1]);
    }
}
