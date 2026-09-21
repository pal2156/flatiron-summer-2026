package graphene.monolayer;

import graphene.lattice.HoneycombLattice;
import graphene.lattice.Vector2;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import javax.imageio.ImageIO;
import org.knowm.xchart.AnnotationLine;
import org.knowm.xchart.AnnotationText;
import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.BitmapEncoder.BitmapFormat;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler.LegendPosition;
import org.knowm.xchart.style.markers.SeriesMarkers;

/**
 * Regenerates the three monolayer band-structure figures described in the
 * README: (a) the full band structure over the two-dimensional Brillouin
 * zone (the centerpiece figure of this project), (b) bands along the
 * conventional Gamma-M-K-Gamma high-symmetry path, and (c) a zoomed-in view
 * of the dispersion near K showing the Dirac cone.
 *
 * <p>All three figures are generated from the same underlying object,
 * {@link MonolayerHamiltonian}: they are just different ways of slicing the
 * same E(k) surface (full 2D zone, a 1D path through it, or a small 2D patch
 * around one Dirac point).
 */
public final class MonolayerBandStructure {

    private MonolayerBandStructure() {
        // Entry point only.
    }

    public static void main(String[] args) throws IOException {
        Path figuresDir = Path.of("figures");
        Files.createDirectories(figuresDir);

        MonolayerHamiltonian model = new MonolayerHamiltonian(MonolayerHamiltonian.DEFAULT_HOPPING_T_EV);

        plotBrillouinZoneBandStructure(model, figuresDir.resolve("monolayer_band_structure_bz.png"));
        plotHighSymmetryPath(model, figuresDir.resolve("monolayer_bands_high_symmetry_path.png"));
        plotDiracConeZoom(model, figuresDir.resolve("monolayer_dirac_cone_zoom.png"));

        System.out.println("Monolayer figures written to " + figuresDir.toAbsolutePath());
    }

    // ------------------------------------------------------------------
    // (a) Full band structure over the 2D Brillouin zone: the centerpiece figure.
    // ------------------------------------------------------------------

    private static void plotBrillouinZoneBandStructure(MonolayerHamiltonian model, Path outputPath) throws IOException {
        int panelPixels = 560;
        double halfWidth = HoneycombLattice.K.norm() * 1.18;
        double vmax = 3.05 * model.hoppingT();
        double vmin = -vmax;

        BufferedImage valencePanel =
                BrillouinZoneHeatmap.renderPanel(model, 0, Vector2.ZERO, halfWidth, panelPixels, vmin, vmax, true);
        BufferedImage conductionPanel =
                BrillouinZoneHeatmap.renderPanel(model, 1, Vector2.ZERO, halfWidth, panelPixels, vmin, vmax, true);

        int leftMargin = 70;
        int panelGap = 60;
        int topMargin = 100;
        int bottomMargin = 80;
        int colorbarGap = 40;
        int colorbarWidth = 34;
        int rightMargin = 110;

        int width = leftMargin + 2 * panelPixels + panelGap + colorbarGap + colorbarWidth + rightMargin;
        int height = topMargin + panelPixels + bottomMargin;

        BufferedImage canvas = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = canvas.createGraphics();
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, width, height);

        // Title.
        g.setColor(new Color(30, 30, 30));
        g.setFont(new Font("SansSerif", Font.BOLD, 22));
        centerText(g, "Monolayer Graphene Band Structure over the Brillouin Zone", width / 2, 40);
        g.setFont(new Font("SansSerif", Font.PLAIN, 14));
        centerText(
                g,
                "Nearest-neighbor tight binding, t = " + model.hoppingT() + " eV, Dirac points K, K' marked",
                width / 2,
                64);

        int leftPanelX = leftMargin;
        int rightPanelX = leftMargin + panelPixels + panelGap;
        int panelY = topMargin;

        g.drawImage(valencePanel, leftPanelX, panelY, null);
        g.drawImage(conductionPanel, rightPanelX, panelY, null);

        drawPanelFrame(g, leftPanelX, panelY, panelPixels, Vector2.ZERO, halfWidth, "Valence band  π(k)");
        drawPanelFrame(g, rightPanelX, panelY, panelPixels, Vector2.ZERO, halfWidth, "Conduction band  π*(k)");

        markDiracPoints(g, leftPanelX, panelY, Vector2.ZERO, halfWidth, panelPixels);
        markDiracPoints(g, rightPanelX, panelY, Vector2.ZERO, halfWidth, panelPixels);

        int colorbarX = rightPanelX + panelPixels + colorbarGap;
        drawColorbar(g, colorbarX, panelY, colorbarWidth, panelPixels, vmin, vmax, "E (eV)");

        g.dispose();
        ImageIO.write(canvas, "png", outputPath.toFile());
    }

    /** Draws the hexagonal Brillouin-zone outline, a light square border, axis labels, and a panel title. */
    private static void drawPanelFrame(
            Graphics2D g, int panelX, int panelY, int panelPixels, Vector2 center, double halfWidth, String title) {

        g.setColor(new Color(60, 60, 60));
        g.setStroke(new BasicStroke(1.2f));
        g.drawRect(panelX, panelY, panelPixels, panelPixels);

        Vector2[] hexagon = HoneycombLattice.brillouinZoneCorners();
        GeneralPath path = new GeneralPath();
        for (int i = 0; i < hexagon.length; i++) {
            double px = panelX + BrillouinZoneHeatmap.pixelX(hexagon[i], center, halfWidth, panelPixels);
            double py = panelY + BrillouinZoneHeatmap.pixelY(hexagon[i], center, halfWidth, panelPixels);
            if (i == 0) {
                path.moveTo(px, py);
            } else {
                path.lineTo(px, py);
            }
        }
        path.closePath();
        g.setColor(new Color(40, 40, 40));
        g.setStroke(new BasicStroke(1.6f));
        g.draw(path);

        g.setFont(new Font("SansSerif", Font.BOLD, 15));
        g.setColor(new Color(20, 20, 20));
        centerText(g, title, panelX + panelPixels / 2, panelY - 12);

        g.setFont(new Font("SansSerif", Font.PLAIN, 12));
        centerText(g, "kₓ (Å⁻¹)", panelX + panelPixels / 2, panelY + panelPixels + 24);
        drawVerticalText(g, "kᵧ (Å⁻¹)", panelX - 44, panelY + panelPixels / 2);
    }

    private static void markDiracPoints(
            Graphics2D g, int panelX, int panelY, Vector2 center, double halfWidth, int panelPixels) {
        markPointGroup(g, panelX, panelY, center, halfWidth, panelPixels, HoneycombLattice.kPoints(), "K");
        markPointGroup(g, panelX, panelY, center, halfWidth, panelPixels, HoneycombLattice.kPrimePoints(), "K'");
    }

    private static void markPointGroup(
            Graphics2D g,
            int panelX,
            int panelY,
            Vector2 center,
            double halfWidth,
            int panelPixels,
            Vector2[] points,
            String label) {
        g.setFont(new Font("SansSerif", Font.BOLD, 13));
        for (Vector2 p : points) {
            double px = panelX + BrillouinZoneHeatmap.pixelX(p, center, halfWidth, panelPixels);
            double py = panelY + BrillouinZoneHeatmap.pixelY(p, center, halfWidth, panelPixels);
            g.setColor(Color.BLACK);
            g.fill(new Ellipse2D.Double(px - 4, py - 4, 8, 8));
            g.setColor(Color.WHITE);
            g.fill(new Ellipse2D.Double(px - 2.2, py - 2.2, 4.4, 4.4));

            g.setColor(new Color(20, 20, 20));
            double labelOffsetX = p.x() >= 0 ? 8 : -8 - g.getFontMetrics().stringWidth(label);
            g.drawString(label, (float) (px + labelOffsetX), (float) (py - 8));
        }
    }

    private static void drawColorbar(
            Graphics2D g, int x, int y, int width, int height, double vmin, double vmax, String label) {
        for (int row = 0; row < height; row++) {
            double value = vmax - (vmax - vmin) * row / (height - 1);
            g.setColor(Colormap.sample(value, vmin, vmax));
            g.fillRect(x, y + row, width, 1);
        }
        g.setColor(new Color(60, 60, 60));
        g.drawRect(x, y, width, height);

        g.setFont(new Font("SansSerif", Font.PLAIN, 12));
        g.setColor(new Color(20, 20, 20));
        int ticks = 5;
        for (int i = 0; i < ticks; i++) {
            double value = vmax - (vmax - vmin) * i / (ticks - 1);
            int ty = y + i * height / (ticks - 1);
            g.drawLine(x + width, ty, x + width + 5, ty);
            g.drawString(String.format("%.1f", value), x + width + 9, ty + 4);
        }

        g.setFont(new Font("SansSerif", Font.BOLD, 13));
        drawVerticalText(g, label, x + width + 62, y + height / 2);
    }

    // ------------------------------------------------------------------
    // (b) Bands along the standard Gamma -> M -> K -> Gamma high-symmetry path.
    // ------------------------------------------------------------------

    private static void plotHighSymmetryPath(MonolayerHamiltonian model, Path outputPath) throws IOException {
        Vector2[] waypoints = {HoneycombLattice.GAMMA, HoneycombLattice.M, HoneycombLattice.K, HoneycombLattice.GAMMA};
        String[] waypointLabels = {"Γ", "M", "K", "Γ"};
        int pointsPerSegment = 200;

        List<Double> distances = new ArrayList<>();
        List<Double> valence = new ArrayList<>();
        List<Double> conduction = new ArrayList<>();
        double[] waypointDistances = new double[waypoints.length];

        double cumulative = 0.0;
        for (int seg = 0; seg < waypoints.length - 1; seg++) {
            Vector2 start = waypoints[seg];
            Vector2 end = waypoints[seg + 1];
            double segmentLength = end.minus(start).norm();
            waypointDistances[seg] = cumulative;

            int steps = seg == waypoints.length - 2 ? pointsPerSegment + 1 : pointsPerSegment;
            for (int i = 0; i < steps; i++) {
                double t = (double) i / pointsPerSegment;
                Vector2 k = start.lerp(end, t);
                double[] energies = model.bandEnergies(k);
                distances.add(cumulative + t * segmentLength);
                valence.add(energies[0]);
                conduction.add(energies[1]);
            }
            cumulative += segmentLength;
        }
        waypointDistances[waypoints.length - 1] = cumulative;

        XYChart chart = new XYChartBuilder()
                .width(980)
                .height(650)
                .title("Monolayer Graphene: Bands along Γ → M → K → Γ")
                .xAxisTitle("Path through the Brillouin zone")
                .yAxisTitle("Energy (eV)")
                .build();

        double axisHeadroom = 3.35 * model.hoppingT();
        chart.getStyler().setLegendPosition(LegendPosition.OutsideE);
        chart.getStyler().setPlotGridLinesVisible(true);
        chart.getStyler().setXAxisTicksVisible(false);
        chart.getStyler().setChartTitleFont(new Font("SansSerif", Font.BOLD, 20));
        chart.getStyler().setAxisTitleFont(new Font("SansSerif", Font.PLAIN, 14));
        chart.getStyler().setChartBackgroundColor(Color.WHITE);
        chart.getStyler().setPlotBackgroundColor(Color.WHITE);
        chart.getStyler().setPlotBorderVisible(true);
        chart.getStyler().setYAxisMin(-axisHeadroom);
        chart.getStyler().setYAxisMax(axisHeadroom);

        double[] xArray = distances.stream().mapToDouble(Double::doubleValue).toArray();

        XYSeries valenceSeries =
                chart.addSeries("Valence band (π)", xArray, valence.stream().mapToDouble(Double::doubleValue).toArray());
        valenceSeries.setMarker(SeriesMarkers.NONE);
        valenceSeries.setLineColor(new Color(33, 102, 172));
        valenceSeries.setLineWidth(2.4f);

        XYSeries conductionSeries = chart.addSeries(
                "Conduction band (π*)", xArray, conduction.stream().mapToDouble(Double::doubleValue).toArray());
        conductionSeries.setMarker(SeriesMarkers.NONE);
        conductionSeries.setLineColor(new Color(178, 24, 43));
        conductionSeries.setLineWidth(2.4f);

        // Fermi level / Dirac-point energy.
        AnnotationLine fermiLevel = new AnnotationLine(0.0, false, false);
        chart.addAnnotation(fermiLevel);

        // Vertical guides at each high-symmetry point, with Greek/roman labels along the top.
        double yTop = 3.15 * model.hoppingT();
        for (int i = 0; i < waypoints.length; i++) {
            if (i != 0 && i != waypoints.length - 1) {
                chart.addAnnotation(new AnnotationLine(waypointDistances[i], true, false));
            }
            chart.addAnnotation(new AnnotationText(waypointLabels[i], waypointDistances[i], yTop, false));
        }

        BitmapEncoder.saveBitmap(chart, outputPath.toString().replace(".png", ""), BitmapFormat.PNG);
    }

    // ------------------------------------------------------------------
    // (c) Zoomed dispersion near K, showing the Dirac cone.
    // ------------------------------------------------------------------

    private static void plotDiracConeZoom(MonolayerHamiltonian model, Path outputPath) throws IOException {
        int panelPixels = 560;
        double halfWidth = 0.55; // angstrom^-1; small enough to stay well inside the zone around K
        double vmax = model.hoppingT() * 3.0 * halfWidth / HoneycombLattice.K.norm() * 1.9;
        double vmin = -vmax;

        BufferedImage conductionPanel =
                BrillouinZoneHeatmap.renderPanel(model, 1, HoneycombLattice.K, halfWidth, panelPixels, vmin, vmax, false);
        BufferedImage valencePanel =
                BrillouinZoneHeatmap.renderPanel(model, 0, HoneycombLattice.K, halfWidth, panelPixels, vmin, vmax, false);

        int leftMargin = 70;
        int panelGap = 60;
        int topMargin = 100;
        int bottomMargin = 80;
        int colorbarGap = 40;
        int colorbarWidth = 34;
        int rightMargin = 110;
        int width = leftMargin + 2 * panelPixels + panelGap + colorbarGap + colorbarWidth + rightMargin;
        int height = topMargin + panelPixels + bottomMargin;

        BufferedImage canvas = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = canvas.createGraphics();
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, width, height);

        g.setColor(new Color(30, 30, 30));
        g.setFont(new Font("SansSerif", Font.BOLD, 22));
        centerText(g, "Dirac Cone: Full Dispersion Zoomed Near K", width / 2, 40);
        g.setFont(new Font("SansSerif", Font.PLAIN, 14));
        centerText(g, "Concentric, circularly-symmetric contours signal the linear |k| dispersion", width / 2, 64);

        int leftPanelX = leftMargin;
        int rightPanelX = leftMargin + panelPixels + panelGap;
        int panelY = topMargin;

        g.drawImage(valencePanel, leftPanelX, panelY, null);
        g.drawImage(conductionPanel, rightPanelX, panelY, null);

        double[][] valenceGrid =
                BrillouinZoneHeatmap.computeEnergyGrid(model, 0, HoneycombLattice.K, halfWidth, panelPixels);
        double[][] conductionGrid =
                BrillouinZoneHeatmap.computeEnergyGrid(model, 1, HoneycombLattice.K, halfWidth, panelPixels);
        drawConeContours(g, leftPanelX, panelY, valenceGrid, vmin, 0.0);
        drawConeContours(g, rightPanelX, panelY, conductionGrid, 0.0, vmax);

        drawZoomPanelFrame(g, leftPanelX, panelY, panelPixels, HoneycombLattice.K, halfWidth, "Valence cone");
        drawZoomPanelFrame(g, rightPanelX, panelY, panelPixels, HoneycombLattice.K, halfWidth, "Conduction cone");

        markSinglePoint(g, leftPanelX, panelY, HoneycombLattice.K, halfWidth, panelPixels, "K");
        markSinglePoint(g, rightPanelX, panelY, HoneycombLattice.K, halfWidth, panelPixels, "K");

        int colorbarX = rightPanelX + panelPixels + colorbarGap;
        drawColorbar(g, colorbarX, panelY, colorbarWidth, panelPixels, vmin, vmax, "E (eV)");

        g.dispose();
        ImageIO.write(canvas, "png", outputPath.toFile());
    }

    /** Draws evenly-spaced iso-energy contour rings on top of a cone heatmap panel. */
    private static void drawConeContours(Graphics2D g, int panelX, int panelY, double[][] grid, double rangeMin, double rangeMax) {
        g.setStroke(new BasicStroke(1.0f));
        g.setColor(new Color(0, 0, 0, 80));
        int levelCount = 8;
        for (int i = 1; i < levelCount; i++) {
            double level = rangeMin + (rangeMax - rangeMin) * i / levelCount;
            for (var segment : ContourTracer.trace(grid, level)) {
                g.draw(new java.awt.geom.Line2D.Double(
                        panelX + segment.getX1(), panelY + segment.getY1(),
                        panelX + segment.getX2(), panelY + segment.getY2()));
            }
        }
    }

    private static void drawZoomPanelFrame(
            Graphics2D g, int panelX, int panelY, int panelPixels, Vector2 center, double halfWidth, String title) {
        g.setColor(new Color(60, 60, 60));
        g.setStroke(new BasicStroke(1.2f));
        g.drawRect(panelX, panelY, panelPixels, panelPixels);

        g.setFont(new Font("SansSerif", Font.BOLD, 15));
        g.setColor(new Color(20, 20, 20));
        centerText(g, title, panelX + panelPixels / 2, panelY - 12);

        g.setFont(new Font("SansSerif", Font.PLAIN, 12));
        centerText(
                g,
                String.format("kₓ - Kₓ (Å⁻¹), range ±%.2f", halfWidth),
                panelX + panelPixels / 2,
                panelY + panelPixels + 24);
        drawVerticalText(g, "kᵧ - Kᵧ (Å⁻¹)", panelX - 44, panelY + panelPixels / 2);
    }

    private static void markSinglePoint(
            Graphics2D g, int panelX, int panelY, Vector2 point, double halfWidth, int panelPixels, String label) {
        double px = panelX + BrillouinZoneHeatmap.pixelX(point, point, halfWidth, panelPixels);
        double py = panelY + BrillouinZoneHeatmap.pixelY(point, point, halfWidth, panelPixels);
        g.setColor(Color.BLACK);
        g.fill(new Ellipse2D.Double(px - 5, py - 5, 10, 10));
        g.setColor(Color.WHITE);
        g.fill(new Ellipse2D.Double(px - 2.5, py - 2.5, 5, 5));
        g.setFont(new Font("SansSerif", Font.BOLD, 14));
        g.setColor(new Color(20, 20, 20));
        g.drawString(label, (float) (px + 10), (float) (py - 10));
    }

    // ------------------------------------------------------------------
    // Small drawing helpers shared by all figures in this class.
    // ------------------------------------------------------------------

    private static void centerText(Graphics2D g, String text, int centerX, int baselineY) {
        int textWidth = g.getFontMetrics().stringWidth(text);
        g.drawString(text, centerX - textWidth / 2, baselineY);
    }

    private static void drawVerticalText(Graphics2D g, String text, int x, int centerY) {
        int textWidth = g.getFontMetrics().stringWidth(text);
        var transform = g.getTransform();
        g.rotate(-Math.PI / 2, x, centerY);
        g.drawString(text, x - textWidth / 2, centerY);
        g.setTransform(transform);
    }
}
