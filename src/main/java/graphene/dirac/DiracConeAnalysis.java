package graphene.dirac;

import graphene.lattice.HoneycombLattice;
import graphene.lattice.Vector2;
import graphene.monolayer.MonolayerHamiltonian;
import java.awt.Color;
import java.awt.Font;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import org.knowm.xchart.AnnotationLine;
import org.knowm.xchart.AnnotationText;
import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.BitmapEncoder.BitmapFormat;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.XYSeries.XYSeriesRenderStyle;
import org.knowm.xchart.style.Styler.LegendPosition;
import org.knowm.xchart.style.lines.SeriesLines;
import org.knowm.xchart.style.markers.SeriesMarkers;

/**
 * Compares the numerically fitted Fermi velocity to the analytic prediction,
 * and regenerates the figure showing where the linear Dirac-cone
 * approximation matches the full tight-binding dispersion, and where it
 * visibly departs from it.
 */
public final class DiracConeAnalysis {

    private DiracConeAnalysis() {
        // Entry point only.
    }

    public static void main(String[] args) throws IOException {
        Path figuresDir = Path.of("figures");
        Files.createDirectories(figuresDir);

        MonolayerHamiltonian hamiltonian = new MonolayerHamiltonian(MonolayerHamiltonian.DEFAULT_HOPPING_T_EV);
        DiracCone cone = new DiracCone(hamiltonian);

        double vFAnalytic = cone.analyticFermiVelocity();
        double vFFit = cone.fitFermiVelocity();
        double relativeDifferencePercent = 100.0 * Math.abs(vFFit - vFAnalytic) / vFAnalytic;

        System.out.println("Dirac cone Fermi velocity, t = " + hamiltonian.hoppingT() + " eV:");
        System.out.printf("  analytic   vF = 3ta/(2hbar) = %.4e m/s%n", vFAnalytic);
        System.out.printf("  fitted     vF (from E(K+q) near K) = %.4e m/s%n", vFFit);
        System.out.printf("  relative difference: %.6f%% (should be ~0: both describe the same linear term)%n",
                relativeDifferencePercent);

        plotConeOverlay(hamiltonian, cone, figuresDir.resolve("dirac_cone_linear_fit_overlay.png"));
        System.out.println("Dirac-cone figure written to " + figuresDir.toAbsolutePath());
    }

    private static void plotConeOverlay(MonolayerHamiltonian hamiltonian, DiracCone cone, Path outputPath)
            throws IOException {
        Vector2 k0 = cone.diracPoint();
        Vector2 towardGamma = HoneycombLattice.GAMMA.minus(k0);
        Vector2 direction = towardGamma.scale(1.0 / towardGamma.norm());

        double qMax = 0.9; // angstrom^-1; |K| ~ 1.7 angstrom^-1, so this stays well short of Gamma
        int samples = 300;
        double[] q = new double[2 * samples + 1];
        double[] tbConduction = new double[q.length];
        double[] tbValence = new double[q.length];

        for (int i = -samples; i <= samples; i++) {
            double t = qMax * i / samples;
            int idx = i + samples;
            Vector2 point = k0.plus(direction.scale(t));
            double[] energies = hamiltonian.bandEnergies(point);
            q[idx] = t;
            tbValence[idx] = energies[0];
            tbConduction[idx] = energies[1];
        }

        // The isotropic linear cone E = slope*|q| is exactly two straight rays meeting
        // at the origin, so 3 points fully describe it — using the same dense sampling
        // as the curved tight-binding bands would make the dashed stroke below restart
        // (and so look solid) at every tiny polyline segment.
        double edgeEnergy = cone.linearApproximation(qMax);
        double[] linearQ = {-qMax, 0.0, qMax};
        double[] linearUpper = {edgeEnergy, 0.0, edgeEnergy};
        double[] linearLower = {-edgeEnergy, 0.0, -edgeEnergy};

        double breakdownRadius = findBreakdownRadius(hamiltonian, cone, k0, direction, qMax);

        XYChart chart = new XYChartBuilder()
                .width(1000)
                .height(680)
                .title("Dirac Cone: Linear Fit vs. Full Tight-Binding Dispersion")
                .xAxisTitle("q = k - K along the K → Γ direction (Å⁻¹)")
                .yAxisTitle("Energy (eV)")
                .build();

        chart.getStyler().setLegendPosition(LegendPosition.OutsideE);
        chart.getStyler().setPlotGridLinesVisible(true);
        chart.getStyler().setChartTitleFont(new Font("SansSerif", Font.BOLD, 20));
        chart.getStyler().setAxisTitleFont(new Font("SansSerif", Font.PLAIN, 14));
        chart.getStyler().setChartBackgroundColor(Color.WHITE);
        chart.getStyler().setPlotBackgroundColor(Color.WHITE);
        chart.getStyler().setPlotBorderVisible(true);
        chart.getStyler().setXAxisMaxLabelCount(9);
        chart.getStyler().setxAxisTickLabelsFormattingFunction(v -> String.format("%.2f", v));

        XYSeries tbConductionSeries = chart.addSeries("Full dispersion, π*", q, tbConduction);
        styleLine(tbConductionSeries, new Color(178, 24, 43), 2.6f, XYSeriesRenderStyle.Line);

        XYSeries tbValenceSeries = chart.addSeries("Full dispersion, π", q, tbValence);
        styleLine(tbValenceSeries, new Color(33, 102, 172), 2.6f, XYSeriesRenderStyle.Line);

        XYSeries linearUpperSeries = chart.addSeries("Linear fit, +vF·ħ|q|", linearQ, linearUpper);
        styleLine(linearUpperSeries, new Color(40, 40, 40), 1.8f, XYSeriesRenderStyle.Line);
        linearUpperSeries.setLineStyle(SeriesLines.DASH_DASH);

        XYSeries linearLowerSeries = chart.addSeries("Linear fit, -vF·ħ|q|", linearQ, linearLower);
        styleLine(linearLowerSeries, new Color(40, 40, 40), 1.8f, XYSeriesRenderStyle.Line);
        linearLowerSeries.setLineStyle(SeriesLines.DASH_DASH);

        chart.addAnnotation(new AnnotationLine(0.0, false, false));
        chart.addAnnotation(new AnnotationLine(breakdownRadius, true, false));
        chart.addAnnotation(new AnnotationLine(-breakdownRadius, true, false));
        chart.addAnnotation(new AnnotationText(
                String.format("~10%% deviation at |q| = %.2f Å⁻¹", breakdownRadius),
                0.0,
                cone.linearApproximation(qMax) * 0.55,
                false));

        BitmapEncoder.saveBitmap(chart, outputPath.toString().replace(".png", ""), BitmapFormat.PNG);
    }

    private static void styleLine(XYSeries series, Color color, float width, XYSeriesRenderStyle style) {
        series.setXYSeriesRenderStyle(style);
        series.setMarker(SeriesMarkers.NONE);
        series.setLineColor(color);
        series.setLineWidth(width);
    }

    /**
     * Finds the smallest |q| (searching outward from K along {@code direction}) at
     * which the full tight-binding conduction band departs from the fitted linear
     * cone by more than 10%, i.e. where the Dirac approximation visibly breaks down.
     */
    private static double findBreakdownRadius(
            MonolayerHamiltonian hamiltonian, DiracCone cone, Vector2 k0, Vector2 direction, double searchMax) {
        int steps = 2000;
        for (int i = 1; i <= steps; i++) {
            double t = searchMax * i / steps;
            double linear = cone.linearApproximation(t);
            double full = hamiltonian.bandEnergies(k0.plus(direction.scale(t)))[1];
            if (Math.abs(full - linear) / linear > 0.10) {
                return t;
            }
        }
        return searchMax;
    }
}
