package graphene.bilayer;

import graphene.lattice.HoneycombLattice;
import graphene.lattice.Vector2;
import graphene.monolayer.MonolayerHamiltonian;
import java.awt.Color;
import java.awt.Font;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import org.knowm.xchart.AnnotationLine;
import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.BitmapEncoder.BitmapFormat;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler.LegendPosition;
import org.knowm.xchart.style.lines.SeriesLines;
import org.knowm.xchart.style.markers.SeriesMarkers;

/**
 * Regenerates the three bilayer graphene figures described in the README:
 * (a) unbiased bands near K, (b) biased bands near K showing the induced
 * gap, and (c) the gap-vs-bias curve: the centerpiece figure demonstrating
 * that, unlike single-layer graphene, Bernal bilayer graphene has a band gap
 * that can be tuned continuously with an applied field.
 */
public final class BilayerBandStructure {

    private static final double GAMMA_0_EV = MonolayerHamiltonian.DEFAULT_HOPPING_T_EV;
    private static final double GAMMA_1_EV = BilayerHamiltonian.DEFAULT_GAMMA_1_EV;

    private BilayerBandStructure() {
        // Entry point only.
    }

    public static void main(String[] args) throws IOException {
        Path figuresDir = Path.of("figures");
        Files.createDirectories(figuresDir);

        plotBandsNearK(
                new BilayerHamiltonian(GAMMA_0_EV, GAMMA_1_EV),
                "Unbiased Bilayer Bands Near K (two touch parabolically, two split at ±γ₁)",
                figuresDir.resolve("bilayer_bands_unbiased.png"));

        double demoBias = 0.10;
        plotBandsNearK(
                new BilayerHamiltonian(GAMMA_0_EV, GAMMA_1_EV, demoBias),
                String.format("Biased Bilayer Bands Near K (V = %.2f eV opens a gap)", demoBias),
                figuresDir.resolve("bilayer_bands_biased.png"));

        plotGapVsBias(figuresDir.resolve("bilayer_gap_vs_bias.png"));

        System.out.println("Bilayer figures written to " + figuresDir.toAbsolutePath());
    }

    // ------------------------------------------------------------------
    // (a) / (b): all four bands along a line through K.
    // ------------------------------------------------------------------

    private static void plotBandsNearK(BilayerHamiltonian model, String title, Path outputPath) throws IOException {
        Vector2 k0 = HoneycombLattice.K;
        Vector2 towardGamma = HoneycombLattice.GAMMA.minus(k0);
        Vector2 direction = towardGamma.scale(1.0 / towardGamma.norm());

        double qMax = 0.3; // angstrom^-1: zoomed in close enough that gamma3/gamma4 physics would matter, but we omit them
        int samples = 300;
        double[] q = new double[2 * samples + 1];
        double[] band0 = new double[q.length];
        double[] band1 = new double[q.length];
        double[] band2 = new double[q.length];
        double[] band3 = new double[q.length];

        for (int i = -samples; i <= samples; i++) {
            double t = qMax * i / samples;
            int idx = i + samples;
            double[] energies = model.bandEnergies(k0.plus(direction.scale(t)));
            q[idx] = t;
            band0[idx] = energies[0];
            band1[idx] = energies[1];
            band2[idx] = energies[2];
            band3[idx] = energies[3];
        }

        XYChart chart = new XYChartBuilder()
                .width(1050)
                .height(680)
                .title(title)
                .xAxisTitle("q = k - K along the K → Γ direction (Å⁻¹)")
                .yAxisTitle("Energy (eV)")
                .build();

        styleCommonChart(chart);
        chart.getStyler().setYAxisMin(-0.7);
        chart.getStyler().setYAxisMax(0.7);

        addLine(chart, "Split band (≈ -γ₁)", q, band0, new Color(120, 120, 120), 2.0f);
        addLine(chart, "Valence (non-dimer)", q, band1, new Color(33, 102, 172), 2.8f);
        addLine(chart, "Conduction (non-dimer)", q, band2, new Color(178, 24, 43), 2.8f);
        addLine(chart, "Split band (≈ +γ₁)", q, band3, new Color(120, 120, 120), 2.0f);

        chart.addAnnotation(new AnnotationLine(0.0, false, false));

        BitmapEncoder.saveBitmap(chart, outputPath.toString().replace(".png", ""), BitmapFormat.PNG);
    }

    // ------------------------------------------------------------------
    // (c) Band gap versus interlayer bias: the centerpiece figure.
    // ------------------------------------------------------------------

    private static void plotGapVsBias(Path outputPath) throws IOException {
        double biasMax = 0.15;
        int samples = 60;
        double[] bias = new double[2 * samples + 1];
        double[] trueGap = new double[bias.length];
        double[] naiveGap = new double[bias.length];

        for (int i = -samples; i <= samples; i++) {
            double v = biasMax * i / samples;
            int idx = i + samples;
            BilayerHamiltonian model = new BilayerHamiltonian(GAMMA_0_EV, GAMMA_1_EV, v);
            bias[idx] = v;
            trueGap[idx] = model.minimumGapNearDiracPoint(HoneycombLattice.K, 0.08, 160, 192);
            naiveGap[idx] = 2.0 * Math.abs(v);
        }

        XYChart chart = new XYChartBuilder()
                .width(950)
                .height(680)
                .title("Bilayer Graphene: Tunable Band Gap from an Interlayer Bias")
                .xAxisTitle("Interlayer bias V (eV)")
                .yAxisTitle("Band gap (eV)")
                .build();

        styleCommonChart(chart);
        chart.getStyler().setYAxisMin(0.0);

        XYSeries trueGapSeries = addLine(chart, "Minimum direct gap (tight binding)", bias, trueGap,
                new Color(178, 24, 43), 3.0f);
        trueGapSeries.setSmooth(true);

        XYSeries naiveGapSeries = addLine(chart, "2|V|, evaluated at K", bias, naiveGap, new Color(90, 90, 90), 1.8f);
        naiveGapSeries.setLineStyle(SeriesLines.DASH_DASH);

        chart.addAnnotation(new AnnotationLine(0.0, true, false));

        BitmapEncoder.saveBitmap(chart, outputPath.toString().replace(".png", ""), BitmapFormat.PNG);
    }

    // ------------------------------------------------------------------
    // Shared chart styling.
    // ------------------------------------------------------------------

    private static void styleCommonChart(XYChart chart) {
        chart.getStyler().setLegendPosition(LegendPosition.OutsideE);
        chart.getStyler().setPlotGridLinesVisible(true);
        chart.getStyler().setChartTitleFont(new Font("SansSerif", Font.BOLD, 18));
        chart.getStyler().setAxisTitleFont(new Font("SansSerif", Font.PLAIN, 14));
        chart.getStyler().setChartBackgroundColor(Color.WHITE);
        chart.getStyler().setPlotBackgroundColor(Color.WHITE);
        chart.getStyler().setPlotBorderVisible(true);
        chart.getStyler().setXAxisMaxLabelCount(9);
        chart.getStyler().setxAxisTickLabelsFormattingFunction(v -> String.format("%.2f", v));
    }

    private static XYSeries addLine(XYChart chart, String name, double[] x, double[] y, Color color, float width) {
        XYSeries series = chart.addSeries(name, x, y);
        series.setMarker(SeriesMarkers.NONE);
        series.setLineColor(color);
        series.setLineWidth(width);
        return series;
    }
}
