package graphene.monolayer;

import java.awt.Color;

/**
 * A small diverging (blue - white - red) colormap used to render the band-energy
 * heatmaps.
 *
 * <p>This is deliberately not physics, just visualization plumbing: graphene's
 * two pi bands are symmetric about zero energy (for the pure nearest-neighbor
 * model), so a colormap that passes through white at the midpoint of the
 * plotted range makes the zero-energy contour, which runs straight through
 * the Dirac points, visually pop out on its own, before any point markers are
 * even drawn.
 *
 * <p>Control points follow the ColorBrewer "RdBu" diverging palette, a standard
 * choice for signed scientific data.
 */
final class Colormap {

    private static final Color[] STOPS = {
        new Color(33, 102, 172),   // deep blue: most negative (valence band bottom)
        new Color(146, 197, 222),  // pale blue
        new Color(247, 247, 247),  // white: zero energy (Dirac points live here)
        new Color(244, 165, 130),  // pale red
        new Color(178, 24, 43)     // deep red: most positive (conduction band top)
    };

    private Colormap() {
        // Static utility.
    }

    /** Maps a value in [min, max] to a color; values outside the range are clamped. */
    static Color sample(double value, double min, double max) {
        double t = max > min ? (value - min) / (max - min) : 0.5;
        t = Math.max(0.0, Math.min(1.0, t));

        double scaled = t * (STOPS.length - 1);
        int lower = (int) Math.floor(scaled);
        int upper = Math.min(lower + 1, STOPS.length - 1);
        double frac = scaled - lower;

        Color a = STOPS[lower];
        Color b = STOPS[upper];
        int r = (int) Math.round(a.getRed() + frac * (b.getRed() - a.getRed()));
        int g = (int) Math.round(a.getGreen() + frac * (b.getGreen() - a.getGreen()));
        int bl = (int) Math.round(a.getBlue() + frac * (b.getBlue() - a.getBlue()));
        return new Color(r, g, bl);
    }
}
