package graphene.lattice;

/**
 * An immutable 2D real vector, used throughout this project for both real-space
 * lattice positions and reciprocal-space wave vectors k = (kx, ky).
 *
 * <p>Kept deliberately minimal: graphene's electronic structure only ever needs
 * vector addition, scaling, dot products, and norms in two dimensions, so a
 * general-purpose linear algebra type would be overkill here.
 */
public record Vector2(double x, double y) {

    public static final Vector2 ZERO = new Vector2(0.0, 0.0);

    public Vector2 plus(Vector2 other) {
        return new Vector2(x + other.x, y + other.y);
    }

    public Vector2 minus(Vector2 other) {
        return new Vector2(x - other.x, y - other.y);
    }

    public Vector2 scale(double factor) {
        return new Vector2(x * factor, y * factor);
    }

    public double dot(Vector2 other) {
        return x * other.x + y * other.y;
    }

    public double norm() {
        return Math.sqrt(x * x + y * y);
    }

    /** Linear interpolation from this vector to {@code other}, at parameter t in [0, 1]. */
    public Vector2 lerp(Vector2 other, double t) {
        return new Vector2(x + t * (other.x - x), y + t * (other.y - y));
    }
}
