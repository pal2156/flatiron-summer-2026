package graphene.lattice;

import java.util.Arrays;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;

/**
 * Diagonalizes the small complex Hermitian Bloch Hamiltonians H(k) that appear
 * throughout this project (the 2x2 monolayer Hamiltonian and the 4x4 Bernal
 * bilayer Hamiltonian).
 *
 * <p>EJML does not ship a native complex-matrix eigensolver, but a Hermitian
 * matrix H = A + iB (A real symmetric, B real antisymmetric with zero diagonal)
 * has exactly the same eigenvalues as the real symmetric matrix
 * <pre>
 *     M = [ A  -B ]
 *         [ B   A ]
 * </pre>
 * with each eigenvalue of H appearing twice in M. This follows because if
 * H v = E v for complex eigenvector v = p + i*q (p, q real), then direct
 * substitution shows that both real vectors (p; q) and (-q; p) are
 * eigenvectors of M with the same real eigenvalue E. Since H is Hermitian its
 * eigenvalues E are guaranteed real, so this "complexification" trick lets us
 * diagonalize H exactly using EJML's real symmetric eigenvalue routine, which
 * is numerically robust and, because it specializes on symmetric input,
 * always returns real eigenvalues rather than the general complex-conjugate
 * pairs a non-symmetric solver would report.
 */
public final class HermitianEigenSolver {

    private HermitianEigenSolver() {
        // Static utility; not instantiable.
    }

    /**
     * Diagonalizes the n x n complex Hermitian matrix H = real + i*imag.
     *
     * @param real the real part of H, an n x n symmetric matrix
     * @param imag the imaginary part of H, an n x n antisymmetric matrix with a zero diagonal
     * @return the n eigenvalues of H, sorted in ascending order
     */
    public static double[] eigenvalues(double[][] real, double[][] imag) {
        int n = real.length;
        DMatrixRMaj doubled = new DMatrixRMaj(2 * n, 2 * n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                doubled.set(i, j, real[i][j]);
                doubled.set(n + i, n + j, real[i][j]);
                doubled.set(i, n + j, -imag[i][j]);
                doubled.set(n + i, j, imag[i][j]);
            }
        }

        EigenDecomposition_F64<DMatrixRMaj> eig = DecompositionFactory_DDRM.eig(2 * n, false, true);
        if (!eig.decompose(doubled)) {
            throw new IllegalStateException("Eigenvalue decomposition failed to converge");
        }

        double[] doubledEigenvalues = new double[2 * n];
        for (int i = 0; i < 2 * n; i++) {
            doubledEigenvalues[i] = eig.getEigenvalue(i).getReal();
        }
        Arrays.sort(doubledEigenvalues);

        // Each eigenvalue of H shows up twice (adjacent after sorting); average
        // the pair to cancel floating-point round-off between the two copies.
        double[] result = new double[n];
        for (int i = 0; i < n; i++) {
            result[i] = 0.5 * (doubledEigenvalues[2 * i] + doubledEigenvalues[2 * i + 1]);
        }
        return result;
    }
}
