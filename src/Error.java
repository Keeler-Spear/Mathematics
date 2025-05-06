import java.util.function.Function;

/**
 * A static class that computes the error between two numerical objects.
 * <p>
 *     This class and all functions within were added after I submitted my application to your institution.
 * </p>
 * @author Keeler Spear
 * @version %I%, %G%
 * @since 1.0
 */

public class Error {

    final static int decimalPrecision = 2;

    private static double round(double val) {
        val = Math.round(val * Math.pow(10, decimalPrecision));
        return val * Math.pow(10, -decimalPrecision);
    }

    /**
     * Computes the absolute error between two values.
     *
     * @param exact The true value.
     * @param approx The approximation of the true value.
     * @return The absolute error between two values.
     */
    public static double absolute(double exact, double approx) {
        return Math.abs(exact - approx);
    }

    /**
     * Computes the absolute error between two matrices.
     *
     * @param exact The true values.
     * @param approx The approximations of the true values.
     * @return The absolute error between two values.
     */
    public static double absolute(Matrix exact, Matrix approx) {
        Matrix sol = LinearAlgebra.subtractMatrices(exact, approx);
        sol = LinearAlgebra.abs(sol);
        return LinearAlgebra.matrixSum(sol);
    }

    /**
     * Computes the relative error between two values.
     *
     * @param exact The true value.
     * @param approx The approximation of the true value.
     * @return The relative error between two values.
     */
    public static double relative(double exact, double approx) {
        return Math.abs(exact - approx) / Math.abs(exact);
    }
}
