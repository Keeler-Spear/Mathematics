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
