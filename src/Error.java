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

    /**
     * Computes the mean squared error between two data sets.
     *
     * @param exact The true data set.
     * @param approx The approximation of the true data set.
     * @return The mean squared error between two data sets.
     * @throws IllegalArgumentException If the data sets provided are not the same length.
     */
    public static double meanSquared(double[] exact, double[] approx) {
        if (exact.length != approx.length) {
            throw new IllegalArgumentException("The data sets must be the same length!");
        }
        int n = exact.length;
        double sum = 0.0;

        for (int i = 0; i < n; i++) {
            sum += Math.pow(exact[i] - approx[i], 2);
        }

        return sum / n;
    }

    /**
     * Computes the mean squared error between two data sets.
     *
     * @param exact The true data set.
     * @param approx The approximation of the true data set.
     * @return The mean squared error between two data sets.
     * @throws IllegalArgumentException If the data sets provided have more than one column.
     * @throws IllegalArgumentException If the data sets provided are not the same length.
     */
    public static double meanSquared(Matrix exact, Matrix approx) {
        if (exact.getCols() != 1 || approx.getCols() != 1) {
            throw new IllegalArgumentException("The data sets must have one column each!");
        }

        if (exact.getRows() != approx.getRows()) {
            throw new IllegalArgumentException("The data sets must be the same length!");
        }

        int n = exact.getRows();
        double sum = 0.0;

        for (int i = 1; i <= n; i++) {
            sum += Math.pow(exact.getValue(i, 1) - approx.getValue(i, 1), 2);
        }

        return sum / n;
    }

    /**
     * Computes the mean squared error between two data sets, where the approximate data set will be created from a
     * polynomial created by the weights provided.
     *
     * @param x The x values of the exact data.
     * @param exact The true data set.
     * @param w The weights of the polynomial from which an approximation will be created.
     * @return The mean squared error between two data sets.
     * @throws IllegalArgumentException If the data sets provided have more than one column.
     * @throws IllegalArgumentException If the data sets provided are not the same length.
     */
    public static double meanSquared(Matrix x, Matrix exact, Matrix w, Function[] fnc) {
        if (exact.getCols() != 1 || x.getCols() != 1) {
            throw new IllegalArgumentException("The data sets must have one column each!");
        }

        if (exact.getRows() != x.getRows()) {
            throw new IllegalArgumentException("The data sets must be the same length!");
        }

        Matrix approx = LinReg.buildFunction(x, w, fnc);

        int n = exact.getRows();
        double sum = 0.0;

        for (int i = 1; i <= n; i++) {
            sum += Math.pow(exact.getValue(i, 1) - approx.getValue(i, 1), 2);
        }

        return sum / n;
    }

    //ToDo: r^2 or coefficient of determination
}
