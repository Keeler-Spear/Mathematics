package Mathematics;

import java.util.function.Function;

/**
 * A static class that performs interpolation operations.
 *
 * @author Keeler Spear
 * @version %I%, %G%
 * @since 1.0
 */
public class Interpolation {

    /**
     * Interpolates n data points with a polynomial of degree n - 1 and evaluates it at the data points provided.
     *
     * @param x The set of x-values.
     * @param y The set of y values.
     * @param xVals The values for the polynomial interpolate to be evaluated at.
     * @return The matrix containing the evaluation of the polynomial interpolate at each x in xVals.
     * @throws IllegalArgumentException If the data sets do not have exactly one column each.
     * @throws IllegalArgumentException If the data sets are not the same length.
     */
    public static Matrix polynomial (Matrix x, Matrix y, Matrix xVals) {
        if (x.getCols() != 1 || y.getCols() != 1) {
            throw new IllegalArgumentException("Each data set must only have one column!");
        }

        if (x.getRows() != y.getRows()) {
            throw new IllegalArgumentException("Each data set must have the same length!");
        }

        int n = x.getRows();
        Matrix fnc = new Matrix(xVals.getRows(), 1);
        Function<Double, Double> tempF;

        for (int i = 1; i <= n; i++) {//For each L
            Matrix tempM = LinearAlgebra.constantMatrix(xVals.getRows(), 1, 1.0);
            for (int j = 1; j <= n; j++) {
                if (i != j) {
                    double xi = x.getValue(i, 1);
                    double xj = x.getValue(j, 1);
                    tempF = (t) -> (t - xj) / (xi - xj);
                    tempM = LinearAlgebra.multiplyValues(tempM, LinearAlgebra.applyFunction(xVals, tempF));
                }
            }

            fnc = LinearAlgebra.addMatrices(fnc, tempM, y.getValue(i, 1));
        }

        return fnc;
    }

    /**
     * Interpolates n data points with a polynomial of degree n - 1 and evaluates it at the data point provided.
     *
     * @param x The set of x-values.
     * @param y The set of y values.
     * @param xVal The value for the polynomial interpolate to be evaluated at.
     * @return The matrix containing the evaluation of the polynomial interpolate at each x in xVals.
     * @throws IllegalArgumentException If the data sets do not have exactly one column each.
     * @throws IllegalArgumentException If the data sets are not the same length.
     */
    public static double polynomial (Matrix x, Matrix y, double xVal) {
        if (x.getCols() != 1 || y.getCols() != 1) {
            throw new IllegalArgumentException("Each data set must only have one column!");
        }

        if (x.getRows() != y.getRows()) {
            throw new IllegalArgumentException("Each data set must have the same length!");
        }

        Matrix val = polynomial(x, y, new Matrix(new double[] {xVal}));

        return val.getValue(1, 1);
    }

    /**
     * Interpolates n data points with a polynomial of degree n - 1.
     *
     * @param x The set of x-values.
     * @param y The set of y values.
     * @return The matrix containing the evaluation of the polynomial interpolate at each x in xVals.
     * @throws IllegalArgumentException If the data sets do not have exactly one column each.
     * @throws IllegalArgumentException If the data sets are not the same length.
     */
    public static Function<Double, Double> polynomial (Matrix x, Matrix y) {
        if (x.getCols() != 1 || y.getCols() != 1) {
            throw new IllegalArgumentException("Each data set must only have one column!");
        }

        if (x.getRows() != y.getRows()) {
            throw new IllegalArgumentException("Each data set must have the same length!");
        }

        Function<Double, Double> fnc = (t) -> {
            double result = 0.0;

            for (int i = 1; i <= x.getRows(); i++) {
                double term = y.getValue(i, 1);
                for (int j = 1; j <= x.getRows(); j++) {
                    if (i != j) {
                        term *= (t - x.getValue(j, 1)) / (x.getValue(i, 1) - x.getValue(j, 1));
                    }
                }
                result += term;
            }

            return result;
        };

        return fnc;
    }

}
