/**
 * A static class that performs calculus on data sets and functions using numerical methods.
 * <p>
 *     This class supports operations such as differentiation and integration.
 * </p>
 *
 * @author Keeler Spear
 * @version %I%, %G%
 * @since 1.0
 */
public class Calculus {
    private static final double TOL = 0.00000001;
    private static final double hDEF = 0.001;
    private static final int BASE_NUM_INTERVALS = 1000;

    //Computes if the provided value is "zero."
    private static boolean isZero(double val) {
        if (Math.abs(val) < TOL) {
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * Numerically calculates a function's derivative at the provided point with O(h^2) error, where h is {@value hDEF}.
     * <p>
     *     The derivative is calculated with the centered difference formula, due to its high accuracy, efficiency, and
     *     general invulnerability to roundoff error.
     * </p>
     * @param function The function that will be differentiated.
     * @param x0 The point the derivative will be evaluated at.
     * @return The function's derivative evaluation at the provided point.
     */
    public static double differentiate(Function function, double x0) {
        return centeredDifference(function.eval(x0 + hDEF), function.eval(x0 - hDEF), hDEF);
    }

    /**
     * Numerically calculates a function's derivative at all data points (assuming the data corresponds to a continuous
     * function) with O(h^2) error, where h is the distance between the x data.
     * <p>
     *     The derivative is calculated with the centered difference and one-sided three point formulas, due to their
     *     high accuracy, efficiency, and general invulnerability to roundoff error.
     * </p>
     *
     * @param x The set of equidistant x points for which f(x) data exists.
     * @param y The set of a function outputs from the provided inputs.
     * @throws IllegalArgumentException If every x point does not have an f(x) point and vice versa.
     * @throws IllegalArgumentException If there are an insufficient number (less than three) of data points.
     * @throws IllegalArgumentException If the x data is not equidistant.
     * @return An array of the function's derivative evaluations at the provided points.
     */
    public static double[] differentiate(double[] x, double[] y) { //Returns derivatives at each data point with O(h^2) error. This method assumes equally spaced data
        if (x.length != y.length) {
            throw new IllegalArgumentException("You must have the same number of data points in x and f(x)!");
        }
        if (x.length < 3) {
            throw new IllegalArgumentException("There are not enough data points to calculate the derivative!");
        }
        double dist = x[1] - x[0];
        for (int i = 1; i < x.length - 1; i++) {
            if (!isZero(dist - (x[i + 1] - x[i]))) {
                throw new IllegalArgumentException("The x points are not equidistant!");
            }
        }

        double[] derivative = new double[x.length];
        double h = x[1] - x[0];

        derivative[0] = leftSidedThreePoint(y[0], y[1], y[2], h);
        for (int i = 1; i < x.length - 1; i++) {
            derivative[i] = centeredDifference(y[i+1], y[i-1], h);
        }
        derivative[derivative.length - 1] = rightSidedThreePoint(y[y.length-1], y[y.length-2], y[y.length-3], h);
        return derivative;
    }

    //Derivative at a point with O(h^2) error. f1 = f(x0+h), fm1 = f(x0-h).
    private static double centeredDifference (double f1, double fm1, double h) {
        return (f1 - fm1) / (2 * h);
    }

    //Derivative at a point with O(h^2) error. f0 = f(x0), f1 = f(x0+h), f2 = f(x0+2h).
    private static double leftSidedThreePoint (double f0, double f1, double f2, double h) {
        return (4*f1 - f2 - 3*f0) / (2 * h);
    }

    //Derivative at a point with O(h^2) error. f0 = f(x0), fm1 = f(x0-h), fm2 = f(x0-2h).
    private static double rightSidedThreePoint (double f0, double fm1, double fm2, double h) {
        return (4*fm1 - fm2 - 3*f0) / (-2 * h);
    }

    /**
     * Numerically calculates a function's integral over the provided interval with O(h^4) error and precision 3,
     * where h is the size of the interval divided by {@value BASE_NUM_INTERVALS}.
     * <p>
     *     The integral is calculated with Simpson's Quadrature, due to its very high accuracy, precision of 3, and
     *     efficiency. Gaussian Quadrature was not chosen due to the computational cost of Legendre Polynomial root-finding.
     * </p>
     *
     * @param function The function that will be integrated.
     * @param a The lower bound of integration.
     * @param b The upper bound of integration.
     * @return The function's integral over the provided interval.
     */
    public static double integrate(Function function, double a, double b) {
        //Simpson's Composite Quadrature Method
        double area = 0.0;
        double h = (b - a) / (BASE_NUM_INTERVALS);
        for (double i = a; i < b; i += 2 * h) {
            area += simpsQuad(function.eval(i), function.eval(i + h), function.eval(i + 2 * h), i, i + 2 * h);
        }
        return area;
    }

    /**
     * Numerically calculates a function's integral over the provided interval with O(h^4) error and precision 3,
     * where h is the size of the interval divided by the number of data points minus one.
     * <p>
     *     The integral is calculated with Simpson's Quadrature, due to its very high accuracy, precision of 3, and
     *     efficiency. Gaussian Quadrature was not chosen due to the computational cost of Legendre Polynomial root-finding.
     * </p>
     *
     * @param y The set of a function outputs from equidistant inputs within the interval.
     * @param a The lower bound of integration.
     * @param b The upper bound of integration.
     * @throws IllegalArgumentException If every x point does not have an f(x) point and vice versa.
     * @throws IllegalArgumentException If there are an insufficient number (less than three) of data points.
     * @throws IllegalArgumentException If the x data lies outside the interval ({@code a} ≤ x ≤ {@code b}).
     * @throws IllegalArgumentException If the x data is not equidistant.
     * @return The function's integral over the provided interval.
     */
    public static double integrate(double[] x, double[] y, double a, double b) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("You must have the same number of data points in x and f(x)!");
        }
        if (y.length < 3) {
            throw new IllegalArgumentException("There are not enough data points to calculate the integral!");
        }
        for (int i = 0; i < x.length; i++) {
            if (a > x[i] || b < x[i]) {
                throw new IllegalArgumentException("The data is not within the interval!");
            }
        }
        double dist = x[1] - x[0];
        for (int i = 1; i < x.length - 1; i++) {
            if (!isZero(dist - (x[i + 1] - x[i]))) {
                throw new IllegalArgumentException("The x points are not equidistant!");
            }
        }
        double area = 0.0;
        double h = (b - a) / (y.length - 1);
        for (int i = 0; i < y.length - 2; i += 2) {
            area += simpsQuad(y[i], y[i + 1], y[i + 2], x[i], x[i + 2]);
        }
        return area;
    }

    //Integrates of the function over the interval with O(h^4) error.
    private static double simpsQuad(double f0, double f1, double f2, double a, double b) {
        return ((b - a) / 6) * (f0 + 4 * f1 + f2);
    }
}
