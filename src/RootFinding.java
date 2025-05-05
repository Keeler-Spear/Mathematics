import java.util.function.Function;

/**
 * A static class that finds the roots of functions via numerical methods.
 *
 * @author Keeler Spear
 * @version %I%, %G%
 * @since 1.0
 */
public class RootFinding {

    final static double TOL = 0.0001;
    final static int MAX_ITERATIONS = 100000;
    final static double FIRST_GUESS = 1.0;
    final static double SECOND_GUESS = 2.0;


    /**
     * Finds a root of the given function using Newton's Method.
     *
     * @param f The function whose root will be calculated.
     * @param fp The derivative of the function.
     * @param x0 The initial guess for the root.
     * @return A root of the provided function.
     */
    public static double newtonsMethod(Function<Double, Double> f, Function<Double, Double> fp, double x0) {
        int i = 0;

        while (Math.abs(f.apply(x0)) > TOL && i < MAX_ITERATIONS) {
            x0 = x0 - f.apply(x0) / fp.apply(x0);
            i++;
        }

        return x0;
    }

    /**
     * Finds a root of the given function using Newton's Method.
     *
     * @param f The function whose root will be calculated.
     * @param fp The derivative of the function.
     * @return A root of the provided function.
     */
    public static double newtonsMethod(Function<Double, Double> f, Function<Double, Double> fp) {
        return newtonsMethod(f, fp, FIRST_GUESS);
    }

    /**
     * Finds a root of the given function using the secant method.
     *
     * @param f The function whose root will be calculated.
     * @param x0 The initial guess for the root.
     * @param x1 The second guess for the root.
     * @return A root of the provided function.
     */
    public static double secantMethod(Function<Double, Double> f, double x0, double x1) {
        int i = 0;
        double f0 = f.apply(x0);
        double f1 = f.apply(x1);
        double temp;

        while (Math.abs(f1) > TOL && i < MAX_ITERATIONS) {
            temp = x1;
            x1 = x1 - f1 * (x1 - x0) / (f1 - f0);
            x0 = temp;
            f0 = f1;
            f1 = f.apply(x1);
            i++;
        }

        return x1;
    }

    /**
     * Finds a root of the given function using the secant method.
     *
     * @param f The function whose root will be calculated.
     * @return A root of the provided function.
     */
    public static double secantMethod(Function<Double, Double> f) {
        return secantMethod(f, FIRST_GUESS, SECOND_GUESS);
    }

}
