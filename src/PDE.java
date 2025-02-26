import java.util.function.Function;

/**
 * A static class that solves partial differential equations via numerical methods.
 * <p>
 *     This class and all functions within were added after I submitted my application to your institution. All code
 *     within this project will be submitted and graded for MTH-49cr Numerical Differential Equations Modeling.
 * </p>
 * @author Keeler Spear
 * @version %I%, %G%
 * @since 1.0
 */
public class PDE {

    //Returns the coefficients of the fourier series on [-L, L] up to order n.
    //f(x) = a0 + a1sin(x) + a2cos(x)
    /**
     * Calculates the coefficients for a Fourier Series on [-L, L].
     * <p>
     * The coefficients will form the following trigonometric function:
     * a + bsin(x) + acos(x) + bsin(2x) + acos(2x) + ... + bsin(nx) + acos(nx)
     * </p>
     *
     * @param fnc The function to converted into a fourier series
     * @param l The length of the interval.
     * @param n The order of the Fourier Series.
     * @return A matrix with the coefficients for a Fourier Series on [-L, L].
     * @throws IllegalArgumentException If n is less than zero.
     */
    public static Matrix fourierSeries(Function fnc, double l, int n) {
        if (n < 0) {
            throw new IllegalArgumentException("The set's order must be at least zero!");
        }

        n = 2 * n; //Adjustment for fourier polynomial indexing

        double[] vals = new double[n + 1];
        Function[] basis = BasisFunctions.fourier(n, l);
        Function<Double, Double> temp = (x) -> (Double) fnc.apply(x) * (Double) basis[0].apply(x);
        double invL = 1.0 / l;

        vals[0] = 0.5 * invL * Calculus.integrate(temp, -l, l);

        for (int i = 1; i < vals.length; i++) {
            int I = i;
            temp = (x) -> (Double) fnc.apply(x) * (Double) basis[I].apply(x);
            vals[i] = invL * Calculus.integrate(temp, -l, l);
        }

        return new Matrix(vals);
    }

}
