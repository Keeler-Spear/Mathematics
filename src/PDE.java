import java.util.function.BiFunction;
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

    /**
     * Solves the transport equation (ut + cux = 0) analytically.
     *
     * @param c The coefficient of the ux term.
     * @param f The initial condition for the PDE.
     * @return A function that is a solution to the transport equation.
     */
    public static BiFunction<Double, Double, Double> solveTrEq(double c, Function<Double, Double> f) {
        BiFunction<Double, Double, Double> sol = (x, t) -> f.apply(x - c * t);
        return sol;
    }

    /**
     * Solves the transport equation (ut + cux = 0) numerically with FTFS.
     *
     * @param c The coefficient of the ux term.
     * @param f The initial condition for the PDE.
     * @param x0 The starting spacial position.
     * @param xMax The ending spacial position.
     * @param t0 The starting temporal position
     * @param tMax The ending temporal position.
     * @param h The space step to be used.
     * @param k The time step to be used.
     * @return A function that is a solution to the transport equation.
     */
    public static Matrix solveTrEqFTFS(double c, Function<Double, Double> f, double x0, double xMax, double t0, double tMax, double h, double k) {
        double lambda = c * k / h;
        int iterations = (int) Math.round((tMax - t0) / k);

        //Building v0
        Matrix v0 = LinearAlgebra.linSpace(x0, xMax, h);
        int size = v0.getRows();
        v0 = LinearAlgebra.applyFunction(v0, f);

        //Building A
        Matrix A = new Matrix(size, size);
        for (int i = 1; i <= size; i++) {
            A.setValue(i, i, 1 + lambda);
        }
        for (int i = 1; i < size; i++) {
            A.setValue(i, i + 1, -lambda);
        }

        //Iterating the system in time
        return iterateSystem(A, v0, iterations);
    }

    /**
     * Solves the transport equation (ut + cux = 0) numerically with FTBS.
     *
     * @param c The coefficient of the ux term.
     * @param f The initial condition for the PDE.
     * @param x0 The starting spacial position.
     * @param xMax The ending spacial position.
     * @param t0 The starting temporal position
     * @param tMax The ending temporal position.
     * @param h The space step to be used.
     * @param k The time step to be used.
     * @return A function that is a solution to the transport equation.
     */
    public static Matrix solveTrEqFTBS(double c, Function<Double, Double> f, double x0, double xMax, double t0, double tMax, double h, double k) {
        double lambda = c * k / h;
        int iterations = (int) Math.round((tMax - t0) / k);

        //Building v0
        Matrix v0 = LinearAlgebra.linSpace(x0, xMax, h);
        int size = v0.getRows();
        v0 = LinearAlgebra.applyFunction(v0, f);

        //Building A
        Matrix A = new Matrix(size, size);
        for (int i = 1; i <= size; i++) {
            A.setValue(i, i, 1 - lambda);
        }
        for (int i = 2; i <= size; i++) {
            A.setValue(i, i - 1, lambda);
        }

        //Iterating the system in time
        return iterateSystem(A, v0, iterations);
    }

    /**
     * Solves the transport equation (ut + cux = 0) numerically with FTCS.
     *
     * @param c The coefficient of the ux term.
     * @param f The initial condition for the PDE.
     * @param x0 The starting spacial position.
     * @param xMax The ending spacial position.
     * @param t0 The starting temporal position
     * @param tMax The ending temporal position.
     * @param h The space step to be used.
     * @param k The time step to be used.
     * @return A function that is a solution to the transport equation.
     */
    public static Matrix solveTrEqFTCS(double c, Function<Double, Double> f, double x0, double xMax, double t0, double tMax, double h, double k) {
        double lambda = c * k / (2 * h);
        int iterations = (int) Math.round((tMax - t0) / k);

        //Building v0
        Matrix v0 = LinearAlgebra.linSpace(x0, xMax, h);
        int size = v0.getRows();
        v0 = LinearAlgebra.applyFunction(v0, f);

        //Building A
        Matrix A = new Matrix(size, size);
        for (int i = 1; i <= size; i++) {
            A.setValue(i, i, 1);
        }
        for (int i = 2; i <= size; i++) {
            A.setValue(i, i - 1, lambda);
        }
        for (int i = 1; i < size; i++) {
            A.setValue(i, i + 1, -lambda);
        }

        //Iterating the system in time
        return iterateSystem(A, v0, iterations);
    }

    private static Matrix iterateSystem(Matrix A, Matrix v, int n) {
        Matrix sol = LinearAlgebra.transpose(v);

        for (int i = 0; i < n; i++) {
            v = LinearAlgebra.multiplyMatrices(A, v);
            sol.addRowBottom(v.getMatrix());
        }

        return sol;
    }

}
