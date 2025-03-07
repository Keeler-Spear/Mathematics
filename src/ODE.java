import java.util.function.BiFunction;
import java.util.function.Function;

/**
 * A static class that solves ordinary differential equations via numerical methods.
 * <p>
 *     This class and all functions within were added after I submitted my application to your institution. All code
 *     within this project will be submitted and graded for MTH-49cr Numerical Differential Equations Modeling.
 * </p>
 * @author Keeler Spear
 * @version %I%, %G%
 * @since 1.0
 */
public class ODE {

    private static final double BASE_VAL = 9998; //The present value matrices are filled with
    private static final int MAX_SIZE = 10000;

    public static Matrix solveIVP(BiFunction<Double, Double, Double> ode, double t0, double y0, double t, double h) {
        return pc(ode, t0, y0, t, h);
    }

    public static Matrix solveIVP(TriFunction<Double, Double, Double, Double> ode, double t0, double y0, double yp0, double t, double h) {
        TriFunction<Double, Double, Double, Double>[] odes = decomposeODE(ode);
        return eulerSystem(odes[0], odes[1], t0, y0, yp0, t, h);
    }

    //Linear shooting method
    public static Matrix solveBVP (TriFunction<Double, Double, Double, Double> ode, double t0, double y0, double t1, double y1, double h) {
        Matrix yp = solveIVP(ode, t0, y0, 0, t1, h);

        TriFunction<Double, Double, Double, Double> homoODE = (x, y, yprime) -> ode.apply(x, y, yprime) - ode.apply(x, 0.0, 0.0);

        Matrix yc = solveIVP(homoODE, t0, 0, 1, t1, h);

        double ypb  = yp.getValue(yp.getRows(), 1);
        double ycb  = yc.getValue(yc.getRows(), 1);

        return LinearAlgebra.addMatrices(yp, yc, (y1 - ypb) / ycb);
    }

    /**
     * Numerically solves the first-order ordinary differential equation over the provided interval using Euler's method
     * with O(h) global error, where h is the step size of t.
     *
     * @param ode The first-order ordinary differential equation f(t, y) to be solved.
     * @param t0 The left endpoint of the interval.
     * @param y0 The y value corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param h The step size of t.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     */
    public static Matrix euler (BiFunction<Double, Double, Double> ode, double t0, double y0, double t, double h) {
        Matrix y = LinearAlgebra.constantMatrix(generateYLength(t0, t, h) + 50, 1, BASE_VAL);
        y.setValue(1, 1, y0);
        int i = 2;

        while (t0 < t) {
            y0 = y0 + h * ode.apply(t0, y0);

            y.setValue(i, 1, y0);
            t0 += h;
            i++;
        }

        return trim(y);
    }

    /**
     * Numerically solves the first-order ordinary differential equation over the provided interval using the
     * Runge–Kutta 4 method with O(h^4) global error, where h is the step size of t.
     *
     * @param ode The first-order ordinary differential equation f(t, y) to be solved.
     * @param t0 The left endpoint of the interval.
     * @param y0 The y value corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param h The step size of t.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     */
    public static Matrix rk4 (BiFunction<Double, Double, Double> ode, double t0, double y0, double t, double h) {
        Matrix y = LinearAlgebra.constantMatrix(generateYLength(t0, t, h) + 50, 1, BASE_VAL);
        y.setValue(1, 1, y0);
        double k1;
        double k2;
        double k3;
        double k4;
        int i = 2;

        while (t0 < t) {
            k1 = ode.apply(t0, y0);
            k2 = ode.apply(t0 + 0.5 * h, y0 + 0.5 * h * k1);
            k3 = ode.apply(t0 + 0.5 * h, y0 + 0.5 * h * k2);
            k4 = ode.apply(t0 + h, y0 + h * k3);

            y0 = y0 + h * ((1.0 / 6) * k1 + (1.0 / 3) * k2 + (1.0 / 3) * k3 + (1.0 / 6) * k4);

            y.setValue(i, 1, y0);
            t0 += h;
            i++;
        }

        return trim(y);
    }

    /**
     * Numerically solves the system of first-order ordinary differential equations over the provided interval using
     * Euler's method with O(h) global error, where h is the step size of t. To solve a second order ordinary
     * differential equation, decompose it into a system of two first-order ordinary differential equations.
     *
     * @param ode1 The first first-order ordinary differential equation.
     * @param ode2 The second first-order ordinary differential equation.
     * @param t0 The left endpoint of the interval.
     * @param y0 The y values corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param h The step size of t.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     */
    public static Matrix eulerSystem (TriFunction<Double, Double, Double, Double> ode1, TriFunction<Double, Double, Double, Double> ode2, double t0, double y0, double yp0, double t, double h) {
        Matrix yVals = LinearAlgebra.constantMatrix(generateYLength(t0, t, h) + 50, 2, BASE_VAL);
        yVals.setValue(1, 1, y0);
        yVals.setValue(1, 2, yp0);
        Matrix yi = new Matrix(new double[] {y0, yp0});
        Matrix yP = new Matrix(new double[] {y0, yp0});
        int i = 2;

        while (t0 < t) {
            //Finding y'(t)
            yP.setValue(1, 1, ode1.apply(t0, yi.getValue(1,1), yi.getValue(2, 1)));
            yP.setValue(2, 1, ode2.apply(t0, yi.getValue(1,1), yi.getValue(2, 1)));
            //finding y(t)
            yi = LinearAlgebra.addMatrices(yi, yP, h);

            yVals.setValue(i, 1, yi.getValue(1,1));
            yVals.setValue(i, 2, yi.getValue(2,1));
            t0 += h;
            i++;
        }

        return trim(yVals);
    }

    /**
     * Numerically solves the system of first-order ordinary differential equations over the provided interval using the
     * Runge–Kutta 4 method with O(h^4) global error, where h is the step size of t. To solve a second order ordinary
     * differential equation, decompose it into a system of two first-order ordinary differential equations.
     *
     * @param ode1 The first first-order ordinary differential equation.
     * @param ode2 The second first-order ordinary differential equation.
     * @param t0 The left endpoint of the interval.
     * @param y0 The y values corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param h The step size of t.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     */
    public static Matrix rk4System (TriFunction<Double, Double, Double, Double> ode1, TriFunction<Double, Double, Double, Double> ode2, double t0, double y0, double yp0, double t, double h) {
        Matrix yVals = LinearAlgebra.constantMatrix(generateYLength(t0, t, h) + 50, 2, BASE_VAL);
        Matrix yi = new Matrix(new double[] {y0, yp0});
        yVals.setValue(1, 1, y0);
        yVals.setValue(1, 2, yp0);
        int i = 2;
        Matrix k1 = new Matrix(2, 1);
        Matrix k2 = new Matrix(2, 1);
        Matrix k3 = new Matrix(2, 1);
        Matrix k4 = new Matrix(2, 1);
        Matrix kTemp; //Used as temporary storage for wi + h/2 ki
        Matrix wkSum;

        while (t0 < t) {
            //Finding the k's
            //k1:
            k1.setValue(1, 1, ode1.apply(t0, yi.getValue(1,1), yi.getValue(2, 1)));
            k1.setValue(2, 1, ode2.apply(t0, yi.getValue(1,1), yi.getValue(2, 1)));
            kTemp = LinearAlgebra.addMatrices(yi, k1, 0.5 * h);
            //k2:
            k2.setValue(1, 1, ode1.apply(t0 + 0.5 * h, kTemp.getValue(1,1), kTemp.getValue(2, 1)));
            k2.setValue(2, 1, ode2.apply(t0 + 0.5 * h, kTemp.getValue(1,1), kTemp.getValue(2, 1)));
            kTemp = LinearAlgebra.addMatrices(yi, k2, 0.5 * h);
            //k3:
            k3.setValue(1, 1, ode1.apply(t0 + 0.5 * h, kTemp.getValue(1,1), kTemp.getValue(2, 1)));
            k3.setValue(2, 1, ode2.apply(t0 + 0.5 * h, kTemp.getValue(1,1), kTemp.getValue(2, 1)));
            kTemp = LinearAlgebra.addMatrices(yi, k3, h);
            //k4:
            k4.setValue(1, 1, ode1.apply(t0 + h, kTemp.getValue(1,1), kTemp.getValue(2, 1)));
            k4.setValue(2, 1, ode2.apply(t0 + h, kTemp.getValue(1,1), kTemp.getValue(2, 1)));

            //Multiplying k's by their respective weights and adding them
            wkSum = LinearAlgebra.addMatrices(LinearAlgebra.scaleMatrix(k1, 1.0/6), k2, 1.0/3);
            wkSum = LinearAlgebra.addMatrices(wkSum, k3, 1.0/3);
            wkSum = LinearAlgebra.addMatrices(wkSum, k4, 1.0/6);

            yi = LinearAlgebra.addMatrices(yi, wkSum, h);

            yVals.setValue(i, 1, yi.getValue(1,1));
            yVals.setValue(i, 2, yi.getValue(2,1));
            t0 += h;
            i++;
        }

        return trim(yVals);
    }

    /**
     * Numerically solves the first-order ordinary differential equation over the provided interval using the
     * Runge–Kutta 4 and Adams-Bashforth 4 methods with O(h^4) global error, where h is the step size of t.
     *
     * @param ode The first-order ordinary differential equation f(t, y) to be solved.
     * @param t0 The left endpoint of the interval.
     * @param y0 The y value corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param h The step size of t.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     */
    public static Matrix adamsBash (BiFunction<Double, Double, Double> ode, double t0, double y0, double t, double h) {
        Matrix rk4 = rk4(ode, t0, y0, t0 + 3 * h, h);
        double[] f = new double[rk4.getRows() - 1];

        for (int i = 0; i < f.length; i++) {
            f[i] = ode.apply(t0, rk4.getValue(i + 1, 1));
            t0 += h;
        }

        Matrix sol = mergeMatrices(rk4,adamsBash(ode, t0, rk4.getValue(rk4.getRows(), 1),t, h, f[2], f[1], f[0]));

        return sol;
    }

    /**
     * Numerically solves the first-order ordinary differential equation over the provided interval using the
     * Adams-Bashforth 4 method with O(h^4) global error, where h is the step size of t.
     *
     * @param ode The first-order ordinary differential equation f(t, y) to be solved.
     * @param t0 The left endpoint of the interval.
     * @param y0 The y value corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param h The step size of t.
     * @param fm1 The derivative of the function at t0 - h.
     * @param fm2 The derivative of the function at t0 - 2h.
     * @param fm3 The derivative of the function at t0 - 3h.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     */
    public static Matrix adamsBash (BiFunction<Double, Double, Double> ode, double t0, double y0, double t, double h, double fm1, double fm2, double fm3) {
        Matrix y = LinearAlgebra.constantMatrix(generateYLength(t0, t, h) + 50, 1, BASE_VAL);
        double[] f = {0, fm3, fm2, fm1}; //I could use a stack instead but I have concerns regarding its efficiency.

        y.setValue(1, 1, ode.apply(t0, y0));
        int i = 2;

        while (t0 < t) {
            //Shifting f array to the left
            for (int j = 0; j < f.length - 1; j++) {
                f[j] = f[j + 1];
            }
            f[3] = ode.apply(t0, y0);

            y0 = y0 + (h / 24.0) * (55 * f[3] - 59 * f[2] + 37 * f[1] - 9 * f[0]);

            y.setValue(i, 1, y0);
            t0 += h;
            i++;
        }

        return trim(y);
    }

    /**
     * Numerically solves the first-order ordinary differential equation over the provided interval using the
     * Runge–Kutta 4 and Predictor-Corrctor 4 methods with O(h^4) global error, where h is the step size of t.
     *
     * @param ode The first-order ordinary differential equation f(t, y) to be solved.
     * @param t0 The left endpoint of the interval.
     * @param y0 The y value corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param h The step size of t.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     */
    public static Matrix pc (BiFunction<Double, Double, Double> ode, double t0, double y0, double t, double h) {
        Matrix rk4 = rk4(ode, t0, y0, t0 + 3 * h, h);
        double[] f = new double[rk4.getRows() - 1];

        for (int i = 0; i < f.length; i++) {
            f[i] = ode.apply(t0, rk4.getValue(i + 1, 1));
            t0 += h;
        }

        return mergeMatrices(rk4,pc(ode, t0, rk4.getValue(rk4.getRows(), 1),t, h, f[2], f[1], f[0]));
    }

    /**
     * Numerically solves the first-order ordinary differential equation over the provided interval using the
     * Predictor-Corrector 4 method with O(h^4) global error, where h is the step size of t.
     *
     * @param ode The first-order ordinary differential equation f(t, y) to be solved.
     * @param t0 The left endpoint of the interval.
     * @param y0 The y value corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param h The step size of t.
     * @param fm1 The derivative of the function at t0 - h.
     * @param fm2 The derivative of the function at t0 - 2h.
     * @param fm3 The derivative of the function at t0 - 3h.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     */
    public static Matrix pc (BiFunction<Double, Double, Double> ode, double t0, double y0, double t, double h, double fm1, double fm2, double fm3) {
        Matrix y = LinearAlgebra.constantMatrix(generateYLength(t0, t, h) + 50, 1, BASE_VAL);
        double[] f = {fm3, fm2, fm1, ode.apply(t0, y0)}; //I could use a stack instead but I have concerns regarding its efficiency.

        y.setValue(1, 1, y0);
        double ym1;
        int i = 2;

        while (t0 < t) {
            f[3] = ode.apply(t0, y0);

            ym1 = y0;
            y0 = y0 + (h / 24.0) * (55 * f[3] - 59 * f[2] + 37 * f[1] - 9 * f[0]);
            t0 += h;

            //Shifting the array to the left
            for (int j = 0; j < f.length - 1; j++) {
                f[j] = f[j + 1];
            }
            f[3] = ode.apply(t0, y0);

            y0 = ym1 + (h / 24.0) * (9 * f[3] + 19 * f[2] - 5 * f[1] + f[0]);

            y.setValue(i, 1, y0);
            i++;
        }

        return trim(y);
    }

    //Input should be y'' = a(x)y' + b(x)y + c. Params: (x, y0, y1);
    public static TriFunction[] decomposeODE(TriFunction<Double, Double, Double, Double> ode) {
        TriFunction<Double, Double, Double, Double>[] fncs = new TriFunction[2];
        fncs[0] = (x, y, yp) -> yp;
        fncs[1] = ode;

        return fncs;
    }

    //Makes a new array [A B] and removes the duplicate left-bound of A
    private static Matrix mergeMatrices(Matrix A, Matrix B) {
        double[] C = new double[A.getRows() + B.getRows() - 1];

        for (int i = 0; i < A.getRows() - 1; i++) {
            C[i] = A.getValue(i + 1, 1);
        }

        for (int i = 0; i < B.getRows(); i++) {
            C[i + A.getRows() - 1] = B.getValue(i + 1, 1);
        }

        return new Matrix(C);
    }

    //I need this because the array generation is janky
    private static int generateYLength(double t0, double t, double h) {
        int n = 1 + (int) ((t - t0)/h);

        if (h > 0.5) {
            n += 1;
        }

        return n;
    }

    //"Trims" the matrix of the rows not being used
    private static Matrix trim(Matrix X) {
        int numRows = 0;
        while (X.getValue(numRows + 1, 1) != BASE_VAL) {
            numRows += 1;
        }

        Matrix A = new Matrix (numRows, X.getCols());

        for (int i = 1; i <= A.getRows(); i++) {
            for (int j = 1; j <= A.getCols(); j++) {
                A.setValue(i, j, X.getValue(i, j));
            }
        }

        return A;
    }
}
