import java.util.function.BiFunction;

/**
 * A static class that solves ordinary differential equations via numerical methods.
 * <p>
 *     This class and all functions within were added after I submitted my application to your institution. All code
 *     within this project will be submitted and graded for MTH49cr Numerical Differential Equations Modeling.
 * </p>
 * @author Keeler Spear
 * @version %I%, %G%
 * @since 1.0
 */
public class ODE {

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
    public static double[] euler (BiFunction ode, double t0, double y0, double t, double h) {
        double[] y = new double[1 + (int) ((t - t0)/h)];
        y[0] = y0;
        int i = 1;

        while (t0 < t) {
            y0 = y0 + h * (double) ode.apply(t0, y0);

            y[i] = y0;
            t0 += h;
            i++;
        }

        return y;
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
    public static double[] rk4 (BiFunction ode, double t0, double y0, double t, double h) {
        double[] y = new double[1 + (int) ((t - t0)/h)];
        y[0] = y0;
        double k1;
        double k2;
        double k3;
        double k4;
        int i = 1;

        while (t0 < t) {
            k1 = (double) ode.apply(t0, y0);
            k2 = (double) ode.apply(t0 + 0.5 * h, y0 + 0.5 * h * k1);
            k3 = (double) ode.apply(t0 + 0.5 * h, y0 + 0.5 * h * k2);
            k4 = (double) ode.apply(t0 + h, y0 + h * k3);

            y0 = y0 + h * ((1.0 / 6) * k1 + (1.0 / 3) * k2 + (1.0 / 3) * k3 + (1.0 / 6) * k4);

            y[i] = y0;
            t0 += h;
            i++;
        }

        return y;
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
    public static double[][] eulerSystem (TriFunction ode1, TriFunction ode2, double t0, double[] y0, double t, double h) {
        double[][] yVals = new double[2][1 + (int) ((t - t0)/h)];
        Matrix yi = new Matrix(y0);
        Matrix yP = new Matrix(y0);
        int i = 1;

        while (t0 < t) {
            //Finding y'(t)
            yP.setValue(1, 1, (double) ode1.apply(t0, yi.getValue(1,1), yi.getValue(2, 1)));
            yP.setValue(2, 1, (double) ode2.apply(t0, yi.getValue(1,1), yi.getValue(2, 1)));
            //finding y(t)
            yi = LinearAlgebra.addMatrices(yi, yP, h);

            yVals[0][i] = yi.getValue(1,1);
            yVals[1][i] = yi.getValue(2,1);
            t0 += h;
            i++;
        }

        return yVals;
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
    public static double[][] rk4System (TriFunction ode1, TriFunction ode2, double t0, double[] y0, double t, double h) {
        double[][] yVals = new double[2][1 + (int) ((t - t0)/h)];
        Matrix yi = new Matrix(y0);
        Matrix yP = new Matrix(y0);
        int i = 1;
        Matrix k1 = new Matrix(2, 1);
        Matrix k2 = new Matrix(2, 1);
        Matrix k3 = new Matrix(2, 1);
        Matrix k4 = new Matrix(2, 1);
        Matrix kTemp = new Matrix(2, 1); //Used as temporary storage for wi + h/2 ki
        Matrix wkSum;

        while (t0 < t) {
            //Finding the k's
            //k1:
            k1.setValue(1, 1, (double) ode1.apply(t0, yi.getValue(1,1), yi.getValue(2, 1)));
            k1.setValue(2, 1, (double) ode2.apply(t0, yi.getValue(1,1), yi.getValue(2, 1)));
            kTemp = LinearAlgebra.addMatrices(yi, k1, 0.5 * h);
            //k2:
            k2.setValue(1, 1, (double) ode1.apply(t0 + 0.5 * h, kTemp.getValue(1,1), kTemp.getValue(2, 1)));
            k2.setValue(2, 1, (double) ode2.apply(t0 + 0.5 * h, kTemp.getValue(1,1), kTemp.getValue(2, 1)));
            kTemp = LinearAlgebra.addMatrices(yi, k2, 0.5 * h);
            //k3:
            k3.setValue(1, 1, (double) ode1.apply(t0 + 0.5 * h, kTemp.getValue(1,1), kTemp.getValue(2, 1)));
            k3.setValue(2, 1, (double) ode2.apply(t0 + 0.5 * h, kTemp.getValue(1,1), kTemp.getValue(2, 1)));
            kTemp = LinearAlgebra.addMatrices(yi, k3, h);
            //k4:
            k4.setValue(1, 1, (double) ode1.apply(t0 + h, kTemp.getValue(1,1), kTemp.getValue(2, 1)));
            k4.setValue(2, 1, (double) ode2.apply(t0 + h, kTemp.getValue(1,1), kTemp.getValue(2, 1)));

            //Multiplying k's by their respective weights and adding them
            wkSum = LinearAlgebra.addMatrices(LinearAlgebra.scaleMatrix(k1, 1.0/6), k2, 1.0/3);
            wkSum = LinearAlgebra.addMatrices(wkSum, k3, 1.0/3);
            wkSum = LinearAlgebra.addMatrices(wkSum, k4, 1.0/6);

            yi = LinearAlgebra.addMatrices(yi, wkSum, h);

            yVals[0][i] = yi.getValue(1,1);
            yVals[1][i] = yi.getValue(2,1);
            t0 += h;
            i++;
        }

        return yVals;
    }

}
