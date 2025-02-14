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
     * Numerically solves the ordinary differential equation over the provided interval using a Runge–Kutta 4 method
     * with O(h^4) global error, where h is the step size of t.
     *
     * @param ode The ordinary differential equation f(t, y) to be solved.
     * @param t0 The left endpoint of the interval.
     * @param y0 The y value corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param h The step size of t.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     */
    public static double[] rk4 (BiFunction ode, double t0, double y0, double t, double h) {
        double[] y = new double[1 + (int) ((t - t0)/h)];
        y[0] = y0;
        int i = 1;
        double k1;
        double k2;
        double k3;
        double k4;

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
     * Numerically solves the ordinary differential equation over the provided interval using Euler's method
     * with O(h) global error, where h is the step size of t.
     *
     * @param ode The ordinary differential equation f(t, y) to be solved.
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

}
