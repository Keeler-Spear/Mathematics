import com.sun.security.jgss.GSSUtil;

import java.util.ArrayList;
import java.util.Arrays;
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

    private static final double BASE_VAL = 9999; //The present value matrices are filled with
    private static final int MAX_ITERATIONS = 10000;
    private static final double EPSILON = 1.e-15;
    private static final double TOL = 0.1;

    /**
     * Numerically solves the first-order ordinary differential equation over the provided interval using the
     * Runge–Kutta 4 and Predictor-Corrector 4 methods with O(h^4) global error, where h is the step size of t.
     *
     * @param ode The first-order ordinary differential equation f(t, y) to be solved.
     * @param t0 The left endpoint of the interval.
     * @param y0 The y value corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param h The step size of t.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     */
    public static Matrix solveIVP(BiFunction<Double, Double, Double> ode, double t0, double y0, double t, double h) {
        return pc(ode, t0, y0, t, h);
    }

    /**
     * Numerically solves the second-order ordinary differential equation with the IVP provided over the provided interval using the
     * Runge–Kutta 4 method with O(h^4) global error, where h is the step size of t.
     *
     * @param ode The second-order ordinary differential equation f(t, y, y') to be solved.
     * @param t0 The left endpoint of the interval.
     * @param y0 The y value corresponding to t0.
     * @param yp0 The y' value corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param h The step size of t.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     */
    public static Matrix solveIVP(TriFunction<Double, Double, Double, Double> ode, double t0, double y0, double yp0, double t, double h) {
        TriFunction<Double, Double, Double, Double>[] odes = decomposeODE(ode);
        return rk4System2(odes[0], odes[1], t0, y0, yp0, t, h);
    }

    /**
     * Numerically solves the second-order ordinary differential equation with the BVP provided over the provided interval using the
     * Runge–Kutta 4 Linear Shooting method with O(h^4) global error, where h is the step size of t.
     *
     * @param ode The second-order ordinary differential equation f(t, y, y') to be solved.
     * @param t0 The left endpoint of the interval.
     * @param y0 The y value corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param y1 The y value corresponding to t.
     * @param h The step size of t.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     */
    public static Matrix solveBVP (TriFunction<Double, Double, Double, Double> ode, double t0, double y0, double t, double y1, double h) {
        return linearShooting(ode, t0, y0, t, y1, h);
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
     * @param y0 The y value corresponding to t0.
     * @param yp0 The y' value corresponding to t0.
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
        t0 += h;
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
     * @param y0 The y value corresponding to t0.
     * @param yp0 The y' value corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param h The step size of t.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     */
    public static Matrix rk4System2 (TriFunction<Double, Double, Double, Double> ode1, TriFunction<Double, Double, Double, Double> ode2, double t0, double y0, double yp0, double t, double h) {
        Matrix yVals = LinearAlgebra.constantMatrix(generateYLength(t0, t, h) + 50, 2, BASE_VAL);
        Matrix yi = new Matrix(new double[] {y0, yp0});
        yVals.setValue(1, 1, y0);
        yVals.setValue(1, 2, yp0);
        int i = 2;
        t0 += h;
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
     * Numerically solves the system of first-order ordinary differential equations over the provided interval using the
     * Runge–Kutta 4 method with O(h^4) global error, where h is the step size of t.
     *
     * @param system The system of first-order ordinary differential equations.
     * @param t0 The left endpoint of the interval.
     * @param y0 The y values corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param h The step size of t.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     * @throws IllegalArgumentException If each ODE does not have an initial condition.
     */
    public static Matrix rk4System (NFunction<Double>[] system, double t0, double[] y0, double t, double h) {
        if (system.length != y0.length) {
            throw new IllegalArgumentException("Each ODE must have an initial condition!");
        }

        //Changed rows by removing buffer and also cut the trim function from the return
        Matrix yVals = LinearAlgebra.constantMatrix(generateYLength(t0, t, h)  + 5, system.length, BASE_VAL);

        Matrix yi = new Matrix(y0.length, 1);
        for (int i = 0; i < y0.length; i++) {
            yi.setValue(i + 1, 1, y0[i]);
            yVals.setValue(1, i + 1, y0[i]);
        }

        int j = 2;
        t0 += h;

        Matrix k1 = new Matrix(system.length, 1);
        Matrix k2 = new Matrix(system.length, 1);
        Matrix k3 = new Matrix(system.length, 1);
        Matrix k4 = new Matrix(system.length, 1);
        Matrix kTemp; //Used as temporary storage for wi + h/2 ki
        Matrix wkSum;
        Double[] perams = new Double[y0.length + 1]; //Parameter values to pass to the ODE

        while (t0 < t) {
            //Finding the k's
            //k1:
            for (int i = 1; i <= system.length; i++) {
                perams[0] = t0;
                for (int k = 1; k <= y0.length; k++) {
                    perams[k] = yi.getValue(k, 1);
                }
                k1.setValue(i, 1, system[i - 1].apply(perams));
            }
            kTemp = LinearAlgebra.addMatrices(yi, k1, 0.5 * h);

            //k2:
            for (int i = 1; i <= system.length; i++) {
                perams[0] = t0 + 0.5 * h;
                for (int k = 1; k <= y0.length; k++) {
                    perams[k] = kTemp.getValue(k, 1);
                }
                k2.setValue(i, 1, system[i - 1].apply(perams));
            }
            kTemp = LinearAlgebra.addMatrices(yi, k2, 0.5 * h);

            //k3:
            for (int i = 1; i <= system.length; i++) {
                perams[0] = t0 + 0.5 * h;
                for (int k = 1; k <= y0.length; k++) {
                    perams[k] = kTemp.getValue(k, 1);
                }
                k3.setValue(i, 1, system[i - 1].apply(perams));
            }
            kTemp = LinearAlgebra.addMatrices(yi, k3, h);

            //k4:
            for (int i = 1; i <= system.length; i++) {
                perams[0] = t0 + h;
                for (int k = 1; k <= y0.length; k++) {
                    perams[k] = kTemp.getValue(k, 1);
                }
                k4.setValue(i, 1, system[i - 1].apply(perams));
            }

            //Multiplying k's by their respective weights and adding them
            wkSum = LinearAlgebra.addMatrices(LinearAlgebra.scaleMatrix(k1, 1.0/6), k2, 1.0/3);
            wkSum = LinearAlgebra.addMatrices(wkSum, k3, 1.0/3);
            wkSum = LinearAlgebra.addMatrices(wkSum, k4, 1.0/6);

            yi = LinearAlgebra.addMatrices(yi, wkSum, h);

            for (int i = 1; i <= system.length; i++) {
                yVals.setValue(j, i, yi.getValue(i,1));
            }

            t0 += h;
            j++;
        }

        return trim(yVals);
    }

    /**
     * Numerically solves the system of first-order ordinary differential equations over the provided interval using the
     * adaptive Runge–Kutta 4 method.
     * <p>
     *     This code was copied and translated from https://github.com/AlejGarcia/NM4P/blob/master/Python/nm4p/rka.py.
     * </p>
     *
     * @param system The system of first-order ordinary differential equations.
     * @param t0 The left endpoint of the interval.
     * @param y0 The y values corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param h The step size of t.
     * @param s1 The first safety factor.
     * @param s2 The second safety factor.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     * @throws IllegalArgumentException If each ODE does not have an initial condition.
     * @throws IllegalArgumentException If the first safety factor is greater than or equal to 1.
     * @throws IllegalArgumentException If the second safety factor is equal to or greater than 1.
     */
    public static Matrix adaptiveRK4System (NFunction<Double>[] system, double t0, double[] y0, double t, double h, double s1, double s2) {
        if (system.length != y0.length) {
            throw new IllegalArgumentException("Each ODE must have an initial condition!");
        }

        if (s1 >= 1 || s1 <= 0) {
            throw new IllegalArgumentException("The first safety factor must be less than 1!");
        }

        if (s2 <= 1) {
            throw new IllegalArgumentException("The second safety factor must be greater than 1!");
        }

        Matrix vals = new Matrix(y0);
        ArrayList<Double> time = new ArrayList<>();
        time.add(t0);

        double tSave;
        double[] ySave;
        Matrix ySmallMat = new Matrix(2, 1);
        double ySmall;
        Matrix yBigMat = new Matrix(2, 1);
        double yBig;
        double yDiff;
        double scale;
        double halfH;
        double oldH;
        double errorRatio = 0.0;
        int j = 0;

        while (t0 < t && j < MAX_ITERATIONS) {


            for (int i = 0; i < 100; i++) {
                //Two small steps
                tSave = t0;
                ySave = y0.clone();
                halfH = 0.5 * h;
                ySmallMat = rk4System(system, tSave, y0, tSave + h, halfH); //Will take 2 half steps
                ySmall = Math.sqrt(Math.pow(ySmallMat.getValue(ySmallMat.getRows(), 1), 2) + Math.pow(ySmallMat.getValue(ySmallMat.getRows(), 2), 2));
                t0 = tSave + halfH;

                //One big step
                yBigMat = rk4System(system, tSave, ySave, tSave + h, h); //Will take one big step
                yBig = Math.sqrt(Math.pow(yBigMat.getValue(yBigMat.getRows(), 1), 2) + Math.pow(yBigMat.getValue(yBigMat.getRows(), 2), 2));

                // Calculate the truncation error
                scale = TOL * (Math.abs(ySmall) + Math.abs(yBig)) / 2;
                yDiff = ySmall - yBig;
                errorRatio = Math.abs(yDiff) / (scale + EPSILON);

                // Adjust the step size
                oldH = h;
                if (errorRatio > 1) {
                    // Reduce step size and retry
                    h = s1 * h * Math.pow(errorRatio, -0.20);
                    h = Math.max(h, oldH / s2);
                }
                else {
                    // Accept the step and increase step size for the next iteration
                    t0 = tSave + h;
                    time.add(t0);

                    vals = LinearAlgebra.augmentMatrix(vals, LinearAlgebra.vectorFromRow(yBigMat, yBigMat.getRows()));

                    for (i = 0; i < y0.length; i++) {
                        y0[i] = yBigMat.getValue(yBigMat.getRows(), i + 1);
                    }

                    h = s1 * h * Math.pow(errorRatio, -0.20);
                    h = Math.min(h, oldH * s2);
                    break;
                }

            }

            //After the for loop, we have our h.
            if (errorRatio < 1) {
                vals = LinearAlgebra.augmentMatrix(vals, LinearAlgebra.vectorFromRow(yBigMat, yBigMat.getRows()));
            }

            else {
                vals = LinearAlgebra.augmentMatrix(vals, LinearAlgebra.vectorFromRow(ySmallMat, ySmallMat.getRows()));

            }

            for (int i = 1; i <= ySmallMat.getCols(); i++) {
                y0[i - 1] = ySmallMat.getValue(ySmallMat.getRows(), i);
            }
            t0 += h;
            time.add(t0);
            j++;
        }

        //Since the time variable does not increase linearly, time will be returned along with the state vector
        double[] times = new double[time.size()];
        for (int i = 0; i < time.size(); i++) {
            times[i] = time.get(i);
        }

        vals = LinearAlgebra.transpose(vals);
        vals.addColRight(times);


        return vals;
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
     * Runge–Kutta 4 and Predictor-Corrector 4 methods with O(h^4) global error, where h is the step size of t.
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

    /**
     * Numerically solves the second-order ordinary differential equation with the BVP provided over the provided interval using the
     * Runge–Kutta 4 Linear Shooting method with O(h^4) global error, where h is the step size of t.
     *
     * @param ode The second-order ordinary differential equation f(t, y, y') to be solved.
     * @param t0 The left endpoint of the interval.
     * @param y0 The y value corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param y1 The y value corresponding to t.
     * @param h The step size of t.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     */
    public static Matrix linearShooting (TriFunction<Double, Double, Double, Double> ode, double t0, double y0, double t, double y1, double h) {
        Matrix yp = solveIVP(ode, t0, y0, 0, t, h);

        TriFunction<Double, Double, Double, Double> homoODE = (x, y, yprime) -> ode.apply(x, y, yprime) - ode.apply(x, 0.0, 0.0);

        Matrix yc = solveIVP(homoODE, t0, 0, 1, t, h);

        double ypb  = yp.getValue(yp.getRows(), 1);
        double ycb  = yc.getValue(yc.getRows(), 1);

        return LinearAlgebra.addMatrices(yp, yc, (y1 - ypb) / ycb);
    }

    /**
     * Numerically solves the second-order ordinary differential equation with the BVP provided over the provided interval using the
     * Finite Differences method with O(h^q) global error, where h is the step size of t.
     *
     * @param ode The second-order ordinary differential equation f(t, y, y') to be solved.
     * @param t0 The left endpoint of the interval.
     * @param y0 The y value corresponding to t0.
     * @param t The right endpoint of the interval.
     * @param y1 The y value corresponding to t.
     * @param h The step size of t.
     * @return The numerical approximation of the ordinary differential equation over the interval [t0, t].
     */
    public static Matrix finiteDifferences (TriFunction<Double, Double, Double, Double> ode, double t0, double y0, double t, double y1, double h) {
        int n = generateYLength(t0, t, h);
        Matrix A = new Matrix(n, n);
        Matrix b = new Matrix(n, 1);
        double ti = t0;

        Function<Double, Double> rx = (x) -> ode.apply(x, 0.0, 0.0);
        Function<Double, Double> fx = (x) -> 2.0 + h * h * (ode.apply(x, 1.0, 0.0) - rx.apply(x)); //g(x)
        Function<Double, Double> gx = (x) -> (h / 2.0) * (ode.apply(x, 0.0, 1.0) - rx.apply(x)); //h(x)
        Function<Double, Double> bx = (x) -> h * h * rx.apply(x) / fx.apply(x);

        //Building A
        A.setValue(1, 1, 1.0);
        for (int i = 2; i < n; i++) {
            ti += h;
            A.setValue(i, i - 1, -(gx.apply(ti) + 1) / fx.apply(ti));
            A.setValue(i, i, 1.0);
            A.setValue(i, i + 1, (gx.apply(ti) - 1) / fx.apply(ti));
        }

        A.setValue(n, n, 1.0);

        //Building B
        ti = t0;
        b.setValue(1, 1, y0);
        for (int i = 2; i < n; i++) {
            ti += h;
            b.setValue(i, 1, -bx.apply(ti));
        }
        b.setValue(n, 1, y1);

        return LinearAlgebra.RREFSolve(LinearAlgebra.augmentMatrix(A, b));
    }

    //Input should be y'' = a(x)y' + b(x)y + c. Params: (x, y0, y1);
    private static TriFunction[] decomposeODE(TriFunction<Double, Double, Double, Double> ode) {
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
        int n = 1 + (int) Math.round(((t - t0) / h));

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
