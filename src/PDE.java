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
     * @return A numerical solution to the transport equation.
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
     * @param g The left boundary condition.
     * @param x0 The starting spacial position.
     * @param xMax The ending spacial position.
     * @param t0 The starting temporal position
     * @param tMax The ending temporal position.
     * @param h The space step to be used.
     * @param k The time step to be used.
     * @return A numerical solution to the transport equation.
     */
    public static Matrix solveTrEqFTBS(double c, Function<Double, Double> f, Function<Double, Double> g, double x0, double xMax, double t0, double tMax, double h, double k) {
        double lambda = c * k / h;
        int iterations = (int) Math.round((tMax - t0) / k);

        //Building v0
        Matrix v0 = LinearAlgebra.linSpace(x0, xMax, h);
        int size = v0.getRows();
        v0 = LinearAlgebra.applyFunction(v0, f);

        //Building A
        Matrix A = new Matrix(size, size);
        A.setValue(1,1, 1);
        for (int i = 2; i <= size; i++) {
            A.setValue(i, i, 1 - lambda);
        }
        for (int i = 2; i <= size; i++) {
            A.setValue(i, i - 1, lambda);
        }

        //Iterating the system in time
        return iterateSystem(A, v0, g, k, iterations);
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
     * @return A numerical solution to the transport equation.
     */
    public static Matrix solveTrEqFTCS(double c, Function<Double, Double> f, Function<Double, Double> g, double x0, double xMax, double t0, double tMax, double h, double k) {
        double lambda = c * k / (2 * h);
        int iterations = (int) Math.round((tMax - t0) / k);

        //Building v0
        Matrix v0 = LinearAlgebra.linSpace(x0, xMax, h);
        int size = v0.getRows();
        v0 = LinearAlgebra.applyFunction(v0, f);

        //Building A
        Matrix A = new Matrix(size, size);
        A.setValue(1,1, 1);
        for (int i = 2; i <= size; i++) {
            A.setValue(i, i, 1);
        }
        for (int i = 2; i < size; i++) {
            A.setValue(i, i - 1, lambda);
            A.setValue(i, i + 1, -lambda);
        }

        //Iterating the system in time
        return iterateSystem(A, v0, g, k, iterations);
    }

    /**
     * Solves the heat equation (ut = cuxx) with Dirichlet-Neumann boundary conditions numerically with FTCS.
     * This method assumes the Neumann BC is homogeneous.
     *
     * @param c The coefficient of the uxx term.
     * @param f The initial condition for the PDE.
     * @param g0 The left boundary condition.
     * @param x0 The starting spacial position.
     * @param xMax The ending spacial position.
     * @param t0 The starting temporal position
     * @param tMax The ending temporal position.
     * @param h The space step to be used.
     * @param k The time step to be used.
     * @return A numerical solution to the wave equation.
     */
    public static Matrix solveHeEqFTCS(double c, Function<Double, Double> f, Function<Double, Double> g0, double x0, double xMax, double t0, double tMax, double h, double k) {
        double lambda = c * k / (h * h);
        int iterations = (int) Math.round((tMax - t0) / k);

        //Building v0
        Matrix v0 = LinearAlgebra.linSpace(x0, xMax, h);
        int size = v0.getRows();
        v0 = LinearAlgebra.applyFunction(v0, f);

        //Building A
        Matrix A = new Matrix(size, size);
        A.setValue(1, 1, 1);
        A.setValue(size, size, 1);
        for (int i = 2; i < size; i++) {
            A.setValue(i, i, 1 - 2 * lambda);
            A.setValue(i, i - 1, lambda);
            A.setValue(i, i + 1, lambda);
        }
        //Neumann BC
        A.setValue(size, size, 1 - 2 * lambda);
        A.setValue(size, size - 1, 2 * lambda);

        //Iterating the system in time
        return iterateSystem(A, v0, g0, k, iterations);
    }

    /**
     * Solves the heat equation (ut = cuxx) with Dirichlet-Neumann boundary conditions numerically with BTCS.
     * This method assumes the Neumann BC is homogeneous.
     *
     * @param c The coefficient of the uxx term.
     * @param f The initial condition for the PDE.
     * @param g0 The left boundary condition.
     * @param x0 The starting spacial position.
     * @param xMax The ending spacial position.
     * @param t0 The starting temporal position
     * @param tMax The ending temporal position.
     * @param h The space step to be used.
     * @param k The time step to be used.
     * @return A numerical solution to the wave equation.
     */
    public static Matrix solveHeEqBTCS(double c, Function<Double, Double> f, Function<Double, Double> g0, double x0, double xMax, double t0, double tMax, double h, double k) {
        double lambda = c * k / (h * h);
        int iterations = (int) Math.round((tMax - t0) / k);

        //Building v0
        Matrix v0 = LinearAlgebra.linSpace(x0, xMax, h);
        int size = v0.getRows();
        v0 = LinearAlgebra.applyFunction(v0, f);

        //Building A
        Matrix A = new Matrix(size, size);
        A.setValue(1, 1, 1);
        A.setValue(size, size, 1);
        for (int i = 2; i < size; i++) {
            A.setValue(i, i, 1 + 2 * lambda);
            A.setValue(i, i - 1, -lambda);
            A.setValue(i, i + 1, -lambda);
        }
        //Neumann BC
        A.setValue(size, size, 1 + 2 * lambda);
        A.setValue(size, size - 1, -2 * lambda);

        //Iterating the system in time
        return iterateAndSolveSystem(A, v0, g0, k, iterations);
    }

    /**
     * Solves the heat equation (ut = cuxx) with Dirichlet-Neumann boundary conditions numerically with Crank-Nicolson.
     * This method assumes the Neumann BC is homogeneous.
     *
     * @param c The coefficient of the uxx term.
     * @param f The initial condition for the PDE.
     * @param g0 The left boundary condition.
     * @param x0 The starting spacial position.
     * @param xMax The ending spacial position.
     * @param t0 The starting temporal position
     * @param tMax The ending temporal position.
     * @param h The space step to be used.
     * @param k The time step to be used.
     * @return A numerical solution to the wave equation.
     */
    public static Matrix solveHeEqCN(double c, Function<Double, Double> f, Function<Double, Double> g0, double x0, double xMax, double t0, double tMax, double h, double k) {
        double lambda = c * k / (h * h);
        int iterations = (int) Math.round((tMax - t0) / k);

        //Building v0
        Matrix v0 = LinearAlgebra.linSpace(x0, xMax, h);
        int size = v0.getRows();
        v0 = LinearAlgebra.applyFunction(v0, f);

        //Building A and B
        Matrix A = new Matrix(size, size);
        Matrix B = new Matrix(size, size);
        A.setValue(1, 1, 1);
        A.setValue(size, size, 1);
        B.setValue(1, 1, 1);
        B.setValue(size, size, 1);
        for (int i = 2; i < size; i++) {
            A.setValue(i, i, 1 + lambda);
            A.setValue(i, i - 1, -lambda / 2);
            A.setValue(i, i + 1, -lambda / 2);

            B.setValue(i, i, 1 - lambda);
            B.setValue(i, i - 1, lambda / 2);
            B.setValue(i, i + 1, lambda / 2);
        }

        //Neumann BC
        A.setValue(size, size, 1 + lambda);
        A.setValue(size, size - 1, -lambda);
        B.setValue(size, size, 1 - lambda);
        B.setValue(size, size - 1, lambda);


        //Iterating the system in time
        return iterateAndSolveSystem(A, B, v0, g0, k, iterations);
    }

    /**
     * Solves the heat equation (ut = cuxx) with Dirichlet-Neumann boundary conditions numerically with Leapfrog.
     * This method assumes the Neumann BC is homogeneous.
     *
     * @param c The coefficient of the uxx term.
     * @param f The initial condition for the PDE.
     * @param g0 The left boundary condition.
     * @param x0 The starting spacial position.
     * @param xMax The ending spacial position.
     * @param t0 The starting temporal position
     * @param tMax The ending temporal position.
     * @param h The space step to be used.
     * @param k The time step to be used.
     * @return A numerical solution to the wave equation.
     */
    public static Matrix solveHeEqLF(double c, Function<Double, Double> f, Function<Double, Double> g0, double x0, double xMax, double t0, double tMax, double h, double k) {
        double lambda = c * k / (h * h);
        int iterations = (int) Math.round((tMax - t0) / k);
        double vmn;
        double vmnm1;
        double vmnp1;

        //Building v0
        Matrix v0 = LinearAlgebra.linSpace(x0, xMax, h);
        int size = v0.getRows();
        v0 = LinearAlgebra.applyFunction(v0, f);

        //Building V1 via Euler's Method
        Matrix v1 = new Matrix(v0.getRows(), 1);
        v1.setValue(1, 1, g0.apply(k));
        for (int i = 2; i < v1.getRows(); i++) {
            vmnm1 = v0.getValue(i - 1, 1);
            vmn = v0.getValue(i, 1);
            vmnp1 = v0.getValue(i + 1, 1);
            v1.setValue(i, 1, vmn + lambda * (vmnp1 - 2 * vmn + vmnm1));
        }
        vmnm1 = v0.getValue(v0.getRows() - 1, 1);
        vmn = v0.getValue(v0.getRows(), 1);
        v1.setValue(v1.getRows(), 1, vmn + lambda * (-2 * vmn + 2 * vmnm1));

        //Building A
        Matrix A = new Matrix(size, size);
        A.setValue(1, 1, 1);
        A.setValue(size, size, 1);
        for (int i = 2; i < size; i++) {
            A.setValue(i, i, -4 *  lambda);
            A.setValue(i, i - 1, 2 * lambda);
            A.setValue(i, i + 1, 2 * lambda);
        }
        //Neumann BC
        A.setValue(size, size, -4 *  lambda);
        A.setValue(size, size - 1, 4 * lambda);

        //Iterating the system in time
        return iterateSystemP(A, v0, v1, g0, k, iterations);
    }

    /**
     * Solves the wave equation (utt = cuxx) with Dirichlet boundary conditions numerically with CTCS.
     *
     * @param c The coefficient of the uxx term.
     * @param f The initial condition for the PDE.
     * @param ft The initial condition's time derivative.
     * @param g0 The left boundary condition.
     * @param gl The right boundary condition.
     * @param x0 The starting spacial position.
     * @param xMax The ending spacial position.
     * @param t0 The starting temporal position
     * @param tMax The ending temporal position.
     * @param h The space step to be used.
     * @param k The time step to be used.
     * @return A numerical solution to the wave equation.
     */
    public static Matrix solveWaEqCTCS(double c, Function<Double, Double> f, Function<Double, Double> ft, Function<Double, Double> g0, Function<Double, Double> gl, double x0, double xMax, double t0, double tMax, double h, double k) {
        double lambda = c * k * k / (h * h);
        int iterations = (int) Math.round((tMax - t0) / k);
        double vmn;
        double vmnm1;
        double vmnp1;

        //Building v0
        Matrix v0 = LinearAlgebra.linSpace(x0, xMax, h);
        int size = v0.getRows();
        v0 = LinearAlgebra.applyFunction(v0, f);

        //Building V1 via a Taylor Series expansion
        Matrix v1 = new Matrix(v0.getRows(), 1);
        v1.setValue(1, 1, g0.apply(k));
        v1.setValue(v1.getRows(), 1, gl.apply(k));
        for (int i = 2; i < v1.getRows(); i++) {
            vmnm1 = v0.getValue(i - 1, 1);
            vmn = v0.getValue(i, 1);
            vmnp1 = v0.getValue(i + 1, 1);
            v1.setValue(i, 1, (1 - lambda) * vmn + 0.5 * lambda * (vmnp1+ vmnm1) + k * ft.apply(h * i));
        }

        //Building A
        Matrix A = new Matrix(size, size);
        for (int i = 2; i < size; i++) {
            A.setValue(i, i, 2 * (1 - lambda));
            A.setValue(i, i - 1, lambda);
            A.setValue(i, i + 1, lambda);
        }
        A.setValue(1, 1, 1);
        A.setValue(size, size, 1);

        //Iterating the system in time
        Matrix sol = LinearAlgebra.transpose(v0);
        sol.addRowBottom(v1.getMatrix());
        Matrix temp;

        return iterateSystem(A, v0, v1, g0, gl, k, iterations);
    }

    /**
     * Solves Laplace's Equation (uxx + uyy = 0) with Dirichlet boundary conditions numerically with CTCS.
     *
     * @param f1 The initial condition u(x, y0) for the PDE.
     * @param f2 The initial condition u(x, yMax) for the PDE.
     * @param g1 The initial condition u(x0, y) for the PDE.
     * @param g2 The initial condition u(xMax, y) for the PDE.
     * @param x0 The starting x spacial position.
     * @param xMax The ending x spacial position.
     * @param y0 The starting y spacial position.
     * @param yMax The ending y spacial position.
     * @param h The x space step to be used.
     * @param k The y space step to be used.
     * @return A numerical solution to Laplace's Equation.
     */
    public static Matrix solveLaEq(Function<Double, Double> f1, Function<Double, Double> f2, Function<Double, Double> g1, Function<Double, Double> g2, double x0, double xMax, double y0, double yMax, double h, double k) {
//        double lambda = (k * k) / (h * h);
//        int iterations = (int) Math.round((tMax - t0) / k);
//        double vmn;
//        double vmnm1;
//        double vmnp1;
//
//        //Building v0
//        Matrix v0 = LinearAlgebra.linSpace(x0, xMax, h);
//        int size = v0.getRows();
//        v0 = LinearAlgebra.applyFunction(v0, f);
//
//        //Building V1 via a Taylor Series expansion
//        Matrix v1 = new Matrix(v0.getRows(), 1);
//        v1.setValue(1, 1, g0.apply(k));
//        v1.setValue(v1.getRows(), 1, gl.apply(k));
//        for (int i = 2; i < v1.getRows(); i++) {
//            vmnm1 = v0.getValue(i - 1, 1);
//            vmn = v0.getValue(i, 1);
//            vmnp1 = v0.getValue(i + 1, 1);
//            v1.setValue(i, 1, (1 - lambda) * vmn + 0.5 * lambda * (vmnp1+ vmnm1) + k * ft.apply(h * i));
//        }
//
//        //Building A
//        Matrix A = new Matrix(size, size);
//        for (int i = 2; i < size; i++) {
//            A.setValue(i, i, 2 * (1 - lambda));
//            A.setValue(i, i - 1, lambda);
//            A.setValue(i, i + 1, lambda);
//        }
//        A.setValue(1, 1, 1);
//        A.setValue(size, size, 1);
//
//        //Iterating the system in time
//        Matrix sol = LinearAlgebra.transpose(v0);
//        sol.addRowBottom(v1.getMatrix());
//        Matrix temp;
//
//        return iterateSystem(A, v0, v1, g0, gl, k, iterations);
    }



    private static Matrix iterateSystem(Matrix A, Matrix v, int n) {
        Matrix sol = LinearAlgebra.transpose(v);

        for (int i = 0; i < n; i++) {
            v = LinearAlgebra.multiplyMatrices(A, v);
            sol.addRowBottom(v.getMatrix());
        }

        return sol;
    }

    private static Matrix iterateSystem(Matrix A, Matrix v, Function<Double, Double> g0, double k, int n) {
        Matrix sol = LinearAlgebra.transpose(v);

        for (int i = 2; i < n + 2; i++) {
            v = LinearAlgebra.multiplyMatrices(A, v);
            //Enforcing BC
            v.setValue(1,1, g0.apply(i * k));
            sol.addRowBottom(v.getMatrix());
        }

        return sol;
    }

    private static Matrix iterateSystem(Matrix A, Matrix v, Function<Double, Double> g0, Function<Double, Double> gl, double k, int n) {
        Matrix sol = LinearAlgebra.transpose(v);

        for (int i = 2; i < n + 2; i++) {
            v = LinearAlgebra.multiplyMatrices(A, v);
            //Enforcing BC
            v.setValue(1,1, g0.apply(i * k));
            v.setValue(v.getRows(), 1, gl.apply(i * k));
            sol.addRowBottom(v.getMatrix());
        }

        return sol;
    }

    private static Matrix iterateAndSolveSystem(Matrix A, Matrix v, Function<Double, Double> g0, double k, int n) {
        Matrix sol = LinearAlgebra.transpose(v);

        for (int i = 2; i < n + 2; i++) {
            v = LinearAlgebra.RREFSolve(LinearAlgebra.augmentMatrix(A, v));
            //Enforcing BC
            v.setValue(1,1, g0.apply(i * k));
            sol.addRowBottom(v.getMatrix());
        }

        return sol;
    }

    private static Matrix iterateAndSolveSystem(Matrix A, Matrix B, Matrix v, Function<Double, Double> g0, double k, int n) {
        Matrix sol = LinearAlgebra.transpose(v);

        for (int i = 2; i < n + 2; i++) {
            v = LinearAlgebra.multiplyMatrices(B, v);
            v = LinearAlgebra.RREFSolve(LinearAlgebra.augmentMatrix(A, v));
            //Enforcing BC
            v.setValue(1,1, g0.apply(i * k));
            sol.addRowBottom(v.getMatrix());
        }

        return sol;
    }

    private static Matrix iterateSystem(Matrix A, Matrix v0, Matrix v1, Function<Double, Double> g0, Function<Double, Double> gl, double k, int n) {
        Matrix sol = LinearAlgebra.transpose(v0);
        sol.addRowBottom(v1.getMatrix());
        Matrix temp;

        for (int i = 2; i < n + 1; i++) {
            temp = v1;
            v1 = LinearAlgebra.multiplyMatrices(A, v1);
            v1 = LinearAlgebra.subtractMatrices(v1, v0);
            //enforcing BC
            v1.setValue(1,1, g0.apply(i * k));
            v1.setValue(v1.getRows(), 1, gl.apply(i * k));
            v0 = temp;
            sol.addRowBottom(v1.getMatrix());
        }

        return sol;
    }

    private static Matrix iterateSystem(Matrix A, Matrix v0, Matrix v1, Function<Double, Double> g0, double k, int n) {
        Matrix sol = LinearAlgebra.transpose(v0);
        sol.addRowBottom(v1.getMatrix());
        Matrix temp;

        for (int i = 2; i < n + 1; i++) {
            temp = v1;
            v1 = LinearAlgebra.multiplyMatrices(A, v1);
            v1 = LinearAlgebra.subtractMatrices(v1, v0);
            //enforcing BC
            v1.setValue(1,1, g0.apply(i * k));
            v0 = temp;
            sol.addRowBottom(v1.getMatrix());
        }

        return sol;
    }

    private static Matrix iterateSystemP(Matrix A, Matrix v0, Matrix v1, Function<Double, Double> g0, double k, int n) {
        Matrix sol = LinearAlgebra.transpose(v0);
        sol.addRowBottom(v1.getMatrix());
        Matrix temp;

        for (int i = 2; i < n + 1; i++) {
            temp = v1;
            v1 = LinearAlgebra.multiplyMatrices(A, v1);
            v1 = LinearAlgebra.addMatrices(v1, v0);
            //enforcing BC
            v1.setValue(1,1, g0.apply(i * k));
            v0 = temp;
            sol.addRowBottom(v1.getMatrix());
        }

        return sol;
    }

    public static Matrix build2DFunction(Matrix x, Matrix y, BiFunction<Double, Double, Double> f) {
        Matrix sol = new Matrix(y.getRows(), x.getRows());
        for (int i = 1; i <= sol.getRows(); i++) {
            for (int j = 1; j <= sol.getCols(); j++) {
                sol.setValue(i, j, f.apply(x.getValue(j, 1), y.getValue(i, 1)));
            }
        }

        return sol;
    }

}
