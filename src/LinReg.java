import java.util.function.Function;

/**
 * A static class that performs the operations essential to machine learning.
 * <p>
 *     This class and all functions within were added after I submitted my application to your institution.
 * </p>
 *
 * @author Keeler Spear
 * @version %I%, %G%
 * @since 1.0
 */
public class LinReg {
    final static double TOL = 0.01;
    final static int MAX_ITERATIONS = 100000;
    final static double LR = 0.001;
    final static double RAND_BOUND = 10;

    /**
     * Calculates the weights for a polynomial regression model based on the data set provided using gradient descent.
     * <p>
     * The weights will form the following polynomial function:
     * w + wx + wy + wz + ... + wx^2 + wx^3 + wx^n + ... + wy^2 + wy^3 + ... + wy^n + wz^2 + wz^3 + ... + wz^n + ...
     * </p>
     *
     * @param x A matrix of data parameters.
     * @param y A vector of data labels.
     * @param n The order of the polynomial.
     * @return A matrix with the weights for a polynomial regression model based on the data set provided.
     * @throws IllegalArgumentException If the data does not have one sample for each label.
     * @throws IllegalArgumentException If n is less than one.
     * @see #gradDes(Matrix, Matrix, Matrix, int, double, Function[])
     */
    public static Matrix polyGradDes(Matrix x, Matrix y, int n) {
        if (x.getRows() != y.getRows()) {
            throw new IllegalArgumentException("The data does not have one sample for each label!");
        }

        if (n < 1) {
            throw new IllegalArgumentException("The polynomial's order must be equal to or greater than one!");
        }

        Matrix w = LinearAlgebra.randMatrix(n * x.getCols() + 1, 1, -RAND_BOUND, RAND_BOUND);
        Function[] fncs = BasisFunctions.polynomials(n);
        double a;

        if (n <= 2) {
            a = LR;
        }
        else {
            a = Math.pow(10, 1 - n);
        }

        return polyGradDes(x, y, w, n, a);
    }

    /**
     * Calculates the weights for a polynomial regression model based on the data set provided using gradient descent.
     * <p>
     * The weights will form the following polynomial function:
     * w + wx + wy + wz + ... + wx^2 + wx^3 + wx^n + ... + wy^2 + wy^3 + ... + wy^n + wz^2 + wz^3 + ... + wz^n + ...
     * </p>
     *
     * @param x A matrix of data parameters.
     * @param y A vector of data labels.
     * @param w A vector of initial weight guesses.
     * @param n The order of the polynomial.
     * @param a The learning rate of the model.
     * @return A matrix with the weights for a polynomial regression model based on the data set provided.
     * @throws IllegalArgumentException If the data does not have one sample for each label.
     * @throws IllegalArgumentException If n is less than one.
     * @see #gradDes(Matrix, Matrix, Matrix, int, double, Function[])
     */
    public static Matrix polyGradDes(Matrix x, Matrix y, Matrix w, int n, double a) {
        if (x.getRows() != y.getRows()) {
            throw new IllegalArgumentException("The data does not have one sample for each label!");
        }

        if (n < 1) {
            throw new IllegalArgumentException("The polynomial's order must be equal to or greater than one!");
        }

        Function[] fncs = BasisFunctions.polynomials(n);

        return gradDes(x, y, w, n, a, fncs);
    }

    //ToDo: Make the model only in the terms of fncs[]
    /**
     * Calculates the weights for a regression model based on the data set provided and the basis functions provided
     * using gradient descent.
     * <p>
     * The weights will form the following polynomial function:
     * w + wx + wy + wz + ... + wx^2 + wx^3 + wx^n + ... + wy^2 + wy^3 + ... + wy^n + wz^2 + wz^3 + ... + wz^n + ...
     * </p>
     *
     * @param xOrg A matrix of data parameters.
     * @param y A vector of data labels.
     * @param w A vector of initial weight guesses.
     * @param n The order of the polynomial.
     * @param a The learning rate of the model.
     * @param fncs The set of basis functions the data will be modeled with.
     * @return A matrix with the weights for a polynomial regression model based on the data set provided.
     * @throws IllegalArgumentException If the data does not have one sample for each label.
     * @throws IllegalArgumentException If there is not one weight for each parameter.
     * @throws IllegalArgumentException If n is less than 1.
     * @throws IllegalArgumentException If the number of basis functions is not equal to n + 1.
     */
    public static Matrix gradDes(Matrix xOrg, Matrix y, Matrix w, int n, double a, Function[] fncs) {
        if (xOrg.getRows() != y.getRows()) {
            throw new IllegalArgumentException("The data does not have one sample for each label!");
        }

        if (w.getRows() != n * xOrg.getCols() + 1) {
            throw new IllegalArgumentException("There must be one weight for each parameter!");
        }

        if (n < 1) {
            throw new IllegalArgumentException("The polynomial's order must be equal to or greater than 1!");
        }

        if (fncs.length != n + 1) {
            throw new IllegalArgumentException("There must be exactly " + (n + 1) + " basis functions!");
        }

        Matrix x = new Matrix (xOrg.getRows(), 1);

        for (int i = 1; i <= x.getRows(); i++) {
            x.setValue(i, 1, (Double) fncs[0].apply(xOrg.getValue(i, 1)));
        }

        int numVars = xOrg.getCols();
        Matrix xn;


        //Adding the basis functions to the matrix
        for (int i = 1; i <= numVars; i++) {//For each variable
            for (int j = 1; j <= n; j++) {//For each function order
                xn = LinearAlgebra.vectorFromColumn(xOrg, i);
                xn = LinearAlgebra.applyFunction(xn, fncs[j]);
                x.addColRight(xn.getMatrix());
            }
        }

        //Should I always include a bias term, even when modeling data without polynomials?
//        Matrix ones = LinearAlgebra.constantMatrix(x.getRows(), 1, 1.0);
//        x.addColLeft(ones.getMatrix());

        Matrix dw = LinearAlgebra.constantMatrix(w.getRows(), 1, 5.0);

        int i = 0;
        while (LinearAlgebra.l2Norm(dw) > TOL && i < MAX_ITERATIONS) {
            dw = findDw(x, y, w);
            w = LinearAlgebra.addMatrices(w, dw, -a); //x = x - a * dw
            i++;
        }

        return w;
    }

    //Calculates the derivative of the loss function analytically.
    private static Matrix findDw(Matrix x, Matrix y, Matrix w) {

        Matrix error = LinearAlgebra.addMatrices(LinearAlgebra.multiplyMatrices(x, w), y, -1.0);
        Matrix dw = LinearAlgebra.multiplyMatrices(LinearAlgebra.transpose(x), error);
        dw = LinearAlgebra.scaleMatrix(dw, 1.0 / x.getRows());

        return dw;
    }

    //CURRENTLY ONLY WORKS FOR FUNCTIONS DEPENDENT ON ONE INDEPENDENT VARIABLE.
    public static Matrix buildPolyFunction(Matrix xVals, Matrix w) {
        Matrix y = new Matrix(xVals.getRows(), 1);
        int n = w.getRows(); //Polynomial order
        double x;
        double fnc;

        for (int i = 1; i <= xVals.getRows(); i++) {
            x = xVals.getValue(i, 1);
            fnc = 0;

            for (int j = 0; j < n; j++) {
                fnc += w.getValue(j + 1, 1) * Math.pow(x, j);
            }

            y.setValue(i, 1, fnc);
        }

        return y;
    }

    //ToDo: GradDesc but instead of x^n for additional parameters I use a method to apply ANY function to a matrix.
    //To do this I will need a set of basis functions in the form of an array of functional interfaces

}