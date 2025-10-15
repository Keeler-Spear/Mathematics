package Mathematics;

import java.util.function.Function;

/**
 * A static class that performs regression operations via numerical methods.
 *
 * @author Keeler Spear
 * @version %I%, %G%
 * @since 1.0
 */
public class Regression {

    final static double TOL = 0.01;
    final static int MAX_ITERATIONS = 500000;
    final static double LR = 0.001;
    final static double RAND_BOUND = 10;
    final static int MAX_CLASSES = 100;
    public static Function<Double, Double> sigmoid = (z) -> 1 / (1 + Math.exp(-z));

    /**
     * Calculates the weights for a polynomial regression model based on the data set provided using gradient descent.
     * <p>
     * The weights will form the following polynomial function:
     * w + wx + wx^2 + wx^3 + ... + wx^n + wy + wy^2 + wy^3 + ...
     * </p>
     *
     * @param x A matrix of data parameters.
     * @param y A vector of data labels.
     * @param n The order of the polynomial.
     * @return A matrix with the weights for a polynomial regression model based on the data set provided.
     * @throws IllegalArgumentException If the data does not have one sample for each label.
     * @throws IllegalArgumentException If n is less than one.
     * @see #linearReg(Matrix, Matrix, Matrix, double, Function[])
     */
    public static Matrix polyFit(Matrix x, Matrix y, int n) {
        if (x.getRows() != y.getRows()) {
            throw new IllegalArgumentException("The data does not have one sample for each label!");
        }

        if (n < 1) {
            throw new IllegalArgumentException("The polynomial's order must be equal to or greater than one!");
        }

        Matrix w = LinearAlgebra.randMatrix(n * x.getCols() + 1, 1, -RAND_BOUND, RAND_BOUND);
        double a;

        if (n <= 6) {
            a = LR;
        }
        else if (n <= 12) {
            a = LR / 10;
        }
        else if (n <= 18) {
            a = LR / 100;
        }
        else {
            a = LR / 1000;
        }

        return polyFit(x, y, w, n, a);
    }

    /**
     * Calculates the weights for a polynomial regression model based on the data set provided using gradient descent.
     * <p>
     * The weights will form the following polynomial function:
     * w + wx + wx^2 + wx^3 + ... + wx^n + wy + wy^2 + wy^3 + ...
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
     * @see #linearReg(Matrix, Matrix, Matrix, double, Function[])
     */
    public static Matrix polyFit(Matrix x, Matrix y, Matrix w, int n, double a) {
        if (x.getRows() != y.getRows()) {
            throw new IllegalArgumentException("The data does not have one sample for each label!");
        }

        if (n < 1) {
            throw new IllegalArgumentException("The polynomial's order must be equal to or greater than one!");
        }

        Function[] fncs = BasisFunctions.poly(n);

        return linearReg(x, y, w, a, fncs);
    }

    /**
     * Calculates the weights for a trigonometric regression model based on the data set provided using gradient descent.
     * <p>
     * The weights will form the following trigonometric function:
     * w + wsin(x) + wcos(x) + wsin(2x) + wcos(2x) + ... + wsin(nx) + wcos(nx) + wsin(y) + wcos(y) + wsin(2y)
     * + wcos(2y) + ... + wsin(ny) + wcos(ny) + ...
     * </p>
     *
     * @param x A matrix of data parameters.
     * @param y A vector of data labels.
     * @param w A vector of initial weight guesses.
     * @param n The order of the trigonometric function.
     * @param a The learning rate of the model.
     * @return A matrix with the weights for a trigonometric regression model based on the data set provided.
     * @throws IllegalArgumentException If the data does not have one sample for each label.
     * @throws IllegalArgumentException If n is less than one.
     * @see #linearReg(Matrix, Matrix, Matrix, double, Function[])
     */
    public static Matrix trigFit(Matrix x, Matrix y, Matrix w, int n, double a) {
        if (x.getRows() != y.getRows()) {
            throw new IllegalArgumentException("The data does not have one sample for each label!");
        }

        if (n < 1) {
            throw new IllegalArgumentException("The trigonometric function's order must be equal to or greater than one!");
        }

        Function[] fncs = BasisFunctions.trig(n);

        return linearReg(x, y, w, a, fncs);
    }



    /**
     * Calculates the weights for a linear regression model based on the data set provided and the basis functions
     * provided using gradient descent.
     * <p>
     * The weights will form the following function:
     * w + w1f(x) + wf2(x) + .. + wfn(x) + w1f(y) + wf2(y) + .. + wfn(y) + ...
     * </p>
     * <p>
     * In its current state, this method does have a few limitations. First, for very large xOrg data, the algorithm
     * will diverge if the learning rate is not small enough. Second, it can only handle data with a single label.
     * In the future, I plan to expand the algorithm to cover this. Third, it cannot handle functions with products,
     * quotients, compositions, etc. of multiple variables, such as z = x * y or w = y * w^2 * cos(x).
     * </p>
     *
     * @param xOrg A matrix of data parameters.
     * @param y A vector of data labels.
     * @param w A vector of initial weight guesses.
     * @param a The learning rate of the model.
     * @param fncs The set of basis functions the data will be modeled with.
     * @return A matrix with the weights for a linear regression model based on the data set provided.
     * @throws IllegalArgumentException If the data does not have one sample for each label.
     * @throws IllegalArgumentException If there is not one weight for each parameter.
     */
    public static Matrix linearReg(Matrix xOrg, Matrix y, Matrix w, double a, Function[] fncs) {
        if (xOrg.getRows() != y.getRows()) {
            throw new IllegalArgumentException("The data does not have one sample for each label!");
        }

        if (w.getRows() != (fncs.length - 1) * xOrg.getCols() + 1) {
            throw new IllegalArgumentException("There must be one weight for each parameter!");
        }

        Matrix x = new Matrix (xOrg.getRows(), 1);

        for (int i = 1; i <= x.getRows(); i++) {
            x.setValue(i, 1, (Double) fncs[0].apply(xOrg.getValue(i, 1)));
        }

        int numVars = xOrg.getCols();
        Matrix xn;

        //Adding the basis functions to the matrix
        for (int i = 1; i <= numVars; i++) {//For each variable
            for (int j = 1; j < fncs.length; j++) {//For each function order
                xn = LinearAlgebra.vectorFromColumn(xOrg, i);
                xn = LinearAlgebra.applyFunction(xn, fncs[j]);
                x.addColRight(xn.getMatrix());
            }
        }

        Matrix dw = LinearAlgebra.constantMatrix(w.getRows(), 1, 5.0);

        int i = 0;
        while (LinearAlgebra.l1Norm(dw) > TOL && i < MAX_ITERATIONS) {
            dw = findDw(x, y, LinearAlgebra.multiplyMatrices(x, w));
            w = LinearAlgebra.addMatrices(w, dw, -a); //x = x - a * dw
            i++;
        }

        return w;
    }

    //Calculates the derivative of the loss function analytically.
    private static Matrix findDw(Matrix x, Matrix y, Matrix h) {

        Matrix error = LinearAlgebra.addMatrices(h, y, -1.0); //h(w) - y
        Matrix dw = LinearAlgebra.multiplyMatrices(LinearAlgebra.transpose(x), error); //x * [h(w) - y]
        dw = LinearAlgebra.scaleMatrix(dw, 1.0 / x.getRows());

        return dw;
    }

    /**
     * Calculates the weights for a logistic regression model based on the data set provided and the basis functions provided
     * using gradient descent.
     * <p>
     * The weights will form the following function:
     * 1 / (1 + exp[w + w1f(x) + wf2(x) + .. + wfn(x) + w1f(y) + wf2(y) + .. + wfn(y) + ...])
     * </p>
     * <p>
     * In its current state, this method does have a few limitations. First, for very large xOrg data, the algorithm
     * will diverge if the learning rate is not small enough. Second, it can only handle data with a single label.
     * In the future, I plan to expand the algorithm to cover this. Third, it cannot handle functions with products,
     * quotients, compositions, etc. of multiple variables, such as z = x * y or w = y * w^2 * cos(x).
     * </p>
     *
     * @param xOrg A matrix of data parameters.
     * @param y A vector of data labels.
     * @param w A vector of initial weight guesses.
     * @param a The learning rate of the model.
     * @param fncs The set of basis functions the data will be modeled with.
     * @return A matrix with the weights for a logistic regression model based on the data set provided.
     * @throws IllegalArgumentException If the data does not have one sample for each label.
     * @throws IllegalArgumentException If there is not one weight for each parameter.
     */
    public static Matrix logisticReg(Matrix xOrg, Matrix y, Matrix w, double a, Function[] fncs, boolean rounding) {
        if (xOrg.getRows() != y.getRows()) {
            throw new IllegalArgumentException("The data does not have one sample for each label!");
        }

        if (w.getRows() != (fncs.length - 1) * xOrg.getCols() + 1) {
            throw new IllegalArgumentException("There must be one weight for each parameter!");
        }

        Matrix x = new Matrix (xOrg.getRows(), 1);

        for (int i = 1; i <= x.getRows(); i++) {
            x.setValue(i, 1, (Double) fncs[0].apply(xOrg.getValue(i, 1)));
        }

        int numVars = xOrg.getCols();
        Matrix xn;

        //Adding the basis functions to the matrix
        for (int i = 1; i <= numVars; i++) {//For each variable
            for (int j = 1; j < fncs.length; j++) {//For each function order
                xn = LinearAlgebra.vectorFromColumn(xOrg, i);
                xn = LinearAlgebra.applyFunction(xn, fncs[j]);
                x.addColRight(xn.getMatrix());
            }
        }

        Matrix dw = LinearAlgebra.constantMatrix(w.getRows(), 1, 5.0);

        //Gradient descent
        int i = 0;
        if (!rounding) {
            while (LinearAlgebra.l1Norm(dw) > TOL && i < MAX_ITERATIONS) {
                dw = findDw(x, y, LinearAlgebra.applyFunction(LinearAlgebra.multiplyMatrices(x, w), sigmoid));
                w = LinearAlgebra.addMatrices(w, dw, -a); //x = x - a * dw
                i++;
            }
        }

        else {
            while (LinearAlgebra.l1Norm(dw) > TOL && i < MAX_ITERATIONS) {
                dw = findDw(x, y, LinearAlgebra.roundMatrix(LinearAlgebra.applyFunction(LinearAlgebra.multiplyMatrices(x, w), sigmoid)));
                w = LinearAlgebra.addMatrices(w, dw, -a); //x = x - a * dw
                i++;
            }
        }

        return w;
    }

    /**
     * Calculates the weights for a classification model based on the data set, basis functions, and activation function provided
     * using gradient descent.
     * <p>
     * The weights will form the following function:
     * f[w + w1f(x) + wf2(x) + .. + wfn(x) + w1f(y) + wf2(y) + .. + wfn(y) + ...]
     * </p>
     * <p>
     * In its current state, this method does have a few limitations. First, for very large xOrg data, the algorithm
     * will diverge if the learning rate is not small enough. Second, it can only handle data with a single label.
     * In the future, I plan to expand the algorithm to cover this. Third, it cannot handle functions with products,
     * quotients, compositions, etc. of multiple variables, such as z = x * y or w = y * w^2 * cos(x).
     * </p>
     *
     * @param xOrg A matrix of data parameters.
     * @param y A vector of data labels.
     * @param w A vector of initial weight guesses.
     * @param a The learning rate of the model.
     * @param activationFunction The activation function to be used.
     * @param fncs The set of basis functions the data will be modeled with.
     * @return A matrix with the weights for a logistic regression model based on the data set provided.
     * @throws IllegalArgumentException If the data does not have one sample for each label.
     * @throws IllegalArgumentException If there is not one weight for each parameter.
     */
    public static Matrix classificationReg(Matrix xOrg, Matrix y, Matrix w, double a, Function<Double, Double> activationFunction, Function[] fncs, boolean rounding) {
        if (xOrg.getRows() != y.getRows()) {
            throw new IllegalArgumentException("The data does not have one sample for each label!");
        }

        if (w.getRows() != (fncs.length - 1) * xOrg.getCols() + 1) {
            throw new IllegalArgumentException("There must be one weight for each parameter!");
        }

        Matrix x = new Matrix (xOrg.getRows(), 1);

        for (int i = 1; i <= x.getRows(); i++) {
            x.setValue(i, 1, (Double) fncs[0].apply(xOrg.getValue(i, 1)));
        }

        int numVars = xOrg.getCols();
        Matrix xn;

        //Adding the basis functions to the matrix
        for (int i = 1; i <= numVars; i++) {//For each variable
            for (int j = 1; j < fncs.length; j++) {//For each function order
                xn = LinearAlgebra.vectorFromColumn(xOrg, i);
                xn = LinearAlgebra.applyFunction(xn, fncs[j]);
                x.addColRight(xn.getMatrix());
            }
        }

        Matrix dw = LinearAlgebra.constantMatrix(w.getRows(), 1, 5.0);

        //Gradient descent
        int i = 0;
        if (!rounding) {
            while (LinearAlgebra.l1Norm(dw) > TOL && i < MAX_ITERATIONS) {
                dw = findDw(x, y, LinearAlgebra.applyFunction(LinearAlgebra.multiplyMatrices(x, w), activationFunction));
                w = LinearAlgebra.addMatrices(w, dw, -a); //x = x - a * dw
                i++;
            }
        }

        else {
            while (LinearAlgebra.l1Norm(dw) > TOL && i < MAX_ITERATIONS) {
                dw = findDw(x, y, LinearAlgebra.roundMatrix(LinearAlgebra.applyFunction(LinearAlgebra.multiplyMatrices(x, w), activationFunction)));
                w = LinearAlgebra.addMatrices(w, dw, -a); //x = x - a * dw
                i++;
            }
        }

        return w;
    }


    /**
     * Calculates the values for a function at the provided points using the functions and weights provided.
     *
     * @param x A matrix of points.
     * @param w A vector weights.
     * @param fncsOrg The set of basis functions that comprise the function.
     * @return A vector of values for a function.
     * @throws IllegalArgumentException If the weights matrix has more than one weight.
     */
    public static Matrix buildFunction(Matrix x, Matrix w, Function[] fncsOrg) {
        if (w.getCols() != 1) {
            throw new IllegalArgumentException("The weight matrix must have exactly one column!");
        }

        Matrix temp = LinearAlgebra.vectorFromColumn(x, 1);
        Matrix y = LinearAlgebra.scaleMatrix(LinearAlgebra.applyFunction(temp, fncsOrg[0]), w.getValue(1, 1)); //y = [1; 1; ...; 1]

        //Removing the bias from basis functions
        Function[] fncs = new Function[fncsOrg.length - 1];
        for (int i = 0; i < fncs.length; i++) {
            fncs[i] = fncsOrg[i + 1];
        }

        for (int i = 0; i < fncs.length; i++) {
            for (int j = 1; j <= x.getCols(); j++) {
                temp = LinearAlgebra.vectorFromColumn(x, j);
                y = LinearAlgebra.addMatrices(y, LinearAlgebra.applyFunction(temp, fncs[i]), w.getValue((j - 1) * (fncs.length) + i + 2, 1));
            }
        }

        return y;
    }

    /**
     * Calculates the values for a function at the provided points using the functions and weights provided.
     *
     * @param x A matrix of points.
     * @param w A vector weights.
     * @param fncsOrg The set of basis functions that comprise the function within the sigmoid function.
     * @return A vector of values for a function.
     * @throws IllegalArgumentException If the weights matrix has more than one weight.
     */
    public static Matrix buildLogisticFunction(Matrix x, Matrix w, Function[] fncsOrg) {
        if (w.getCols() != 1) {
            throw new IllegalArgumentException("The weight matrix must have exactly one column!");
        }

        Matrix y = buildFunction(x, w, fncsOrg);

        //Only works for binary classification
        return LinearAlgebra.applyFunction(y, sigmoid);
    }
}
