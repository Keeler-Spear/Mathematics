/**
 * A static class that performs the operations essential to machine learning.
 * <p>
 *     This class and all functions within were added after I submitted my application to your institution.
 * </p>
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
     *
     * @param x A matrix of data parameters.
     * @param y A vector of data labels.
     * @param n The order of the polynomial.
     * @return A matrix with the weights for a polynomial regression model based on the data set provided.
     * @throws IllegalArgumentException If the data does not have one sample for each label.
     * @throws IllegalArgumentException If n is less than one.
     * @see #polyGradDes(Matrix, Matrix, Matrix, int, double)
     */
    public static Matrix polyGradDes(Matrix x, Matrix y, int n) {
        if (x.getRows() != y.getRows()) {
            throw new IllegalArgumentException("The data does not have one sample for each label!");
        }

        if (n < 1) {
            throw new IllegalArgumentException("The polynomial's order must be equal to or greater than one!");
        }

        Matrix w = LinearAlgebra.randMatrix(n * x.getCols() + 1, 1, -RAND_BOUND, RAND_BOUND);

        double a = Math.pow(10, 1 - n);

        return polyGradDes(x, y, w, n, a);
    }

    /**
     * Calculates the weights for a polynomial regression model based on the data set provided using gradient descent.
     * This code only works for one input variable.
     *
     * @param x A matrix of data parameters.
     * @param y A vector of data labels.
     * @param w A vector of initial weight guesses.
     * @param n The order of the polynomial.
     * @param a The learning rate of the model.
     * @return A matrix with the weights for a polynomial regression model based on the data set provided.
     * @throws IllegalArgumentException If the data does not have one sample for each label.
     * @throws IllegalArgumentException If there is not one weight for each parameter.
     * @throws IllegalArgumentException If n is less than 1.
     */
    public static Matrix polyGradDes(Matrix x, Matrix y, Matrix w, int n, double a) {
        if (x.getRows() != y.getRows()) {
            throw new IllegalArgumentException("The data does not have one sample for each label!");
        }

        if (w.getRows() != n * x.getCols() + 1) {
            throw new IllegalArgumentException("There must be one weight for each parameter!");
        }

        if (n < 1) {
            throw new IllegalArgumentException("The polynomial's order must be equal to or greater than 1!");
        }

        //Setting up x^n's
        int numVars = x.getCols();
        Matrix xn;
        Matrix xTemp;

        for (int i = 1; i <= numVars; i++) {//For each variable
            for (int j = 2; j <= n; j++) {//For each polynomial order
                xn = LinearAlgebra.vectorFromColumn(x, i);
                xTemp = xn;
                for (int k = 1; k < j; k++) {
                    xn = LinearAlgebra.multiplyValues(xn, xTemp);
                }
                x.addColRight(xn.getMatrix());
            }
        }

        //For bias term
        Matrix ones = LinearAlgebra.constantMatrix(x.getRows(), 1, 1.0);
        x.addColLeft(ones.getMatrix());

        Matrix dw = LinearAlgebra.constantMatrix(w.getRows(), 1, 5.0);

        int i = 0;
        while (LinearAlgebra.l2Norm(dw) > TOL && i < MAX_ITERATIONS) {
            dw = findDw(x, y, w);
            w = LinearAlgebra.addMatrices(w, dw, -a); //x = x - a * dw
            i++;
        }

        //For some reason the bias is half the size it should be so the line below fixes that.
        w.setValue(1, 1, w.getValue(1, 1) * 2);

        return w;
    }

    //Calculates the derivative of the loss function analytically.
    private static Matrix findDw(Matrix x, Matrix y, Matrix w) {

        Matrix error = LinearAlgebra.addMatrices(LinearAlgebra.multiplyMatrices(x, w), y, -1.0);
        Matrix dw = LinearAlgebra.multiplyMatrices(LinearAlgebra.transpose(x), error);
        dw = LinearAlgebra.scaleMatrix(dw, 1.0 / x.getRows());

        return dw;
    }

    //    /**
//     * Calculates the weights for a linear regression model based on the data set provided using gradient descent.
//     * While doing so, the path of x and y is recorded.
//     *
//     * @param x A matrix of data parameters.
//     * @param y A vector of data labels.
//     * @return A matrix with the weights for a linear regression model based on the data set provided and the weights'
//     * paths.
//     * @throws IllegalArgumentException If the data does not have one sample for each label.
//     * @see #gradDes(Matrix, Matrix, Matrix, double, double)
//     */
//    public static Matrix[] gradDesTrack(Matrix x, Matrix y) {
//        if (x.getRows() != y.getRows()) {
//            throw new IllegalArgumentException("The data does not have one sample for each label!");
//        }
//
//        //Adjust x here too
//
//        Matrix w = LinearAlgebra.randMatrix(x.getCols(), 1, -RAND_BOUND, RAND_BOUND);
//        Random rand = new Random();
//        double b = rand.nextDouble(RAND_BOUND);
//
//        return gradDesTrack(x, y, w, b, LR);
//    }

//    /**
//     * Calculates the weights for a linear regression model based on the data set provided using gradient descent.
//     * While doing so, the path of x and y is recorded.
//     *
//     * @param x A matrix of data parameters.
//     * @param y A vector of data labels.
//     * @param w A vector of initial weight guesses.
//     * @param b A scalar guess for the bias.
//     * @param a The learning rate of the model.
//     * @return A matrix with the weights for a linear regression model based on the data set provided and the weights'
//     * paths.
//     * @throws IllegalArgumentException If the data does not have one sample for each label.
//     * @throws IllegalArgumentException If there is not one weight for each parameter.
//     */
//    public static Matrix[] gradDesTrack(Matrix x, Matrix y, Matrix w, double b, double a) {
//        if (x.getRows() != y.getRows()) {
//            throw new IllegalArgumentException("The data does not have one sample for each label!");
//        }
//
//        if (w.getRows() != x.getCols()) {
//            throw new IllegalArgumentException("There must be one weight for each parameter!");
//        }
//
//        Matrix dw = findDw(x, y, w, b);
//        Matrix wPath = new Matrix(1, MAX_ITERATIONS);
//        Matrix bPath = new Matrix(1, MAX_ITERATIONS);
//        int i = 0;
//
//        while (LinearAlgebra.l2Norm(dw) > TOL && i < MAX_ITERATIONS) {
//            b = b - a * findDb(x, y, w, b);
//            bPath.setValue(1, i + 1, b);
//            dw = findDw(x, y, w, b);
//            w = LinearAlgebra.addMatrices(w, dw, -a); //x = x - a * dw
//            wPath.setValue(1, i + 1, w.getValue(1, 1));
//            i++;
//        }
//
//        Matrix weights = new Matrix(w.getRows() + 1, 1);
//        weights.setValue(1, 1, b);
//
//        for (int j = 1; j <= w.getRows(); j++) {
//            weights.setValue(j + 1, 1, w.getValue(j, 1));
//        }
//
//        Matrix bVals = new Matrix(1, i);
//        Matrix wVals = new Matrix(1, i);
//        for (int j = 1; j <= i; j++) {
//            bVals.setValue(1, j, bPath.getValue(1, j));
//            wVals.setValue(1, j, wPath.getValue(1, j));
//        }
//
//        Matrix[] temp = {weights, bVals, wVals};
//        return temp;
//    }

}