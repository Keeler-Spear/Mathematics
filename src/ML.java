import java.util.Random;

/**
 * A static class that performs the operations essential to machine learning.
 * <p>
 *     This class and all functions within were added after I submitted my application to your institution.
 * </p>
 * @author Keeler Spear
 * @version %I%, %G%
 * @since 1.0
 */
public class ML {
    final static double TOL = 0.01;
    final static int MAX_ITERATIONS = 15000;
    final static double LR = 0.001;
    final static double RAND_MAX = 50;


    /**
     * Calculates the weights for a linear regression model based on the data set provided using gradient descent.
     *
     * @param x A matrix of data parameters.
     * @param y A vector of data labels.
     * @return A matrix with the weights for a linear regression model based on the data set provided.
     * @throws IllegalArgumentException If the data does not have one sample for each label.
     * @see #gradDes(Matrix, Matrix, Matrix, double, double)
     */
    public static Matrix gradDes(Matrix x, Matrix y) {
        if (x.getRows() != y.getRows()) {
            throw new IllegalArgumentException("The data does not have one sample for each label!");
        }

        Matrix w = LinearAlgebra.randMatrix(x.getCols(), 1, RAND_MAX);
        Random rand = new Random();
        double b = rand.nextDouble(RAND_MAX);

        return gradDes(x, y, w, b, LR);
    }

    /**
     * Calculates the weights for a linear regression model based on the data set provided using gradient descent.
     *
     * @param x A matrix of data parameters.
     * @param y A vector of data labels.
     * @param w A vector of initial weight guesses.
     * @param b A scalar guess for the bias.
     * @param a The learning rate of the model.
     * @return A matrix with the weights for a linear regression model based on the data set provided.
     * @throws IllegalArgumentException If the data does not have one sample for each label.
     * @throws IllegalArgumentException If there is not one weight for each parameter.
     */
    public static Matrix gradDes(Matrix x, Matrix y, Matrix w, double b, double a) {
        if (x.getRows() != y.getRows()) {
            throw new IllegalArgumentException("The data does not have one sample for each label!");
        }

        if (w.getRows() != x.getCols()) {
            throw new IllegalArgumentException("There must be one weight for each parameter!");
        }

        Matrix dw = findDw(x, y, w, b);
        int i = 0;

        while (LinearAlgebra.l2Norm(dw) > TOL && i < MAX_ITERATIONS) {
            b = b - a * findDb(x, y, w, b);
            dw = findDw(x, y, w, b);
            w = LinearAlgebra.addMatrices(w, dw, -a); //x = x - a * dw
            i++;
        }

        Matrix weights = new Matrix(w.getRows() + 1, 1);
        weights.setValue(1, 1, b);

        for (int j = 1; j <= w.getRows(); j++) {
            weights.setValue(j + 1, 1, w.getValue(j, 1));
        }
        return weights;
    }

    //Calculates the derivative of the loss function analytically.
    private static Matrix findDw(Matrix x, Matrix y, Matrix w, double b) {
        double val;
        Matrix dw = new Matrix(w.getRows(), 1);
        Matrix xj;


        for (int i = 1; i <= w.getRows(); i++) {
            val = 0.0;
            for (int j = 1; j <= x.getRows(); j++) {
                xj = LinearAlgebra.vectorFromRow(x, j);
                val += (LinearAlgebra.dotProduct(w, xj) - y.getValue(j, 1) + b) * xj.getValue(i, 1); //Replace one with j
            }
            dw.setValue(i, 1, val / x.getRows());
        }
        return dw;
    }

    //Calculates the derivative of the bias analytically.
    private static double findDb(Matrix x, Matrix y, Matrix w, double b) {
        Matrix xi;
        double val = 0.0;

        for (int i = 1; i <= x.getRows(); i++) {
            xi = LinearAlgebra.vectorFromRow(x, i);
            val += (LinearAlgebra.dotProduct(w, xi) - y.getValue(i, 1) + b);
        }

        return val / x.getRows();
    }

}
