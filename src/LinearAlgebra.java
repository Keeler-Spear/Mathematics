import java.util.Random;
import java.util.function.Function;

/**
 * A static class that performs mathematical operations on one or more matrices via analytic and numerical methods.
 * <p>
 *     This class supports matrix operations such as addition, multiplication, Gaussian Elimination, and more.
 * </p>
 *
 * @author Keeler Spear
 * @version %I%, %G%
 * @since 1.0
 */
public class LinearAlgebra {
    final static double TOL = 0.000001;
    final static int MAX_ITERATIONS = 1000;

    //Computes if the provided value is "zero."
    private static boolean isZero(double val) {
        if (Math.abs(val) < TOL) {
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * Creates an nxn identity matrix.
     *
     * @param n The size that the identity matrix will be.
     * @return An nxn identity matrix.
     * @throws IllegalArgumentException If n less than or equal to zero.
     */
    public static Matrix identityMatrix(int n) {
        if (n <= 0) {
            throw new IllegalArgumentException("The matrix's size must be at least 1x1!");
        }

        double[][] identity = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                identity[i][j] = 0.0;
                if (i == j) {
                    identity[i][j] = 1;
                }
            }
        }

        return new Matrix(identity);
    }

    /**
     * Creates an nxn identity matrix and then multiplies it by a constant.
     *
     * @param n The size that the identity matrix will be.
     * @param c constant that the identity matrix will be multiplied by.
     * @return An nxn identity matrix multiplied by a given constant.
     * @throws IllegalArgumentException If n less than or equal to zero.
     * @see #scaleMatrix(Matrix, double)
     */
    public static Matrix constantIdentityMatrix(int n, double c) {
        if (n <= 0) {
            throw new IllegalArgumentException("The matrix's size must be at least 1x1!");
        }

        Matrix identity = identityMatrix(n);

        return scaleMatrix(identity, c);
    }

    /**
     * Creates an mxn matrix with cells holding the value provided. This method was added after I submitted my
     * application to your institution.
     *
     * @param rows The number of rows that the matrix will have.
     * @param cols The number of columns that the matrix will have.
     * @param c constant that matrix will contain.
     * @return An mxn matrix with cells holding the value provided.
     * @throws IllegalArgumentException If rows or cols is less than or equal to zero.
     */
    public static Matrix constantMatrix(int rows, int cols, double c) {
        if (rows <= 0 || cols <= 0) {
            throw new IllegalArgumentException("The matrix's size must be at least 1x1!");
        }

        Matrix constant = new Matrix(rows, cols);

        for (int i = 1; i <= rows; i++) {
            for (int j = 1; j <= cols; j++) {
                constant.setValue(i, j, c);
            }
        }

        return constant;
    }

    /**
     * Creates a zero matrix object of the size provided.
     *
     * @param rows The number of rows the matrix will have.
     * @param cols The number of columns the matrix will have.
     * @throws IllegalArgumentException If the matrix's size is less than 1.
     * @return A zero matrix of the size provided.
     */
    public static Matrix zeroMatrix (int rows, int cols) {
        return constantMatrix(rows, cols, 0.0);
    }

    /**
     * Creates a matrix of the size provided filled with random values. This method was added after I submitted my
     * application to your institution.
     *
     * @param rows The number of rows the matrix will have.
     * @param cols The number of columns the matrix will have.
     * @param min The lower bound of random number generation.
     * @param max The upper bound of random number generation.
     * @throws IllegalArgumentException If the matrix's size is less than 1.
     * @throws IllegalArgumentException If the minimum value is greater than or equal to the maximum value.
     * @return a matrix of the size provided filled with random values
     */
    public static Matrix randMatrix (int rows, int cols, double min, double max) {
        if (rows <= 0 || cols <= 0) {
            throw new IllegalArgumentException("The matrix's size must be at least 1x1!");
        }
        if (min >= max) {
            throw new IllegalArgumentException("The minimum value must be less than the maximum value! Min = " + min + ", Max = " + max + "!");
        }

        Random rand = new Random();
        Matrix matrix = new Matrix(rows, cols);

        for (int i = 1; i <= rows; i++) {
            for (int j = 1; j <= cols; j++) {
                matrix.setValue(i, j, min + (max - min) * rand.nextDouble());
            }
        }

        return matrix;
    }

    /**
     * Creates a matrix of the size provided filled with random values. This method was added after I submitted my
     * application to your institution.
     *
     * @param rows The number of rows the matrix will have.
     * @param cols The number of columns the matrix will have.
     * @param min The lower bound of random number generation.
     * @param max The upper bound of random number generation.
     * @param seed The seed use for random number generation.
     * @throws IllegalArgumentException If the matrix's size is less than 1.
     * @throws IllegalArgumentException If the minimum value is greater than or equal to the maximum value.
     * @return a matrix of the size provided filled with random values
     */
    public static Matrix randMatrix (int rows, int cols, double min, double max, long seed) {
        if (rows <= 0 || cols <= 0) {
            throw new IllegalArgumentException("The matrix's size must be at least 1x1!");
        }
        if (min >= max) {
            throw new IllegalArgumentException("The minimum value must be less than the maximum value!");
        }

        Random rand = new Random();
        rand.setSeed(seed);
        Matrix matrix = new Matrix(rows, cols);

        for (int i = 1; i <= rows; i++) {
            for (int j = 1; j <= cols; j++) {
                matrix.setValue(i, j, min + (max - min) * rand.nextDouble());
            }
        }

        return matrix;
    }

    /**
     * Creates a vector filled with the values between min and max, inclusive. This method was added after I submitted
     * my application to your institution.
     *
     * @param min The minimum or initial value of the matrix.
     * @param max The maximum or final value of the matrix.
     * @param inc The step size between values.
     * @throws IllegalArgumentException If the minimum value is greater than or equal to the maximum value.
     * @return A vector filled with the values between min and max, inclusive.
     */
    public static Matrix linSpace(double min, double max, double inc) {
        if (min >= max) {
            throw new IllegalArgumentException("The minimum value must be less than the maximum value!");
        }

        double[] vals = new double[1 + (int) ((max - min) / inc)];

        for (int i = 0; i < vals.length; i++) {
            vals[i] = min;
            min += inc;
        }

        return new Matrix(vals);
    }

    /**
     * Creates a matrix whose values are those of the input matrix's after a function has been applied to them. This
     * method was added after I submitted my application to your institution.
     *
     * @param matrixOrg The matrix that the function will be applied to.
     * @param fnc The function that will be applied to the matrix.
     * @return A matrix whose values are those of the input matrix's after a function has been applied to them.
     */
    public static Matrix applyFunction(Matrix matrixOrg, Function<Double, Double> fnc) {
        Matrix A = new Matrix(matrixOrg.getRows(), matrixOrg.getCols());

        for (int i = 1; i <= A.getRows(); i++) {
            for(int j = 1; j <= A.getCols(); j++) {
                A.setValue(i, j, fnc.apply(matrixOrg.getValue(i, j)));
            }
        }

        return A;
    }

    /**
     * Creates a matrix that is the transpose of the one provided.
     *
     * @param A The matrix that will be copied and transposed.
     * @return A copy of the provided matrix, transposed.
     */
    public static Matrix transpose(Matrix A) {
        Matrix matrix = new Matrix (A.getCols(), A.getRows());

        for (int i = 1; i <= A.getCols(); i++) {
            for (int j = 1; j <= A.getRows(); j++) {
                matrix.setValue(i, j, A.getValue(j, i));
            }
        }

        return matrix;
    }

    /**
     * Creates a matrix that is augmented on the right with another matrix.
     *
     * @param A The matrix that will be augmented.
     * @param B The matrix that will be augmented onto A's right side.
     * @return A copy of the matrix A, augmented with Matrix B.
     * @throws IllegalArgumentException If the matrices' row counts are not the same.
     */
    public static Matrix augmentMatrix(Matrix A, Matrix B) {
        if (A.getRows() != B.getRows()) {
            throw new IllegalArgumentException("The matrices must have the same number of rows!");
        }

        Matrix matrix = new Matrix (A.getRows(), A.getCols() + B.getCols());

        for (int i = 1; i <= A.getCols(); i++) {
            for (int j = 1; j <= A.getRows(); j++) {
                matrix.setValue(j, i, A.getValue(j, i));
            }
        }

        for (int i = 1; i <= B.getCols(); i++) {
            for (int j = 1; j <= A.getRows(); j++) {
                matrix.setValue(j, i + A.getCols(), B.getValue(j, i));
            }
        }

        return matrix;
    }

    /**
     * Creates a matrix whose values are [A; B]. This method was added after I submitted my application to your
     * institution.
     *
     * @param A The matrix that will be stacked on top of the other.
     * @param B The matrix that will be below the other.
     * @return A matrix whose values are [A; B].
     * @throws IllegalArgumentException If the matrices' column counts are not the same.
     */
    public static Matrix stackMatrices(Matrix A, Matrix B) {
        if (A.getCols() != B.getCols()) {
            throw new IllegalArgumentException("The matrices must have the same number of columns!");
        }

        Matrix C = new Matrix(A.getRows() + B.getRows(), A.getCols());

        for (int i = 1; i <= A.getRows(); i++) {
            for (int j = 1; j <= A.getCols(); j++) {
                C.setValue(i, j, A.getValue(i, j));
            }
        }

        for (int i = 1; i <= B.getRows(); i++) {
            for (int j = 1; j <= B.getCols(); j++) {
                C.setValue(i + A.getRows(), j, B.getValue(i, j));
            }
        }

        return C;
    }

    /**
     * Creates a matrix whose entries will be rounded to the nearest integer. This method was added after I submitted my
     * application to your institution.
     *
     * @param matrixOrg The matrix that will be copied and rounded.
     * @return A copy of matrix matrixOrg, rounded to the nearest integer.
     */
    public static Matrix roundMatrix(Matrix matrixOrg) {
        Matrix A = new Matrix(matrixOrg.getRows(), matrixOrg.getCols());

        for (int i = 1; i <= A.getRows(); i++) {
            for (int j = 1; j <= A.getCols(); j++) {
                A.setValue(i, j, Math.round(matrixOrg.getValue(i, j)));
            }
        }

        return A;
    }

    /**
     * Creates a matrix whose entries will be rounded to the nth decimal place. This method was added after I submitted
     * my application to your institution.
     *
     * @param matrixOrg The matrix that will be copied and rounded.
     * @param n The number of decimal places the matrix's values will be rounded to.
     * @return A copy of matrix matrixOrg, rounded.
     */
    public static Matrix roundMatrix(Matrix matrixOrg, int n) {
        Matrix A = new Matrix(matrixOrg.getRows(), matrixOrg.getCols());
        double val;

        for (int i = 1; i <= A.getRows(); i++) {
            for (int j = 1; j <= A.getCols(); j++) {
                val = Math.round(matrixOrg.getValue(i, j) * Math.pow(10, n));
                A.setValue(i, j, val * Math.pow(10, -n));
            }
        }

        return A;
    }

    /**
     * Returns the sum of values in the row indicated. This method was added after I submitted
     * my application to your institution.
     *
     * @param A The matrix that will be used.
     * @param row The row index containing the row that will be summed.
     * @return The sum of values in the row indicated.
     * @throws IllegalArgumentException If the row index does not exist.
     */
    public static double rowSum(Matrix A, int row) {
        if (row < 1 || row > A.getRows()) {
            throw new IllegalArgumentException("This row index does not exist!");
        }

        double sum = 0;

        for (int i = 1; i <= A.getCols(); i++) {
            sum += A.getValue(row, i);
        }

        return sum;
    }

    /**
     * Returns the sum of values in the column indicated. This method was added after I submitted
     * my application to your institution.
     *
     * @param A The matrix that will be used.
     * @param col The column index containing the column that will be summed.
     * @return The sum of values in the column indicated.
     * @throws IllegalArgumentException If the column index does not exist.
     */
    public static double colSum(Matrix A, int col) {
        if (col < 1 || col > A.getCols()) {
            throw new IllegalArgumentException("This column index does not exist!");
        }

        double sum = 0;

        for (int i = 1; i <= A.getRows(); i++) {
            sum += A.getValue(i, col);
        }

        return sum;
    }

    /**
     * Returns the sum of values on the matrix's diagonal. This method was added after I submitted
     * my application to your institution.
     *
     * @param A The matrix that will be used.
     * @return The sum of values on the matrix's diagonal.
     * @throws IllegalArgumentException If the matrix is not square
     */
    public static double diagSum(Matrix A) {
        if (!A.isSquare()) {
            throw new IllegalArgumentException("This matrix is not square!");
        }

        double sum = 0;

        for (int i = 1; i <= A.getRows(); i++) {
            sum += A.getValue(i, i);
        }

        return sum;
    }

    /**
     * Returns the sum of values in the matrix provided. This method was added after I submitted
     * my application to your institution.
     *
     * @param A The matrix that will be used.
     * @return The sum of values in the matrix provided.
     */
    public static double matrixSum(Matrix A) {

        double sum = 0;

        for (int i = 1; i <= A.getRows(); i++) {
            sum += rowSum(A, i);
        }

        return sum;
    }

    /**
     * Creates a matrix that is multiplied by the scalar provided.
     *
     * @param matrixOrg The matrix that will be copied and scaled.
     * @param scalar The value that the copied matrix will be multiplied by.
     * @return A copy of matrix matrixOrg, multiplied by the provided scalar.
     */
    public static Matrix scaleMatrix(Matrix matrixOrg, double scalar) {
        Matrix A;

        try {
            A = (Matrix) matrixOrg.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        for (int i = 1; i <= A.getRows(); i++) {
            for (int j = 1; j <= A.getCols(); j++) {
                A.setValue(i, j, A.getValue(i, j) * scalar);
            }
        }

        return A;
    }

    /**
     * Creates a matrix whose entries are the absolute of the one provided.
     *
     * @param matrixOrg The matrix that will be copied and made absolute.
     * @return A copy of matrix matrixOrg, whose values are turned absolute.
     */
    public static Matrix abs(Matrix matrixOrg) {
        Matrix A = new Matrix(matrixOrg.getRows(), matrixOrg.getCols());

        for (int i = 1; i <= matrixOrg.getRows(); i++) {
            for (int j = 1; j <= matrixOrg.getCols(); j++) {
                A.setValue(i, j, Math.abs(matrixOrg.getValue(i, j)));
            }
        }

        return A;
    }

    /**
     * Adds one matrix to another matrix that is multiplied by a scalar. Matrix sum = matrix A + scalar * matrix B.
     *
     * @param A The matrix that will be added to the other.
     * @param B The matrix that will be multiplied by a scalar and added to the other.
     * @param scalar The value the second matrix (B) will be multiplied by in matrix addition.
     * @return The matrix that is the sum of matrix A and matrix B after B is multiplied by a scalar.
     * @throws IllegalArgumentException If the matrices A and B do not have the same dimensions.
     */
    public static Matrix addMatrices(Matrix A, Matrix B, double scalar) {
        if (A.getRows() != B.getRows() || A.getCols() != B.getCols()) {
            throw new IllegalArgumentException("The matrices must be of the same size!");
        }

        Matrix sum = new Matrix (A.getRows(), A.getCols());

        for (int i = 1; i <= A.getRows(); i++) {
            for (int j = 1; j <= A.getCols(); j++) {
                sum.setValue(i, j, A.getValue(i , j) + scalar * B.getValue(i, j));
            }
        }

        return sum;
    }

    /**
     * Adds one matrix to another matrix. Matrix sum = matrix A + matrix B.
     *
     * @param A The matrix that will be added to the other.
     * @param B The matrix that will be added to the other.
     * @return The matrix that is the sum of matrix A and matrix B.
     * @throws IllegalArgumentException If the matrices A and B do not have the same dimensions.
     */
    public static Matrix addMatrices(Matrix A, Matrix B) {
        if (A.getRows() != B.getRows() || A.getCols() != B.getCols()) {
            throw new IllegalArgumentException("The matrices must be of the same size!");
        }

        return addMatrices(A, B, 1.0);
    }

    /**
     * Subtracts one matrix to another matrix. Matrix difference = matrix A - matrix B.
     *
     * @param A The matrix that will be subtracted from.
     * @param B The matrix that will be subtracted from the other.
     * @return The matrix that is the difference of matrix A and matrix B.
     * @throws IllegalArgumentException If the matrices A and B do not have the same dimensions.
     */
    public static Matrix subtractMatrices(Matrix A, Matrix B) {
        if (A.getRows() != B.getRows() || A.getCols() != B.getCols()) {
            throw new IllegalArgumentException("The matrices must be of the same size!");
        }

        return addMatrices(A, B, -1.0);
    }

    /**
     * Multiplies the values of two matrices together. This method was added after I submitted my application to your
     * institution.
     *
     * @param A The first matrix.
     * @param B The second matrix.
     * @return A matrix whose values are the product of those in the matrices provided.
     * @throws IllegalArgumentException If the provided matrices are not the same size.
     */
    public static Matrix multiplyValues(Matrix A, Matrix B) {
        if (A.getRows() != B.getRows() || A.getCols() != B.getCols()) {
            throw new IllegalArgumentException("The matrices must be of the same size!");
        }

        Matrix mul = new Matrix (A.getRows(), A.getCols());

        for (int i = 1; i <= A.getRows(); i++) {
            for (int j = 1; j <= A.getCols(); j++) {
                mul.setValue(i, j, A.getValue(i, j) * B.getValue(i, j));
            }
        }

        return mul;
    }

    /**
     * Divides the values of two matrices together. Matrix quotient(i, j) = A(i, j) / B(i, j)
     *
     * @param A The first matrix.
     * @param B The second matrix.
     * @return A matrix whose values are the quotient of those in the matrices provided.
     * @throws IllegalArgumentException If the provided matrices are not the same size.
     */
    public static Matrix divideValues(Matrix A, Matrix B) {
        if (A.getRows() != B.getRows() || A.getCols() != B.getCols()) {
            throw new IllegalArgumentException("The matrices must be of the same size!");
        }

        Matrix quot = new Matrix (A.getRows(), A.getCols());

        for (int i = 1; i <= A.getRows(); i++) {
            for (int j = 1; j <= A.getCols(); j++) {
                quot.setValue(i, j, A.getValue(i, j) / B.getValue(i, j));
            }
        }

        return quot;
    }

    /**
     * Right multiplies a provided matrix by another. Matrix product = matrix A * matrix B.
     *
     * @param A The matrix that will be multiplied by the other on the right.
     * @param B The matrix that will be multiplied by the other on the left.
     * @return The matrix that is the product of matrix A being right multiplied by matrix B.
     * @throws IllegalArgumentException If matrix A does not have the same number of columns as matrix B has rows.
     */
    public static Matrix multiplyMatrices(Matrix A, Matrix B) {
        if (A.getCols() != B.getRows()) {
            throw new IllegalArgumentException("The second matrix must have a number of rows equal to the number of columns of the first! " + A.getRows() + " != " +  B.getRows() + "!");
        }

        Matrix C = new Matrix(A.getRows(), B.getCols());
        double sum;

        for (int i = 1; i <= A.getRows(); i++) {
            for (int j = 1; j <= B.getCols(); j++) {
                sum = 0.0;
                for (int k = 1; k <= A.getCols(); k++) {
                    sum += A.getValue(i, k) * B.getValue(k, j);
                }
                C.setValue(i, j, sum);
            }
        }

        return C;
    }


    /**
     * Multiplies a matrix by itself n times.
     *
     * @param A The matrix that will be multiplied n times.
     * @param n The number of time a matrix will be multiplied by itself.
     * @return The matrix that is the product of matrix A being multiplied by itself n times.
     * @throws IllegalArgumentException If the matrix provided is not square.
     * @throws IllegalArgumentException If the power is negative.
     */
    public static Matrix matrixPow(Matrix A, int n) {
        if (!A.isSquare()) {
            throw new IllegalArgumentException("The matrix provided must be square!");
        }

        if (n < 0) {
            throw new IllegalArgumentException("The power must be non-negative!");
        }

        if (n == 0) {
            return identityMatrix(A.getRows());
        }

        if (n == 1) {
            return A;
        }

        else {
            Matrix product = multiplyMatrices(A, A);

            for (int i = 2; i < n; i++) {
                product = multiplyMatrices(product, A);
            }

            return product;
        }
    }

    /**
     * Multiplies the values of the matrix by themselves n times.
     *
     * @param A The matrix whose values will be multiplied by themselves n times.
     * @param n The number of times the matrix's values will be multiplied by themselves.
     * @return A matrix whose values have been multiplied by themselves n times.
     * @throws IllegalArgumentException If the power is negative.
     */
    public static Matrix valuesPow(Matrix A, int n) {
        if (n < 0) {
            throw new IllegalArgumentException("The power must be non-negative!");
        }

        Matrix B = new Matrix(A.getRows(), A.getCols());

        for (int i = 1; i <= A.getRows(); i++) {
            for (int j = 1; j <= A.getCols(); j++) {
                B.setValue(i, j, Math.pow(A.getValue(i, j), n));
            }
        }

        return B;
    }

    /**
     * Constructs an mx1 matrix from the column specified.
     *
     * @param A The matrix that the vector will be created from.
     * @param col The column index that the vector will be created from.
     * @return An mx1 matrix with the values from the column specified.
     * @throws IllegalArgumentException If the specified column index is outside the matrix's bounds.
     */
    public static Matrix vectorFromColumn (Matrix A, int col) {
        if (A.getCols() < col || col <= 0) {
            throw new IllegalArgumentException("The row index is of bounds");
        }

        return new Matrix(A.getCol(col));
    }

    /**
     * Constructs an nx1 matrix from the column specified. This method was added after I submitted my application to your institution.
     *
     * @param A The matrix that the vector will be created from.
     * @param row The row index that the vector will be created from.
     * @return An nx1 matrix with the values from the row specified.
     * @throws IllegalArgumentException If the specified row index is outside the matrix's bounds.
     */
    public static Matrix vectorFromRow (Matrix A, int row) {
        if (A.getRows() < row || row <= 0) {
            throw new IllegalArgumentException("Row " + row + " does not exist in a matrix of " + A.getRows() + " row(s)!");
        }

        Matrix vector = new Matrix(A.getCols(), 1);

        for (int i = 1; i <= A.getCols(); i++) {
            vector.setValue(i, 1, A.getValue(row, i));
        }

        return vector;
    }

    //Multiplies the given mxn matrix by an nx1 vector.
    private static Matrix matrixVectorMultiplication (Matrix A, Matrix b) {
        if (b.getCols() != 1 || A.getCols() != b.getRows()) {
            throw new IllegalArgumentException();
        }

        Matrix newMatrix = new Matrix(A.getRows(), 1);

        for (int i = 1; i <= A.getRows(); i++) {
            double val = 0.0;
            for (int j = 1; j <= b.getRows(); j++) {
                val += A.getValue(i, j) * b.getValue(j, 1);
            }
            newMatrix.setValue(i, 1, val);
        }

        return newMatrix;
    }

    /**
     * Row-reduces the provided matrix to row reduced echelon form via Gaussian Elimination.
     *
     * @param matrixOrg The matrix that will be row reduced.
     * @return The matrix that is a copy of A that has been reduced to reduced row echelon form.
     */
    public static Matrix RREF(Matrix matrixOrg) {
        Matrix A = gaussianElimination(matrixOrg);
        A = backSolve(A);
        return A;
    }

    /**
     * Row-reduces the provided matrix to row reduced echelon form via Gaussian Elimination.
     *
     * @param matrixOrg The matrix that will be row reduced.
     * @return The vector that is the solution to the augmented system.
     */
    public static Matrix RREFSolve(Matrix matrixOrg) {
        Matrix A = gaussianElimination(matrixOrg);
        A = backSolve(A);

        Matrix x = vectorFromColumn(A, A.getCols());

        return x;
    }

    /**
     * Row-reduces the provided matrix to row echelon form via Gaussian Elimination.
     *
     * @param matrixOrg The matrix that will be row reduced.
     * @return The matrix that is a copy of A that has been reduced to row echelon form.
     */
    public static Matrix gaussianElimination(Matrix matrixOrg) {
        Matrix A;

        try {
            A = (Matrix) matrixOrg.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        A.setAugmentation(true);

        for (int i = 1; i <= A.getRows(); i++) { //i = current row
            partialPivoting(A, i); //Partially pivots the A.
            double[] vals = findPivot(A, i);
            double pivot = vals[0];
            int pivotCol = (int) vals[1];
            for (int j = i; j < A.getRows(); j++) { //Getting 0's bellow the pivot
                double ratio = -A.getValue(j + 1, pivotCol)/pivot; //S = -R1/R2 at the pivot column
                A.addRows(j + 1, i, ratio); //R1 = R1 + SR2
            }
        }

        return A;
    }

    //Pivots a matrix using scaled partial pivoting and bubble sort.
    private static void partialPivoting(Matrix A, int row) {
        for ( int i = row - 1; i < A.getRows(); i++ ) {
            for (int j = row; j < A.getRows() - i; j++) {
                if (A.compareScaledRows(j, j+1) < 0) {
                    A.swapRows(j, j+1);
                }
            }
        }
    }

    //Finds the pivot value in a row along with its column index, returned as [pivot value, col index].
    private static double[] findPivot (Matrix A, int row) {
        double pivot = A.getValue(row, 1);
        int i = 2;

        while (pivot == 0 && i <= A.getCols()) {
            pivot = A.getValue(row, i);
            i++;
        }

        return new double[]{pivot, --i};
    }

    //Finds the column index of a pivot value.
    private static int findPivotCol (Matrix A, int row) {
        double pivot = A.getValue(row, 1);
        int i = 1;

        while (pivot == 0 && i < A.getCols()) {
            pivot = A.getValue(row, i);
            if (pivot == 0 && i <= A.getCols()) {
                i++;
            }
        }

        return i;
    }

    /**
     * Back-solves an upper-triangular matrix, converting it to row reduced echelon form.
     *
     * @param matrixOrg The upper-triangular matrix that will be back-solved.
     * @return The matrix that is the row reduced echelon form version of the one provided.
     * @throws IllegalArgumentException If the provided matrix is not upper-triangular.
     */
    public static Matrix backSolve(Matrix matrixOrg) {
        if (!isZero(lowerTriSum(matrixOrg))) {
            throw new IllegalArgumentException("The provided matrix is not upper-triangular!");
        }

        Matrix A;

        try {
            A = (Matrix) matrixOrg.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        for (int i = A.getRows(); i >= 1; i--) { // Gets zero's above the current pivot
            int pivotCol = findPivotCol(A, i);
            double pivot = A.getValue(i, pivotCol);
            if (!isZero(pivot)) {
                A.scaleRow(i, 1 / pivot); //Scales the current row so the pivot = 1.0.
            }
            for (int j = i; j > 1; j--) {
                A.addRows(j - 1, i, -A.getValue(j - 1, pivotCol));
            }
        }

        return A;
    }

    //Calculates the sum of the lower-triangular values in a matrix.
    private static double lowerTriSum(Matrix A) {
        double sum = 0.0;

        for (int i = 1; i < A.getCols(); i++) {
            for (int j = i + 1; j <= A.getRows(); j++) {
                sum += Math.abs(A.getValue(j, i));
            }
        }

        return sum;
    }

    /**
     * Forward-solves an upper-triangular matrix, converting it to row reduced echelon form.
     *
     * @param matrixOrg The upper-triangular matrix that will be forward-solved.
     * @return The matrix that is the row reduced echelon form version of the one provided.
     * @throws IllegalArgumentException If the provided matrix is not lower-triangular.
     */
    public static Matrix forwardSolve (Matrix matrixOrg) {
        if (!isZero(upperTriSum(matrixOrg))) {
            throw new IllegalArgumentException("The provided matrix is not lower-triangular!");
        }

        Matrix A;

        try {
            A = (Matrix) matrixOrg.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        for (int i = 1; i < A.getRows(); i++) { // Gets zero's bellow the current pivot
            int pivotCol = findPivotCol(A, i);
            double pivot = A.getValue(i, pivotCol);

            if (!isZero(pivot)) {
                A.scaleRow(i, 1 / pivot); //Scales the current row so the pivot = 1.0.
            }

            for (int j = i; j < A.getRows(); j++) {
                A.addRows(j + 1, i, -A.getValue(j+1, pivotCol));
            }
        }

        return A;
    }

    //Calculates the sum of the upper-triangular values in a matrix.
    private static double upperTriSum(Matrix A) {
        double sum = 0.0;

        for (int i = 1; i <= A.getRows(); i++) {
            for (int j = i + 1; j <= A.getCols(); j++) {
                sum += Math.abs(A.getValue(i, j));
            }
        }

        return sum;
    }

    //The absolute sum of the off-diagonal terms of A.
    private static double nonDiagSum(Matrix A) {
        double sum = 0.0;

        for (int i = 1; i <= A.getCols(); i++) {
            for (int j = 1; j <= A.getRows(); j++) {
                if (i != j) {
                    sum += Math.abs(A.getValue(i, j));
                }
            }
        }

        return sum;
    }

    /**
     * Factors a matrix into permutation, lower-triangular, and upper-triangular matrices. Matrix P * matrixOrg =
     * matrix L * matrix U.
     *
     * @param matrixOrg The matrix that will be factored.
     * @return An array of matrices containing:
     *         <ul>
     *             <li> A lower-triangular matrix</li>
     *             <li> An upper-triangular matrix</li>
     *             <li> A permutation matrix which the provided one must be multiplied by for correct factorization</li>
     *         </ul>
     * @throws IllegalArgumentException If the provided matrix is not square.
     */
    public static Matrix[] LUDecomp(Matrix matrixOrg) { //Returns L and U
        if (!matrixOrg.isSquare()) {
            throw new IllegalArgumentException("The provided matrix must be square!");
        }

        Matrix U;
        Matrix P;
        int h = 1;
        int cols = matrixOrg.getCols();

        try {
            U = (Matrix) matrixOrg.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        Matrix L = identityMatrix(U.getRows());

        U = augmentMatrix(U, identityMatrix(U.getRows())); //Setup to get P
        partialPivoting(U, 1);

        try {
            P = (Matrix) U.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        for (int i = cols + 1; i <= 2 * cols; i++) { //Removing P from UP
            U = createSubMatrix(U, cols + 1);
        }

        U.setAugmentation(false);

        for (int j = cols; j >= 1; j--) { //Removing U from UP
            P = createSubMatrix(P, 1);
        }

        for (int i = 1; i <= U.getRows(); i++) { //Gaussian Elimination but L is created alongside U
            double[] vals = findPivot(U, i);
            double pivot = vals[0];
            int pivotCol = (int) vals[1];
            for (int k = h; k <= U.getRows(); k++) { //The columns of L will be the scaled pivot columns of U
                L.setValue(k, i, U.getValue(k, pivotCol)/pivot); //L[k,i] = U[k,i] / pivot
            }
            h++;
            for (int j = i; j < U.getRows(); j++) { //Getting 0's bellow the pivot
                double ratio = -U.getValue(j + 1, pivotCol)/pivot; //S = -R1/R2 at the pivot column
                U.addRows(j + 1, i, ratio); //R1 = R1 + SR2
            }
        }

        return new Matrix[]{L,U, P};
    }

    //Creates a copy of the provided matrix, with one column removed.
    private static Matrix createSubMatrix(Matrix A, int col) {
        Matrix subMatrix;

        try {
            subMatrix = (Matrix) A.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        subMatrix.removeCol(col);

        return subMatrix;
    }

    /**
     * Solves the system provided (matrixOrg * x = b) through forward and back solving its LU matrices.
     *
     * @param matrixOrg The matrix that will be solved for based on the provided vector.
     * @param b The matrix vector for which the system will be solved for.
     * @return A matrix vector which is the solution to the system.
     * @throws IllegalArgumentException If the provided matrix is not square.
     * @throws IllegalArgumentException If the provided matrix vector's number of columns is not one or if the
     * vector's number of rows is not equal to the matrix's number of columns.
     * @see #LUDecomp(Matrix)
     */
    public static Matrix LUSolve(Matrix matrixOrg, Matrix b) {
        if (!matrixOrg.isSquare()) {
            throw new IllegalArgumentException("The provided matrix must be square!");
        }

        if (b.getCols() != 1 || b.getRows() != matrixOrg.getCols()) {
            throw new IllegalArgumentException("The provided b matrix is not valid!");
        }

        matrixOrg.setAugmentation(false);

        Matrix[] LU = LUDecomp(matrixOrg);
        Matrix L = LU[0];
        Matrix U = LU[1];
        Matrix P = LU[2];

        P = matrixInverse(P);
        L = multiplyMatrices(P, L);

        Matrix Ly = RREF(augmentMatrix(L,b)); //Ly = b; Needed because P^-1 L is not triangular
        Matrix y = vectorFromColumn(Ly, Ly.getCols()); //Separates the solution from the matrix

        Matrix Ux = backSolve(augmentMatrix(U, y)); //Ux=y
        Ux.setAugmentation(false);
        partialPivoting(Ux, 1); //Gets the vector in order
        y = vectorFromColumn(Ux, Ux.getCols());

        return  y;
    }

    /**
     * Inverts the provided matrix, given that it is invertible.
     *
     * @param matrixOrg The square matrix that will be inverted.
     * @return The inverse of the matrix provided.
     * @throws IllegalArgumentException If the provided matrix is not square.
     * @throws IllegalArgumentException If the provided matrix is singular.
     */
    public static Matrix matrixInverse (Matrix matrixOrg) {
        if (!matrixOrg.isSquare()) {
            throw new IllegalArgumentException("The A must be square!");
        }

        if (isZero(determinant(matrixOrg))) {
            throw new IllegalArgumentException("The provided matrix is singular!");
        }

        Matrix A;

        try {
            A = (Matrix) matrixOrg.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        int cols = A.getCols();


        Matrix I = identityMatrix(A.getRows());
        A = augmentMatrix(A, I);

        A = RREF(A); //[A I] ~ [I A^-1]

        for (int i = 1; i <= cols; i++) { //Making [I A^-1] into [A^-1]
            A.removeCol(1);
        }

        return A;
    }

    /**
     * Right "divides" a provided matrix by another. Matrix quotient = matrix A * matrix B^(-1)
     *
     * @param A The matrix that will be "divided" by the other on the right.
     * @param B The matrix that will be "divided" by the other on the left.
     * @return The matrix that is the quotient of matrix A being right divided by matrix B.
     */
    public static Matrix divideMatrices(Matrix A, Matrix B) {
        return LinearAlgebra.multiplyMatrices(A, matrixInverse(B));
    }

    /**
     * Calculates the determinant of a square matrix.
     *
     * @param A The square matrix whose determinant will be calculated.
     * @return The determinant of the provided matrix.
     * @throws IllegalArgumentException If the provided matrix is not square.
     */
    public static double determinant(Matrix A) {
        if (!A.isSquare()) {
            throw new IllegalArgumentException("Matrix is not square");
        }

        return determinantRec(A);
    }

    //Recursively calculates the determinant of matrix A.
    private static double determinantRec(Matrix A) {
        double determinant = 0.0;

        if (A.getRows() == 2) {
            determinant = A.getValue(1, 1) * A.getValue(2, 2) - A.getValue(1, 2) * A.getValue(2, 1);
        }
        else {
            for (int i = 1; i <= A.getCols(); i++) {
                determinant += A.getValue(1, i) * Math.pow(-1, i - 1) * determinantRec(createSubMatrix(A, 1, i));
            }
        }

        return determinant;
    }

    //Creates a copy of the provided matrix, with one row and one column removed.
    private static Matrix createSubMatrix(Matrix A, int row, int col) {
        if (row <= 0 || col <= 0 || row > A.getRows() || col > A.getCols()) {
            throw new IllegalArgumentException("The row, column, or both are out of bounds!");
        }

        Matrix subMatrix;

        try {
            subMatrix = (Matrix) A.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        subMatrix.removeRow(row);
        subMatrix.removeCol(col);

        return subMatrix;
    }

    /**
     * Calculates the L1 norm of a matrix.
     *
     * @param A The matrix whose L1 norm will be calculated.
     * @return One of the following:
     *         <ul>
     *             <li> The sum of the matrix's absolute values if it has one column</li>
     *             <li> The largest absolute column sum of the matrix otherwise</li>
     *         </ul>
     */
    public static double l1Norm(Matrix A) {
        double l;

        if (A.getCols() == 1) {
            l = l1Vector(A, 1);
        }
        else {
            l = l1Matrix(A);
        }

        return l;
    }

    //The sum of the matrix's absolute values.
    private static double l1Vector (Matrix b, int col) {
        double l = 0.0;

        for (int i = 1; i <= b.getRows(); i++) {
            l += Math.abs(b.getValue(i, col));
        }

        return l;
    }

    //The largest absolute column sum of the matrix.
    private static double l1Matrix (Matrix A) {
        double l = 0.0;
        double temp;

        for (int i = 1; i <= A.getCols(); i++) {
            temp = l1Vector(A, i);
            if (temp > l) {
                l = temp;
            }
        }

        return l;
    }

    /**
     * Calculates the L2 norm of a matrix.
     *
     * @param A The matrix whose L2 norm will be calculated.
     * @return One of the following:
     *         <ul>
     *             <li> The length of the matrix if it has one column</li>
     *             <li> The square root of the matrix's spectral radius otherwise</li>
     *         </ul>
     * @throws UnsupportedOperationException If the provided matrix is not a vector nor square.
     */
    public static double l2Norm (Matrix A) {
        if (A.getCols() == 1) {
            return l2Vector(A);
        }
        else if (A.isSquare()) {
            return l2Matrix(A);
        }
        else {
            throw new UnsupportedOperationException("L2Norm of non-square matrices is not yet supported!");
        }
    }

    //The length of a vector.
    private static double l2Vector (Matrix b) {
        double l = 0.0;

        for (int i = 1; i <= b.getRows(); i++) {
            l += Math.pow(b.getValue(i, 1), 2);
        }

        return Math.sqrt(l);
    }

    //The square root of the matrix's spectral radius.
    private static double l2Matrix (Matrix A) {
        return Math.sqrt(maxEig(A));
    }

    /**
     * Calculates the LInfinite norm of a matrix.
     *
     * @param A The matrix whose LInfinite norm will be calculated.
     * @return One of the following:
     *         <ul>
     *             <li> The largest absolute value in the matrix if it has one column</li>
     *             <li> The largest absolute row sum of the matrix otherwise</li>
     *         </ul>
     */
    public static double lInfNorm (Matrix A) {
        if (A.getCols() == 1) {
            return lInfVector(A);
        }
        else {
            return lInfMatrix(A);
        }
    }

    //The largest absolute value in the vector.
    private static double lInfVector (Matrix b) {
        double l = 0.0;
        double temp;

        for (int i = 1; i <= b.getRows(); i++) {
            temp = Math.abs(b.getValue(i, 1));
            if (temp > l) {
                l = temp;
            }
        }

        return l;
    }

    //The largest absolute row sum of the matrix.
    private static double lInfMatrix (Matrix A) {
        double l = 0.0;
        double temp;

        for (int i = 1; i <= A.getRows(); i++) {
            temp = 0.0;
            for (int j = 1; j <= A.getCols(); j++) {
                temp += Math.abs(A.getValue(i, j));
            }
            if (temp > l) {
                l = temp;
            }
        }

        return l;
    }

    //The l2 norm of the residual error of the provided matrices (Ax - b).
    private static double rError2 (Matrix A, Matrix x, Matrix b) {
        return l2Norm(addMatrices(b, multiplyMatrices(A, x), -1));
    }

    //The l1 norm of the difference between the two provided vectors.
    private static double xError1 (Matrix A, Matrix B) {
        return l1Norm(addMatrices(A, B, -1));
    }

    //The l2 norm of the difference between the two provided vectors.
    private static double xError2 (Matrix x, Matrix b) {
        return l2Norm(addMatrices(x, b, -1));
    }

    /**
     * Normalizes a matrix such that each column has an L2 norm of 1.
     *
     * @param A The matrix that will be normalized.
     * @return The matrix with normalized columns.
     */
    public static Matrix normalize(Matrix A) {
        Matrix ci;

        for (int i = 1; i <= A.getCols(); i++) {
            ci = normalizeVector(vectorFromColumn(A, i));
            A.setCol(i, ci.getMatrix());
        }

        return A;
    }

    private static Matrix normalizeVector (Matrix b) {
        double l2 = l2Norm(b);

        for (int i = 1; i <= b.getRows(); i++) {
            b.scaleRow(i, 1.0/l2);
        }

        return b;
    }

    /**
     * Calculates the dot product of two vectors.
     *
     * @param a The first vector.
     * @param b The second vector.
     * @return The dot product of the two provided vectors.
     * @throws IllegalArgumentException If the provided matrices are not vectors.
     * @throws IllegalArgumentException If the provided matrices do not have the same number of rows
     */
    public static double dotProduct (Matrix a, Matrix b) {
        if (a.getCols() != 1 || b.getCols() != 1) {
            throw new IllegalArgumentException("The matrices must be vectors!");
        }

        if (a.getRows() != b.getRows()) {
            throw new IllegalArgumentException("The vectors do not have the same number of rows!");
        }

        double dp = 0.0;

        for (int i = 1; i <= a.getRows(); i++) {
            dp += a.getValue(i, 1) * b.getValue(i, 1);
        }

        return dp;
    }

    /**
     * Calculates the cross product of two vectors.
     *
     * @param a The first vector.
     * @param b The second vector.
     * @return The cross product of the two provided vectors.
     * @throws IllegalArgumentException If the provided matrices are not vectors.
     * @throws IllegalArgumentException If the provided matrices do not have the same number of rows.
     * @throws IllegalArgumentException If the provided matrices have less than two rows.
     */
    public static Matrix crossProduct (Matrix a, Matrix b) {
        if (a.getCols() != 1 || b.getCols() != 1) {
            throw new IllegalArgumentException("The matrices must be vectors!");
        }

        if (a.getRows() != b.getRows()) {
            throw new IllegalArgumentException("The vectors do not have the same number of rows!");
        }

        if (a.getRows() < 2) {
            throw new IllegalArgumentException("The provided object must be vectors, not scalars!");
        }

        if (a.getRows() == 2) { //Converting the 2d vectors into 3d ones
            a.addRowBottom(new double[] {0.0});
            b.addRowBottom(new double[] {0.0});
        }

        Matrix cp = new Matrix(a.getRows(), 1);
        Matrix detMatrix = transpose(a);
        detMatrix.addRowBottom(b.getMatrix());
        double det;

        for (int i = 1; i <= a.getRows(); i++) {
            det = Math.pow(-1.0, i + 1) * determinant(createSubMatrix(detMatrix, i));
            cp.setValue(i, 1, det);
        }

        return cp;
    }

    /**
     * Calculates the projection of y onto u. (Y dot U / U dot U) * U
     *
     * @param y The vector to be projected.
     * @param u The vector to be projected onto.
     * @return The projection of y onto u.
     * @throws IllegalArgumentException If the provided matrices are not vectors.
     * @throws IllegalArgumentException If the provided matrices do not have the same number of rows
     */
    public static Matrix proj(Matrix y, Matrix u) {
        if (y.getCols() != 1 || u.getCols() != 1) {
            throw new IllegalArgumentException("The matrices must be vectors!");
        }

        if (y.getRows() != u.getRows()) {
            throw new IllegalArgumentException("The vectors not have the same number of rows!");
        }

        Matrix proj;

        try {
            proj = (Matrix) u.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        double scalar = dotProduct(y, u) / dotProduct(u, u);

        for (int i = 1; i <= u.getRows(); i++) {
            proj.scaleRow(i, scalar);
        }

        return proj;
    }

    /**
     * Creates an orthogonal basis of vectors in R^n from the provided mxn matrix.
     *
     * @param matrixOrg The set of n vectors to be used to form an orthogonal basis in R^n.
     * @return A set of orthogonal basis vectors in R^n in the form a matrix.
     */
    public static Matrix GramSchmidt (Matrix matrixOrg) {
        Matrix orthMatrix;
        Matrix vn;
        Matrix ci;
        Matrix xn;
        Matrix temp;

        try {
            orthMatrix = (Matrix) matrixOrg.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        for (int i = 2; i <= orthMatrix.getCols(); i++) { //Start at 2 b/c v1=x1
            temp = zeroMatrix(matrixOrg.getRows(), 1);
            xn = vectorFromColumn(matrixOrg, i);

            for (int j = 1; j < i; j++) {
                ci = vectorFromColumn(orthMatrix, j);
                temp = addMatrices(temp, proj(xn, ci), 1);
            }
            vn = addMatrices(xn, temp, -1 );
            orthMatrix.setCol(i, vn.getMatrix());
        }

        return orthMatrix;
    }

    /**
     * Numerically calculates the dominant eigenvalue of a square matrix at a linear rate via the power method, provided
     * the matrix has a unique dominant eigenvalue.
     *
     * @param A The square matrix whose dominant eigenvalue will be calculated.
     * @param x The initial guess for the matrix's dominant eigenvector.
     * @param t The tolerance this method will use.
     * @return The dominant eigenvalue of the provided matrix.
     * @throws IllegalArgumentException If the provided matrix is not square.
     */
    public static double powerMethod (Matrix A, Matrix x, double t) {
        if (!A.isSquare()) {
            throw new IllegalArgumentException("The matrix must be square!");
        }

        Matrix xk = normalize(multiplyMatrices(A, x)); //xk = Ax(k-1)/||x(k-1)||
        int iterations = 0;

        while (xError2(xk, x) > t && iterations < MAX_ITERATIONS) {
            x = xk;
            xk = normalize(multiplyMatrices(A, xk)); //xk = Ax(k-1)/||x(k-1)||
            iterations++;
        }

        Matrix eigVal = multiplyMatrices(A, xk); //eig = xTk*A*xk
        eigVal = multiplyMatrices(transpose(xk), eigVal);

        return eigVal.getValue(1, 1);
    }

    /**
     * Numerically calculates the dominant eigenvalue of a square matrix at a linear rate via the power method, provided
     * the matrix has a unique dominant eigenvalue.
     *
     * @param A The square matrix whose dominant eigenvalue will be calculated.
     * @return The dominant eigenvalue of the provided matrix.
     * @throws IllegalArgumentException If the provided matrix is not square.
     */
    public static double maxEig (Matrix A) {
        if (!A.isSquare()) {
            throw new IllegalArgumentException("The matrix must be square!");
        }

        Matrix x = new Matrix(A.getRows(), 1);

        for (int i = 1; i <= x.getRows(); i++) {
            x.setValue(i,1, 1);
        }

        return powerMethod(A, x, TOL);
    }

    /**
     * Factors a matrix into the product of an orthogonal and upper-triangular matrix. Matrix A = matrix Q * matrix R.
     *
     * @param A The matrix that will be factored.
     * @return An array of matrices containing:
     *         <ul>
     *             <li> An orthogonal matrix</li>
     *             <li> An upper-triangular matrix</li>
     *         </ul>
     */
    public static Matrix[] QRFactorization (Matrix A) {
        Matrix[] QR = new Matrix[2];

        Matrix Q = normalize(GramSchmidt(A));
        Matrix R = multiplyMatrices(transpose(Q), A);
        QR[0] = Q;
        QR[1] = R;

        return QR;
    }

    /**
     * Numerically calculates the eigenvalues of a square matrix via the QR method, provided all eigenvalues are distinct
     * and real.
     *
     * @param A The square matrix whose eigenvalues will be calculated.
     * @param t The tolerance this method will use.
     * @return A matrix with eigenvalues along its diagonal.
     * @throws IllegalArgumentException If the provided matrix is not square.
     */
    public static Matrix QREig (Matrix A, double t) {
        if (A.getRows() != A.getCols()) {
            throw new IllegalArgumentException("The matrix must be square!");
        }

        Matrix[] QR = QRFactorization(A);
        Matrix Q = QR[0];
        Matrix R = QR[1];
        Matrix eig = multiplyMatrices(R, Q);
        int iterations = 0;

        while (lowerTriSum(eig) > t && iterations < MAX_ITERATIONS) {
            QR = QRFactorization(eig);
            Q = QR[0];
            R = QR[1];
            eig = multiplyMatrices(R, Q); //E(k) = R(k)Q(k)
            iterations++;
        }

        return eig;
    }

    /**
     * Numerically calculates the eigenvalues of a square matrix via the QR shifting method, provided all eigenvalues
     * are distinct and real. The eigenvalue corresponding to the provided guess will be the most accurate, with all
     * others being less accurate.
     *
     * @param A The square matrix whose eigenvalues will be calculated.
     * @param row The row coordinate of the guess for the desired eigenvalue.
     * @param col The column coordinate of the guess for the desired eigenvalue.
     * @param t The tolerance this method will use.
     * @return A matrix with eigenvalues along its diagonal.
     * @throws IllegalArgumentException If the provided matrix is not square.
     */
    public static Matrix QREigShift (Matrix A, int  row,int col,double t) {
        if (A.getRows() != A.getCols()) {
            throw new IllegalArgumentException("The matrix must be square!");
        }

        double shift = A.getValue(row, col);
        A = addMatrices(A, constantIdentityMatrix(A.getRows(), shift), -1); //~A(k) = A(k) - shift*I
        Matrix[] QR = QRFactorization(A);
        Matrix Q = QR[0];
        Matrix R = QR[1];
        Matrix eig = addMatrices(multiplyMatrices(R, Q), constantIdentityMatrix(A.getRows(), shift), 1); //E1 = RQ + shift*I
        int iterations = 0;

        while (lowerTriSum(eig) > t && iterations < MAX_ITERATIONS) {
            shift = eig.getValue(row, col);
            eig = addMatrices(eig, constantIdentityMatrix(A.getRows(), shift), -1); //~E(k) = E(k) - shift*I
            QR = QRFactorization(eig);
            Q = QR[0];
            R = QR[1];
            eig = addMatrices(multiplyMatrices(R, Q), constantIdentityMatrix(A.getRows(), shift), 1);  //E(k) = R(k)Q(k) + shift*I
            iterations++;
        }

        return eig;
    }

    /**
     * Numerically calculates the eigenvalues of a square matrix via the QR method, provided all eigenvalues are distinct
     * and real.
     * <p>
     *     The eigenvalues are calculated using the regular QR method due to its output of all of a matrix's eigenvalues
     *     with similar accuracy.
     * </p>
     *
     * @param A The square matrix whose eigenvalues will be calculated.
     * @return An array of the matrix's eigenvalues.
     * @throws IllegalArgumentException If the provided matrix is not square.
     * @see #QREig(Matrix, double)
     */
    public static double[] eig(Matrix A) {
        if (!A.isSquare()) {
            throw new IllegalArgumentException("The matrix must be square!");
        }

        Matrix eigenVals =  QREig(A, TOL);

        return diagVals(eigenVals);
    }

    //Creates a 1d array of the values along a matrix's diagonal.
    private static double[] diagVals(Matrix A) {
        if (!(A.isSquare())) {
            throw new IllegalArgumentException("The matrix must be square!");
        }

        double[] vals = new double[A.getCols()];

        for (int i = 1; i <= A.getRows(); i++) {
            vals[i-1] = A.getValue(i, i);
        }

        return vals;
    }
}