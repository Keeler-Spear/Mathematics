public class LinearAlgebraCalculator {

    public static double determinant(Matrix matrix) {
        if (!matrix.isSquare()) {
            throw new IllegalArgumentException("Matrix is not square");
        }
        return determinantRec(matrix);
    }

    private static double determinantRec(Matrix matrix) {
        double determinant = 0.0;
        if (matrix.getRows() == 2) { //Base Case
            determinant = matrix.getValue(1, 1) * matrix.getValue(2, 2) - matrix.getValue(1, 2) * matrix.getValue(2, 1);
        }
        else {
            for (int i = 1; i <= matrix.getCols(); i++) {
                determinant += matrix.getValue(1, i) * Math.pow(-1, i - 1) * determinantRec(createSubMatrix(matrix, 1, i));
            }
        }
        return determinant;
    }

    private static Matrix createSubMatrix(Matrix matrix, int row, int col) {
        Matrix subMatrix;
        try {
            subMatrix = (Matrix) matrix.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }
        subMatrix.removeRow(row);
        subMatrix.removeCol(col);
        return subMatrix;
    }

    //ToDo: Fix
//    private static Matrix vectorFromColumn (Matrix tempMatrix, int col) {
//        double[][] vector = new double[tempMatrix.getRows()][1];
//        for (int i = 0; i < tempMatrix.getRows(); i++) {
//            vector[i][0] = tempMatrix.getValue(i, col - 1);
//        }
//        return vector;
//    }

    //ToDo Finish multiplication method
    //Calle * Peram Matrices
//    public Matrix multiplyMatrices(Matrix matrix2) {
//        if (cols != matrix2.getRows()) {
//            throw new IllegalArgumentException("The second matrix must have a number of rows equal to the number of columns of the first! " + cols + " != " +  matrix2.getRows() + "!");
//        }
//    }

    //Todo Finish matrix-vector multiplication method
//    private Matrix matrixVectorMultiplication (Matrix matrix1, Matrix matrix2) {
//
//    }
}
