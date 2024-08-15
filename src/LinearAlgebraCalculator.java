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

    private static Matrix vectorFromColumn (Matrix tempMatrix, int col) {
        Matrix vector = new Matrix(tempMatrix.getRows(), 1);
        for (int i = 1; i <= tempMatrix.getRows(); i++) {
            vector.setValue(i, 1, tempMatrix.getValue(i, col));
        }
        return vector;
    }

    public static Matrix multiplyMatrices(Matrix matrix1, Matrix matrix2) {
        if (matrix1.getCols() != matrix2.getRows()) {
            throw new IllegalArgumentException("The second matrix must have a number of rows equal to the number of columns of the first! " + matrix1.getRows() + " != " +  matrix2.getRows() + "!");
        }
        Matrix product = new Matrix(matrix1.getRows(), matrix2.getCols());
        for (int i = 1; i <= matrix2.getCols(); i++) {
            Matrix tempVector = vectorFromColumn(matrix2, i);
            Matrix tempCol = matrixVectorMultiplication(matrix1, tempVector);
            for (int j = 1; j <= product.getRows(); j++) {
                product.setValue(j, i, tempCol.getValue(j, 1));
            }
        }
        return product;
    }

    private static Matrix matrixVectorMultiplication (Matrix matrix, Matrix vector) {
        if (vector.getCols() != 1 || matrix.getCols() != vector.getRows()) {
            throw new IllegalArgumentException();
        }
        Matrix newMatrix = new Matrix(matrix.getRows(), 1);
        for (int i = 1; i <= matrix.getRows(); i++) {
            double val = 0.0;
            for (int j = 1; j <= vector.getRows(); j++) {
                val += matrix.getValue(i, j) * vector.getValue(j, 1);
            }
            newMatrix.setValue(i, 1, val);
            val = 0.0;
        }
        return newMatrix;
    }
}
