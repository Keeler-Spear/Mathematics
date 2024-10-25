public class LinearAlgebraCalculator {

    //Determinant and helper methods
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

    private static Matrix identityMatrix(int n) {
        double[][] matrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matrix[i][j] = 0.0;
                if (i == j) {
                    matrix[i][j] = 1;
                }
            }
        }
        return new Matrix(matrix);
    }

    //Matrix Multiplication and helper methods
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



    //Gaussian Elimination, LU decomposition, and helper methods
    public static Matrix augmentMatrix(Matrix matrix1, Matrix matrix2) {
        if (matrix1.getRows() != matrix2.getRows()) {
            throw new IllegalArgumentException("The matrices must have the same number of rows!");
        }
        Matrix matrix = new Matrix (matrix1.getRows(), matrix1.getCols() + matrix2.getCols());
        for (int i = 1; i <= matrix1.getCols(); i++) {
            for (int j = 1; j <= matrix1.getRows(); j++) {
                matrix.setValue(j, i, matrix1.getValue(j, i));
            }
        }
        for (int i = 1; i <= matrix2.getCols(); i++) {
            for (int j = 1; j <= matrix1.getRows(); j++) {
                matrix.setValue(j, i + matrix1.getCols(), matrix2.getValue(j, i));
            }
        }
        return matrix;
    }

    private static void partialPivoting(Matrix matrix, int row) {//Can optimize further if I use quick sort
        for ( int i = row - 1; i < matrix.getRows(); i++ ) {
            for (int j = row; j < matrix.getRows(); j++) {
                if (matrix.compareRows(j, j+1, 1) < 0) {
                    matrix.swapRows(j, j+1);
                }
            }
        }
    }

    private static double[] findPivot (Matrix matrix, int row) {
        double pivot = matrix.getValue(row, 1);
        int i = 2;
        while (pivot == 0 && i <= matrix.getCols()) {
            pivot = matrix.getValue(row, i);
            i++;
        }
        return new double[]{pivot, --i};
    }

    private static int findPivotCol (Matrix matrix, int row) {
        double pivot = matrix.getValue(row, 1);
        int i = 1;
        while (pivot == 0 && i < matrix.getCols()) {
            pivot = matrix.getValue(row, i);
            if (pivot == 0 && i <= matrix.getCols()) {
                i++;
            }
        }
        return i;
    }

    public static Matrix gaussianElimination(Matrix matrixOrg) {
        Matrix matrix;
        try { //Needs to be used so the code in not altering the original matrix A
            matrix = (Matrix) matrixOrg.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        for (int i = 1; i <= matrix.getRows(); i++) { //i = current row
            partialPivoting(matrix, i); //Partially pivots the matrix. While it includes prior rows, the 0's in the columns below the pivot positions prevents any jankyness.
            double[] vals = findPivot(matrix, i);
            double pivot = vals[0];
            int pivotCol = (int) vals[1];


            for (int j = i; j < matrix.getRows(); j++) {//Getting 0's bellow the pivot
                double ratio = -matrix.getValue(j + 1, pivotCol)/pivot; //S = -R1/R2 at the pivot column
                matrix.addRows(j + 1, i, ratio); //R1 = R1 + SR2
            }
        }
        return matrix;
    }

    private static Matrix backSolve(Matrix matrixOrg) {
        Matrix matrix;
        try {
            matrix = (Matrix) matrixOrg.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        for (int i = matrix.getRows(); i >= 1; i--) {// Gets zero's above the current pivot
            int pivotCol = findPivotCol(matrix, i);
            double pivot = matrix.getValue(i, pivotCol);
            if (pivot != 0) {
                matrix.scaleRow(i, 1 / pivot); //Scales the current row so the pivot = 1.0.
            }
            for (int j = i; j > 1; j--) {
                matrix.addRows(j - 1, i, -matrix.getValue(j - 1, pivotCol));
            }
        }
        partialPivoting(matrix, 1);
        return matrix;
    }

    private static Matrix forwardSolve (Matrix matrixOrg) {
        Matrix matrix;
        try {
            matrix = (Matrix) matrixOrg.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        for (int i = 1; i < matrix.getRows(); i++) {// Gets zero's bellow the current pivot
            int pivotCol = findPivotCol(matrix, i);
            double pivot = matrix.getValue(i, pivotCol);
            if (pivot != 0) {
                matrix.scaleRow(i, 1 / pivot); //Scales the current row so the pivot = 1.0.
            }
            for (int j = i; j < matrix.getRows(); j++) {
                matrix.addRows(j + 1, i, -matrix.getValue(j+1, pivotCol));
            }
        }
        return matrix;
    }

    public static Matrix RREF(Matrix matrixOrg) {
        Matrix matrix;
        matrix = gaussianElimination(matrixOrg); //Turns a matrix into row echelon form
        System.out.println(matrix);
        matrix = backSolve(matrix); // backsolves for the solution.
        Matrix x = vectorFromColumn(matrix, matrix.getCols());
        return x;
    }

    public static Matrix matrixInverse (Matrix matrixOrg) { //Row Reduce {A, I}
        if (!matrixOrg.isSquare()) {
            throw new IllegalArgumentException("The matrix must be square!");
        }

        if (determinant(matrixOrg) == 0) { //Is A ~ I? via Invertible Matrix Theorem
            throw new IllegalArgumentException("The matrix is singular!");
        }

        Matrix matrix;

        try {
            matrix = (Matrix) matrixOrg.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        for (int i = 0; i < matrix.getRows(); i++) { //Making [A I]
            double[] e = new double [matrix.getRows()];
            for (int j = 0; j < matrix.getRows(); j++) {
                e[j] = 0.0;
                if (i == j) {
                    e[j] = 1.0;
                }
            }
            matrix.addCol(e);
        }

        matrix = RREF(matrix); //[A I] ~ [I A^-1]

        for (int i = 1; i <= matrix.getRows(); i++) {//Making [I A^-1] into [A^-1]
            matrix.removeCol(1);
        }

        return matrix;
    }

    public static Matrix[] LUDecomp(Matrix matrixOrg) { //Returns L and U
        Matrix U; //U will be an mxn matrix

        try { //Needs to be used so the code in not altering the original matrix A
            U = (Matrix) matrixOrg.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        Matrix L = identityMatrix(U.getRows()); //L will be an mxm matrix


        //Gaussian Elimination but L is created alongside U
        for (int i = 1; i <= U.getRows(); i++) { //i = current row
            double[] vals = findPivot(U, i);
            double pivot = vals[0];
            int pivotCol = (int) vals[1];

            for (int k = pivotCol; k <= U.getRows(); k++) { //The columns of L will be the scaled pivot columns of U
                L.setValue(k, i, U.getValue(k, pivotCol)/pivot); //L[k,i] = U[k,i] / pivot
            }


            for (int j = i; j < U.getRows(); j++) {//Getting 0's bellow the pivot
                double ratio = -U.getValue(j + 1, pivotCol)/pivot; //S = -R1/R2 at the pivot column
                U.addRows(j + 1, i, ratio); //R1 = R1 + SR2
            }
        }

        return new Matrix[]{L,U}; //Returns a list of L and U
    } //Returns L and U

    public static Matrix LUSolve(Matrix matrixOrg, Matrix b) {
        Matrix[] LU = LUDecomp(matrixOrg);
        Matrix L = LU[0];
        Matrix U = LU[1];
        System.out.println("U");
        System.out.println(U);
        Matrix Ly = forwardSolve(augmentMatrix(L, b)); //Ly=b
        Matrix y = vectorFromColumn(Ly, Ly.getCols()); //Separates the solution from the matrix

        Matrix Ux = backSolve(augmentMatrix(U, y)); //Ux=y
        return vectorFromColumn(Ux, Ux.getCols()); //Separates the solution from the matrix
    }
}
