public class LinearAlgebraCalculator {

    final static double tol = 0.000001;

    public static Matrix addMatrices(Matrix A, Matrix B, double sign) {//1 for add, -1 for sub\
        //A + B or A - B
        if (A.getRows() != B.getRows() || A.getCols() != B.getCols()) {
            throw new IllegalArgumentException("The matrices must be of the same size!");
        }

        Matrix result = new Matrix (A.getRows(), A.getCols());

        for (int i = 1; i <= A.getRows(); i++) {
            for (int j = 1; j <= A.getCols(); j++) {
                result.setValue(i, j, A.getValue(i , j) + sign * B.getValue(i, j));
            }
        }
        return result;
    }
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

    private static Matrix createSubMatrix(Matrix matrix, int col) {
        Matrix subMatrix;
        try {
            subMatrix = (Matrix) matrix.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }
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

    private static Matrix constantIdentityMatrix(int n, double c) {
        double[][] matrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matrix[i][j] = 0.0;
                if (i == j) {
                    matrix[i][j] = c;
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
            for (int j = row; j < matrix.getRows() - i; j++) {
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
            partialPivoting(matrix, i); //Partially pivots the matrix.
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

    public static Matrix RREF(Matrix matrixOrg) { //Returns RREF(A)
        Matrix matrix;
        matrix = gaussianElimination(matrixOrg); //Turns a matrix into row echelon form
        partialPivoting(matrix, 1); //Pivots to ensure we can backsolve
        matrix = backSolve(matrix); // backsolves for the solution.
        return matrix;
    }

    public static Matrix RREFSolve(Matrix matrixOrg) { //Returns X
        Matrix matrix;
        matrix = gaussianElimination(matrixOrg); //Turns a matrix into row echelon form
        partialPivoting(matrix, 1); //Pivots to ensure we can backsolve
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

        int cols = matrix.getCols();


        Matrix I = identityMatrix(matrix.getRows());
        matrix = augmentMatrix(matrix, I);

        matrix = RREF(matrix); //[A I] ~ [I A^-1]

        for (int i = 1; i <= cols; i++) {//Making [I A^-1] into [A^-1]
            matrix.removeCol(1);
        }

        return matrix;
    }

    //TODO: Multiply by p^-1
    public static Matrix[] LUDecomp(Matrix matrixOrg) { //Returns L and U
        Matrix U; //U will be an mxn matrix
        Matrix P; //Permutation Matrix
        int h = 1;
        int cols = matrixOrg.getCols();

        try { //Needs to be used so the code in not altering the original matrix A
            U = (Matrix) matrixOrg.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        Matrix L = identityMatrix(U.getRows()); //L will be an mxm matrix
        U = augmentMatrix(U, identityMatrix(U.getRows())); //Setup to get P
        partialPivoting(U, 1); //Partially pivoting U for clean Gaussian Elimination

        try { //Needs to be used so the code in not altering the original matrix
            P = (Matrix) U.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        for (int i = cols + 1; i <= 2 * cols; i++) { //Removing P from UP
            U = createSubMatrix(U, cols + 1);
        }

        for (int j = cols; j >= 1; j--) {//Removing U from UP
            P = createSubMatrix(P, 1);
        }

        //Gaussian Elimination but L is created alongside U
        for (int i = 1; i <= U.getRows(); i++) { //i = current row
            double[] vals = findPivot(U, i);
            double pivot = vals[0];
            int pivotCol = (int) vals[1];

            for (int k = h; k <= U.getRows(); k++) { //The columns of L will be the scaled pivot columns of U. h is used to ensure L is lower-triangular.
                L.setValue(k, i, U.getValue(k, pivotCol)/pivot); //L[k,i] = U[k,i] / pivot
            }
            h++;

            for (int j = i; j < U.getRows(); j++) {//Getting 0's bellow the pivot
                double ratio = -U.getValue(j + 1, pivotCol)/pivot; //S = -R1/R2 at the pivot column
                U.addRows(j + 1, i, ratio); //R1 = R1 + SR2
            }
        }
        return new Matrix[]{L,U, P}; //Returns a list of L, U, and P
    } //Returns L and U

    public static Matrix LUSolve(Matrix matrixOrg, Matrix b) { //Returns to solution to the system
        Matrix[] LU = LUDecomp(matrixOrg);
        Matrix L = LU[0];
        Matrix U = LU[1];
        Matrix P = LU[2];
        P = matrixInverse(P);
        L = multiplyMatrices(P, L); //L = P^-1 L
        Matrix Ly = RREF(augmentMatrix(L,b)); //Ly = b; Needed because P^-1 L is not triangular
        Matrix y = vectorFromColumn(Ly, Ly.getCols()); //Separates the solution from the matrix

        Matrix Ux = backSolve(augmentMatrix(U, y)); //Ux=y
        partialPivoting(Ux, 1); //Gets the vector in order
        y = vectorFromColumn(Ux, Ux.getCols()); //Separates the solution from the matrix
        return  y;
    }

    //Eigenvalues

    public static Matrix transpose(Matrix matrixOrg) { //r1 = c1, r2 = c2
        Matrix matrix = new Matrix (matrixOrg.getCols(), matrixOrg.getRows());
            for (int i = 1; i <= matrixOrg.getCols(); i++) {
                for (int j = 1; j <= matrixOrg.getRows(); j++) {
                    matrix.setValue(i, j, matrixOrg.getValue(j, i));
                }
            }

        return matrix;
    }

    private static double l2Vector (Matrix matrix) {
        double l = 0.0;
        for (int i = 1; i <= matrix.getRows(); i++) {
            l += Math.pow(matrix.getValue(i, 1), 2);
        }
        return Math.sqrt(l);
    }

    //Todo
    private static double l2Matrix (Matrix matrix) {
        return 0.0;
    }

    public static double l2Norm (Matrix matrix) {
        if (matrix.getCols() == 1) {
            return l2Vector(matrix);
        }
        else {
            return l2Matrix(matrix);
        }
    }

    public static double l1Norm(Matrix matrix) {
        double a = 0.0;
        double b;
        for (int i = 1; i <= matrix.getCols(); i++) {
            b = 0.0;
            for (int j = 1; j <= matrix.getRows(); j++) {
                b += Math.abs(matrix.getValue(j, i));
            }
            if (b > a) {
                a = b;
            }
        }
        return a;
    }

    private static double xError1 (Matrix A, Matrix B) {
        return l1Norm(addMatrices(A, B, -1));
    }
    private static double rError2 (Matrix A, Matrix x, Matrix b) {
        return l2Norm(addMatrices(b, multiplyMatrices(A, x), -1));
    }
    private static double xError2 (Matrix x, Matrix b) {
        return l2Norm(addMatrices(x, b, -1));
    }

    public static Matrix normalize(Matrix x) {
        Matrix x1;
        for (int i = 1; i <= x.getCols(); i++) {
            x1 = normalizeVector(vectorFromColumn(x, i));
            x.setCol(i, x1.getMatrix());
        }
        return x;
    }

    private static Matrix normalizeVector (Matrix x) {
        double l2 = l2Norm(x);
        for (int i = 1; i <= x.getRows(); i++) {
            x.scaleRow(i, 1.0/l2);
        }
        return x;
    }

    public static double powerMethod (Matrix A, Matrix x, double t) { //Returns dominant eigenvalue and associated eigenvector
        if (A.getCols() != A.getRows()) {
            throw new IllegalArgumentException("The matrix must be square!");
        }
        Matrix xk = normalize(multiplyMatrices(A, x)); //xk = Ax(k-1)/||x(k-1)||
        while (xError2(xk, x) > t) {
            x = xk;
            xk = normalize(multiplyMatrices(A, xk)); //xk = Ax(k-1)/||x(k-1)||
        }
        Matrix eigVal = multiplyMatrices(A, xk); //eig = xTk*A*xk
        eigVal = multiplyMatrices(transpose(xk), eigVal);
        return eigVal.getValue(1, 1);
    }

    public static double powerMethod (Matrix A) {
        Matrix x = new Matrix(A.getRows(), 1);
        for (int i = 1; i <= x.getRows(); i++) {
            x.setValue(i,1, 1);
        }
        return powerMethod(A, x, tol);
    }

    //Vectors
    public static double dotProduct (Matrix a, Matrix b) {
        if (a.getCols() != 1 || b.getCols() != 1) {
            throw new IllegalArgumentException("The matrices must be vectors!");
        }
        if (a.getRows() != b.getRows()) {
            throw new IllegalArgumentException("The vectors not have the same number of rows!");
        }
        double dp = 0.0;

        for (int i = 1; i <= a.getRows(); i++) {
            dp += a.getValue(i, 1) * b.getValue(i, 1);
        }
        return dp;
    }

    public static Matrix proj(Matrix y, Matrix u) { //(Y*U)/(U*U) U
        if (y.getCols() != 1 || u.getCols() != 1) {
            throw new IllegalArgumentException("The matrices must be vectors!");
        }
        if (y.getRows() != u.getRows()) {
            throw new IllegalArgumentException("The vectors not have the same number of rows!");
        }
        Matrix proj;
        try { //Needs to be used so the code in not altering the original matrix A
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

    public static Matrix zeroMatrix (int rows, int cols) {
        Matrix matrix = new Matrix(rows, cols);
        for (int i = 1; i <= rows; i++) {
            for (int j = 1; j <= cols; j++) {
                matrix.setValue(i, j, 0);
            }
        }
        return matrix;
    }

    private static double[] diagVals(Matrix matrix) {
        if (matrix.getCols() != matrix.getRows()) {
            throw new IllegalArgumentException("The matrix must be square!");
        }
        double[] vals = new double[matrix.getCols()];
        for (int i = 1; i <= matrix.getRows(); i++) {
            vals[i-1] = matrix.getValue(i, i);
        }
        return vals;
    }

    public static Matrix GramSchmidt (Matrix  x) {
        Matrix v;
        Matrix vn;
        Matrix vnTemp;
        Matrix xn;
        Matrix temp;
        try { //Needs to be used so the code in not altering the original matrix A
             v = (Matrix) x.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        for (int i = 2; i <= v.getCols(); i++) { //Start at 2 b/c v1=x1
            temp = zeroMatrix(x.getRows(), 1);
            xn = vectorFromColumn(x, i);

            for (int j = 1; j < i; j++) {
                vnTemp = vectorFromColumn(v, j);
                temp = addMatrices(temp, proj(xn, vnTemp), 1);
            }

            vn = addMatrices(xn, temp, -1 );
            v.setCol(i, vn.getMatrix());
        }

        return v;
    }

    public static Matrix[] QRFactorization (Matrix A) {
        Matrix[] QR = new Matrix[2];
        Matrix Q = normalize(GramSchmidt(A));
        Matrix R = multiplyMatrices(transpose(Q), A);
        QR[0] = Q;
        QR[1] = R;
        return QR;
    }

    public static Matrix QREig (Matrix A, double t) {
        Matrix eig0 = A;
        Matrix[] QR = QRFactorization(A);
        Matrix Q = QR[0];
        Matrix R = QR[1];
        Matrix eig1 = multiplyMatrices(R, Q); //E1 = RQ
        while (xError1(eig1, eig0) > t) { //Runs while ||E1 - E0||1 > tolerance
            QR = QRFactorization(eig1);
            Q = QR[0];
            R = QR[1];
            eig0 = eig1;
            eig1 = multiplyMatrices(R, Q); //E1(k) = R(k)Q(k)
        }
        return eig1;
    }

    public static Matrix QREigShift (Matrix A, int  i,int j,double t) { //Shifting is based on element at (i, j)
        Matrix eig0 = A;
        double shift = A.getValue(i, j);
        A = addMatrices(A, constantIdentityMatrix(A.getRows(), shift), -1); //~A(k) = A(k) - shift*I
        Matrix[] QR = QRFactorization(A);
        Matrix Q = QR[0];
        Matrix R = QR[1];
        Matrix eig1 = addMatrices(multiplyMatrices(R, Q), constantIdentityMatrix(A.getRows(), shift), 1); //E1 = RQ + shift*I
        while (xError1(eig1, eig0) > t) { //Runs while ||E1 - E0||1 > tolerance
            eig0 = eig1;
            shift = eig1.getValue(i, j);
            eig1 = addMatrices(eig1, constantIdentityMatrix(A.getRows(), shift), -1); //~E1(k) = E1(k) - shift*I
            QR = QRFactorization(eig1);
            Q = QR[0];
            R = QR[1];
            eig1 = addMatrices(multiplyMatrices(R, Q), constantIdentityMatrix(A.getRows(), shift), 1);  //E1(k) = R(k)Q(k) + shift*I
        }
        return eig1;
    }


}