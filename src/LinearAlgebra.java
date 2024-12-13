import java.sql.SQLOutput;

public class LinearAlgebra {

    final static double tol = 0.000001;

    public static Matrix addMatrices(Matrix A, Matrix B, double sign) {//1 for add, -1 for sub
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
    public static double determinant(Matrix A) {
        if (!A.isSquare()) {
            throw new IllegalArgumentException("Matrix is not square");
        }
        return determinantRec(A);
    }

    private static double determinantRec(Matrix A) {
        double determinant = 0.0;
        if (A.getRows() == 2) { //Base Case
            determinant = A.getValue(1, 1) * A.getValue(2, 2) - A.getValue(1, 2) * A.getValue(2, 1);
        }
        else {
            for (int i = 1; i <= A.getCols(); i++) {
                determinant += A.getValue(1, i) * Math.pow(-1, i - 1) * determinantRec(createSubMatrix(A, 1, i));
            }
        }
        return determinant;
    }

    private static Matrix createSubMatrix(Matrix A, int row, int col) {
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

    public static Matrix identityMatrix(int n) {
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

    private static Matrix constantIdentityMatrix(int n, double c) {
        double[][] identity = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                identity[i][j] = 0.0;
                if (i == j) {
                    identity[i][j] = c;
                }
            }
        }
        return new Matrix(identity);
    }

    //Matrix Multiplication and helper methods
    private static Matrix vectorFromColumn (Matrix A, int col) {
        Matrix vector = new Matrix(A.getRows(), 1);
        for (int i = 1; i <= A.getRows(); i++) {
            vector.setValue(i, 1, A.getValue(i, col));
        }
        return vector;
    }


    public static Matrix multiplyMatrices(Matrix A, Matrix B) {
        if (A.getCols() != B.getRows()) {
            throw new IllegalArgumentException("The second matrix must have a number of rows equal to the number of columns of the first! " + A.getRows() + " != " +  B.getRows() + "!");
        }
        Matrix product = new Matrix(A.getRows(), B.getCols());
        for (int i = 1; i <= B.getCols(); i++) {
            Matrix tempVector = vectorFromColumn(B, i);
            Matrix tempCol = matrixVectorMultiplication(A, tempVector);
            for (int j = 1; j <= product.getRows(); j++) {
                product.setValue(j, i, tempCol.getValue(j, 1));
            }
        }
        return product;
    }

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
            val = 0.0;
        }
        return newMatrix;
    }

    //Gaussian Elimination, LU decomposition, and helper methods
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

    private static void partialPivoting(Matrix A, int row) {//Can optimize further if I use quick sort
        for ( int i = row - 1; i < A.getRows(); i++ ) {
            for (int j = row; j < A.getRows() - i; j++) {
                if (A.compareScaledRows(j, j+1) < 0) {
                    A.swapRows(j, j+1);
                }
            }
        }
    }

    private static double[] findPivot (Matrix A, int row) {
        double pivot = A.getValue(row, 1);
        int i = 2;
        while (pivot == 0 && i <= A.getCols()) {
            pivot = A.getValue(row, i);
            i++;
        }
        return new double[]{pivot, --i};
    }

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

    public static Matrix gaussianElimination(Matrix matrixOrg) {
        Matrix A;
        try { //Needs to be used so the code in not altering the original A A
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
            
            for (int j = i; j < A.getRows(); j++) {//Getting 0's bellow the pivot
                double ratio = -A.getValue(j + 1, pivotCol)/pivot; //S = -R1/R2 at the pivot column
                A.addRows(j + 1, i, ratio); //R1 = R1 + SR2
            }
        }
        return A;
    }

    private static Matrix backSolve(Matrix matrixOrg) {
        Matrix A;
        try {
            A = (Matrix) matrixOrg.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        for (int i = A.getRows(); i >= 1; i--) {// Gets zero's above the current pivot
            int pivotCol = findPivotCol(A, i);
            double pivot = A.getValue(i, pivotCol);

            if (Math.abs(pivot) > tol) {
                A.scaleRow(i, 1 / pivot); //Scales the current row so the pivot = 1.0.
            }
            for (int j = i; j > 1; j--) {
                A.addRows(j - 1, i, -A.getValue(j - 1, pivotCol));
            }
        }
        return A;
    }

    private static Matrix forwardSolve (Matrix matrixOrg) {
        Matrix A;
        try {
            A = (Matrix) matrixOrg.clone();
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }

        for (int i = 1; i < A.getRows(); i++) {// Gets zero's bellow the current pivot
            int pivotCol = findPivotCol(A, i);
            double pivot = A.getValue(i, pivotCol);
            if (Math.abs(pivot) > tol) {
                A.scaleRow(i, 1 / pivot); //Scales the current row so the pivot = 1.0.
            }
            for (int j = i; j < A.getRows(); j++) {
                A.addRows(j + 1, i, -A.getValue(j+1, pivotCol));
            }
        }
        return A;
    }

    public static Matrix RREF(Matrix matrixOrg) { //Returns RREF(A)
        Matrix A = gaussianElimination(matrixOrg); //Turns A into row echelon form
        //partialPivoting(A, 1); //Pivots to ensure we can backsolve
        A = backSolve(A); // backsolves for the solution.
        return A;
    }

    public static Matrix RREFSolve(Matrix matrixOrg) { //Returns X
        Matrix A = gaussianElimination(matrixOrg); //Turns a matrix into row echelon form
        //partialPivoting(A, 1); //Pivots to ensure we can backsolve
        A = backSolve(A); // backsolves for the solution.
        Matrix x = vectorFromColumn(A, A.getCols());
        return x;
    }

    public static Matrix matrixInverse (Matrix matrixOrg) { //Row Reduce {A, I}
        if (!matrixOrg.isSquare()) {
            throw new IllegalArgumentException("The A must be square!");
        }

        if (determinant(matrixOrg) == 0) { //Is A ~ I? via Invertible Matrix Theorem
            throw new IllegalArgumentException("The A is singular!");
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

        for (int i = 1; i <= cols; i++) {//Making [I A^-1] into [A^-1]
            A.removeCol(1);
        }

        return A;
    }

    //PA = LU, A = P^-1 * LU
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
        U.setAugmentation(false);

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
    }

    public static Matrix LUSolve(Matrix matrixOrg, Matrix b) { //Returns to solution to the system
        matrixOrg.setAugmentation(false);
        Matrix[] LU = LUDecomp(matrixOrg);
        Matrix L = LU[0];
        Matrix U = LU[1];
        Matrix P = LU[2];
        P = matrixInverse(P);
        L = multiplyMatrices(P, L); //L = P^-1 L
        Matrix Ly = RREF(augmentMatrix(L,b)); //Ly = b; Needed because P^-1 L is not triangular
        Matrix y = vectorFromColumn(Ly, Ly.getCols()); //Separates the solution from the matrix

        Matrix Ux = backSolve(augmentMatrix(U, y)); //Ux=y
        Ux.setAugmentation(false);
        partialPivoting(Ux, 1); //Gets the vector in order
        y = vectorFromColumn(Ux, Ux.getCols()); //Separates the solution from the matrix
        return  y;
    }

    //Eigenvalues

    public static Matrix transpose(Matrix A) { //r1 = c1, r2 = c2
        Matrix matrix = new Matrix (A.getCols(), A.getRows());
            for (int i = 1; i <= A.getCols(); i++) {
                for (int j = 1; j <= A.getRows(); j++) {
                    matrix.setValue(i, j, A.getValue(j, i));
                }
            }

        return matrix;
    }

    private static double l2Vector (Matrix b) {
        double l = 0.0;
        for (int i = 1; i <= b.getRows(); i++) {
            l += Math.pow(b.getValue(i, 1), 2);
        }
        return Math.sqrt(l);
    }

    //Todo
    private static double l2Matrix (Matrix A) {
        return 0.0;
    }

    public static double l2Norm (Matrix A) {
        if (A.getCols() == 1) {
            return l2Vector(A);
        }
        else {
            return l2Matrix(A);
        }
    }

    public static double l1Norm(Matrix A) {
        double a = 0.0;
        double b;
        for (int i = 1; i <= A.getCols(); i++) {
            b = 0.0;
            for (int j = 1; j <= A.getRows(); j++) {
                b += Math.abs(A.getValue(j, i));
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

    public static Matrix GramSchmidt (Matrix matrixOrg) {
        Matrix orthMatrix;
        Matrix vn;
        Matrix ci;
        Matrix xn;
        Matrix temp;
        try { //Needs to be used so the code in not altering the original matrix A
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

    public static Matrix[] QRFactorization (Matrix A) {
        Matrix[] QR = new Matrix[2];
        Matrix Q = normalize(GramSchmidt(A));
        Matrix R = multiplyMatrices(transpose(Q), A);
        QR[0] = Q;
        QR[1] = R;
        return QR;
    }

    public static Matrix QREig (Matrix A, double t) {
        if (A.getRows() != A.getCols()) {
            throw new IllegalArgumentException("The matrix must be square!");
        }
        Matrix[] QR = QRFactorization(A);
        Matrix Q = QR[0];
        Matrix R = QR[1];
        Matrix eig = multiplyMatrices(R, Q); //E1 = RQ
        while (nonDiagSum(eig) > t) { //Runs while the sum of the non-diagonal entries is greater than zero
            QR = QRFactorization(eig);
            Q = QR[0];
            R = QR[1];
            eig = multiplyMatrices(R, Q); //E1(k) = R(k)Q(k)
        }
        return eig;
    }

    public static Matrix QREigShift (Matrix A, int  i,int j,double t) { //Shifting is based on element at (i, j)
        if (A.getRows() != A.getCols()) {
            throw new IllegalArgumentException("The matrix must be square!");
        }
        double shift = A.getValue(i, j);
        A = addMatrices(A, constantIdentityMatrix(A.getRows(), shift), -1); //~A(k) = A(k) - shift*I
        Matrix[] QR = QRFactorization(A);
        Matrix Q = QR[0];
        Matrix R = QR[1];
        Matrix eig = addMatrices(multiplyMatrices(R, Q), constantIdentityMatrix(A.getRows(), shift), 1); //E1 = RQ + shift*I
        while (nonDiagSum(eig) > t) { //Runs while the sum of the non-diagonal entries is greater than zero
            shift = eig.getValue(i, j);
            eig = addMatrices(eig, constantIdentityMatrix(A.getRows(), shift), -1); //~E1(k) = E1(k) - shift*I
            QR = QRFactorization(eig);
            Q = QR[0];
            R = QR[1];
            eig = addMatrices(multiplyMatrices(R, Q), constantIdentityMatrix(A.getRows(), shift), 1);  //E1(k) = R(k)Q(k) + shift*I
        }
        return eig;
    }

    public static double[] eig(Matrix A) {
        double[] eVals = new double[A.getCols()];
        Matrix eigs =  QREig(A, tol);
        for (int i = 1; i <= eVals.length; i++) {
            eVals[i-1] = eigs.getValue(i, i);
        }
        return eVals;
    }
}