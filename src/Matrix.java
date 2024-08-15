public class Matrix {
    private double matrix[][];
    private int rows;
    private int cols;
    private boolean square = false;
    private static final int MIN_SIZE = 2;

    public Matrix(double[][] vals) {
        if (vals.length < MIN_SIZE || vals[0].length < MIN_SIZE) {
            throw new IllegalArgumentException("Your matrix must have at least " + MIN_SIZE + " rows and " + MIN_SIZE + " columns");

        }
        this.matrix = vals;
        this.rows = vals.length;
        this.cols = vals[0].length;
        if (rows == cols){
            square = true;
        }
    }

    public Matrix(int rows, int cols) {
        if (rows < MIN_SIZE || cols < MIN_SIZE) {
            throw new IllegalArgumentException("Your matrix must have at least " + MIN_SIZE + " rows and " + MIN_SIZE + " columns");

        }
        matrix = new double[rows][cols];
        this.rows = rows;
        this.cols = cols;
        if (rows == cols){
            square = true;
        }
    }

    private Matrix createSubMatrix(Matrix baseMatrix, int col) {
        double[][] tempMatrix1 = new double[baseMatrix.getRows() - 1][baseMatrix.getCols()];
        tempMatrix1[0] = baseMatrix.getMatrix()[1];
        for (int i = 1; i < baseMatrix.getRows() - 1; i++) {
            tempMatrix1[i] = baseMatrix.getMatrix()[i + 1];
        }
        double[][] tempMatrix2 = new double[baseMatrix.getRows() - 1][baseMatrix.getCols() - 1];
        for (int i = 0; i < baseMatrix.getRows() - 1; i++) {
            for (int j = 0; j < col; j++) {
                tempMatrix2[i][j] = tempMatrix1[i][j];
            }
        }
        for (int i = 0; i < baseMatrix.getRows() - 1; i++) {
            for (int j = col; j < baseMatrix.getCols() - 1; j++) {
                tempMatrix2[i][j] = tempMatrix1[i][j+1];
            }
        }
        return new Matrix(tempMatrix2);
    }

    public void setValue (int row, int col, double value) {
        if (row < 0 || row > rows || col < 0 || col > cols) {
            throw new ArrayIndexOutOfBoundsException();
        }
        matrix[row-1][col-1] = value;
    }

    public double[][] getMatrix() {
        return matrix;
    }

    public double getValue(int row, int col) {
        if (row < 0 || row > rows || col < 0 || col > cols) {
            throw new ArrayIndexOutOfBoundsException();
        }
        return matrix[row - 1][col - 1];
    }

    public int getRows() {
        return rows;
    }

    public int getCols() {
        return cols;
    }

//    public void addRow(double[] vals) {
//        if (vals.length != cols) {
//            throw new ArrayIndexOutOfBoundsException();
//        }
//        rows++;
//        double[][] tempMatrix = new double[rows][cols];
//        for (int i = 0; i < rows - 1; i++) {
//            for (int j = 0; j < cols; j++) {
//                tempMatrix[i][j] = matrix[i][j];
//            }
//        }
//        for (int i = 0; i < cols; i++) {
//            tempMatrix[rows - 1][i] = vals[i];
//        }
//        matrix = tempMatrix;
//    }
//
//    public void addCol(double[] vals) {
//        if (vals.length != rows) {
//            throw new ArrayIndexOutOfBoundsException();
//        }
//        cols++;
//        double[][] tempMatrix = new double[rows][cols];
//        for (int i = 0; i < rows; i++) {
//            for (int j = 0; j < cols - 1; j++) {
//                tempMatrix[i][j] = matrix[i][j];
//            }
//        }
//        for (int i = 0; i < rows; i++) {
//            tempMatrix[i][cols - 1] = vals[i];
//        }
//        matrix = tempMatrix;
//    }

    public void removeRow(int index) {
        if (rows < index || index <= 0) {
            throw new ArrayIndexOutOfBoundsException();
        }
        index = index - 1;
        rows--;
        double[][] tempMatrix = new double[rows][cols];
        for (int i = 0; i < index; i++) {
            tempMatrix[i] = matrix[i];
        }
        for (int i = index; i < rows; i++) {
            tempMatrix[i] = matrix[i+1];
        }
        matrix = tempMatrix;
    }

    public void removeCol(int index) {
        if (cols < index || index <= 0) {
            throw new ArrayIndexOutOfBoundsException();
        }
        index = index - 1;
        cols--;
        double[][] tempMatrix = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < index; j++) {
                tempMatrix[i][j] = matrix[i][j];
            }
        }
        for (int i = 0; i < rows; i++) {
            for (int j = index; j < cols; j++) {
                tempMatrix[i][j] = matrix[i][j+1];
            }
        }
        matrix = tempMatrix;
    }

    private void addRows(int R1, int R2, double scalar) { // R1 = R1 + R2 * Scalar
        for (int i = 0; i < cols; i++) {
            matrix[R1][i] += matrix[R2][i] * scalar;
        }
    }

    private void swapRows(int R1, int R2) {
        double[] temRow = matrix[R1];
        matrix[R1] = matrix[R2];
        matrix[R2] = temRow;

    }

    private void scaleRow(int R, double scalar) {
        for (int i = 0; i < cols; i++) {
            matrix[R][i] = matrix[R][i] * scalar;
        }
    }

    private double largestColVal(int col) {
        double maxVal = 0;
        for (int i = 0; i < rows; i++) {
            if (maxVal < matrix[i][col]) {
                maxVal = matrix[i][col];
            }
        }
        return maxVal;
    }

    private double[][] vectorFromColumn (Matrix tempMatrix, int col) {
        double[][] vector = new double[tempMatrix.getRows()][1];
        for (int i = 0; i < tempMatrix.getRows(); i++) {
            vector[i][0] = tempMatrix.getValue(i, col - 1);
        }
        return vector;
    }

    //Todo Finish matrix-vector multiplication method
//    private Matrix matrixVectorMultiplication (Matrix matrix1, Matrix matrix2) {
//
//    }

    //ToDo Finish multiplication method
    //Calle * Peram Matrices
//    public Matrix multiplyMatrices(Matrix matrix2) {
//        if (cols != matrix2.getRows()) {
//            throw new IllegalArgumentException("The second matrix must have a number of rows equal to the number of columns of the first! " + cols + " != " +  matrix2.getRows() + "!");
//        }
//    }

    public double determinant(Matrix subMatrix) {
        if (!square || (subMatrix.getRows() == subMatrix.getCols() && subMatrix.getRows() < 2)) {
            throw new IllegalArgumentException("Matrix is not square");
        }
        double determinant = 0.0;
        //Base case
        if (subMatrix.getRows() == 2) {
            determinant += subMatrix.getValue(0,0) * subMatrix.getValue(1,1) - subMatrix.getValue(0,1) * subMatrix.getValue(1,0);
            return determinant;
        }
        for (int i = 0; i < subMatrix.getRows(); i++) {
            determinant += subMatrix.getValue(0, i)  * Math.pow(-1, i) * determinant(createSubMatrix(subMatrix, i));
        }
        return determinant;
    }

    @Override
    public String toString() {
        String s = "";
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (j == 0) {
                    s += "[";
                }
                s += matrix[i][j];
                if (j == cols - 1) {
                    s += "]\n";
                }
                else {
                    s += " ";
                }
            }
        }
        return s;
    }
}