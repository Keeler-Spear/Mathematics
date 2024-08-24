public class Matrix {
    private double matrix[][];
    private int rows;
    private int cols;
    private boolean square = false;
    private static final int MIN_ROWS = 2;
    private static final int MIN_COLS = 1;


    public Matrix(double[][] vals) {
        if (vals.length < MIN_ROWS || vals[0].length < MIN_COLS) {
            throw new IllegalArgumentException("Your matrix must have at least " + MIN_ROWS + " row(s) and " + MIN_COLS + " column(s)!");
        }
        this.matrix = vals;
        this.rows = vals.length;
        this.cols = vals[0].length;
        if (rows == cols){
            square = true;
        }
    }

    public Matrix(int rows, int cols) {
        if (rows < MIN_ROWS || cols < MIN_COLS) {
            throw new IllegalArgumentException("Your matrix must have at least " + MIN_ROWS + " row(s) and " + MIN_COLS + " column(s)!");
        }
        matrix = new double[rows][cols];
        this.rows = rows;
        this.cols = cols;
        if (rows == cols){
            square = true;
        }
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
            throw new ArrayIndexOutOfBoundsException("Row: " + row + ", Col: " + col + " is out of bounds!");
        }
        return matrix[row - 1][col - 1];
    }

    public int getRows() {
        return rows;
    }

    public int getCols() {
        return cols;
    }

    public boolean isSquare() {
        return square;
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

    public void addRows(int R1, int R2, double scalar) { // R1 = R1 + R2 * Scalar
        for (int i = 0; i < cols; i++) {
            matrix[R1 - 1][i] += matrix[R2 - 1][i] * scalar;
        }
    }

    public void swapRows(int R1, int R2) {
        double[] temRow = matrix[R1 - 1];
        matrix[R1 - 1] = matrix[R2 - 1];
        matrix[R2 - 1] = temRow;

    }

    public void scaleRow(int R, double scalar) {
        for (int i = 0; i < cols; i++) {
            matrix[R - 1][i] = matrix[R - 1][i] * scalar;
        }
    }

    public double largestColVal(int col) {
        double maxVal = 0;
        for (int i = 0; i < rows; i++) {
            if (maxVal < matrix[i][col]) {
                maxVal = matrix[i][col];
            }
        }
        return maxVal;
    }

    public int compareRows(int R1, int R2, int col) {
        int result = 0;

        if (Math.abs(matrix[R1 - 1][col - 1]) > Math.abs(matrix[R2 - 1][col - 1]) ){
            result = 1;
        }
        else if (Math.abs(matrix[R1 - 1][col - 1])  < Math.abs(matrix[R2 - 1][col -1]) ){
            result = -1;
        }
        else if (Math.abs(matrix[R1 - 1][col - 1])  == Math.abs(matrix[R2 - 1][col - 1])  && col != cols){ //Base case
            result = compareRows(R1, R2, col + 1);
        }
        return result;
    }

    @Override
    public String toString() {
        String s = "";
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (j == 0) {
                    s += "[";
                }
                s += String.format("%.3f", matrix[i][j]);
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

    @Override
    protected Object clone() throws CloneNotSupportedException {
        Matrix newMatrix = new Matrix(matrix);
        return newMatrix;
    }
}