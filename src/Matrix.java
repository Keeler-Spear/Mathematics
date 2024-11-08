public class Matrix {
    private double[][] matrix;
    private int rows;
    private int cols;
    private boolean square = false;
    private static final int MIN_ROWS = 1;
    private static final int MIN_COLS = 1;
    private static final double h = 0.00000000001;


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
        return rows == cols;
    }

    public void addRow(double[] vals) {
        if (vals.length != cols) {
            throw new ArrayIndexOutOfBoundsException();
        }
        rows++;
        double[][] tempMatrix = new double[rows][cols];
        for (int i = 0; i < rows - 1; i++) {
            for (int j = 0; j < cols; j++) {
                tempMatrix[i][j] = matrix[i][j];
            }
        }
        for (int i = 0; i < cols; i++) {
            tempMatrix[rows - 1][i] = vals[i];
        }
        matrix = tempMatrix;
    }

    public void addCol(double[] vals) {
        if (vals.length != rows) {
            throw new ArrayIndexOutOfBoundsException();
        }
        cols++;
        double[][] tempMatrix = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols - 1; j++) {
                tempMatrix[i][j] = matrix[i][j];
            }
        }
        for (int i = 0; i < rows; i++) {
            tempMatrix[i][cols - 1] = vals[i];
        }
        matrix = tempMatrix;
    }

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

    public void swapRows(int R1, int R2) { //Maybe use permutation matrix and left multiply it by the current matrix
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
        double val1 = Math.abs(matrix[R1 - 1][col - 1]);
        double val2 = Math.abs(matrix[R2 - 1][col - 1]);

        if (Math.abs(val1  - val2) < h && col != cols){ //Base case
            result = compareRows(R1, R2, col + 1);
        }
        else if (val1 - val2 > h ){ //0.00  and -0.00 are triggered here
            result = 1;
        }
        else if (val1 - val2 < -1.0 * h){
            result = -1;
        }

        return result;
    }

    public void setCol(int col, double[][] v) {
        col = col - 1;
        if (col < 0 || col >= cols) {
            throw new ArrayIndexOutOfBoundsException();
        }
        if (v.length != rows || v[0].length != 1) {
            throw new IllegalArgumentException("The given column is not the correct size!");
        }
        for (int i = 0; i < rows; i++) {
            matrix[i][col] = v[i][0];
        }
    }

    @Override
    public String toString() {
        String s = "";
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (j == 0) {
                    s += "[";
                }
                s += String.format("%.4f", matrix[i][j]);
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
    protected Object clone() throws CloneNotSupportedException { //Copies the array and then makes a new matrix with the copy
        double[][] vals = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                vals[i][j] = matrix[i][j];
            }
        }
        Matrix newMatrix = new Matrix(vals);

        return newMatrix;
    }
    public boolean equals(Matrix matrix2) {
        boolean equal = true;
        if (rows == matrix2.rows && cols == matrix2.cols) {
            for (int i = 1; i <= rows; i++) {
                for (int j = 1; j <= cols; j++) {
                    if (getValue(i, j) != matrix2.getValue(i, j)) {
                        equal = false;
                    }
                }
            }
        }
        else {
            equal = false;
        }
        return equal;
    }
}