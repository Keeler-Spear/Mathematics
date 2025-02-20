/**
 * Represents a mathematical matrix, providing support for internal matrix operations.
 * <p>
 *     This class supports matrix operations such as elementary row operations, row comparison, and more.
 * </p>
 * <p>
 *     For mathematical consistency, 1-based indexing is used for accessing and setting elements.
 * </p>
 *
 * @author Keeler Spear
 * @version %I%, %G%
 * @since 1.0
*/
public class Matrix {
    private static final int MIN_ROWS = 1;
    private static final int MIN_COLS = 1;
    private static final double tol = 0.00000000001;

    private double[][] matrix;
    private int rows;
    private int cols;
    private boolean augmented = false; //Needed to perform scaled partial pivoting

    //Computes if the provided value is "zero."
    private static boolean isZero(double val) {
        if (Math.abs(val) < tol) {
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * Creates a matrix object from the data provided.
     *
     * @param vals A 2d array of doubles representing the matrix entries. Each inner array corresponds to a row.
     * @throws IllegalArgumentException If the matrix does not meet the minimum required dimensions
     * ({@value MIN_ROWS} row(s) and {@value MIN_COLS} column(s)).
     */
    public Matrix(double[][] vals) {
        if (vals.length < MIN_ROWS || vals[0].length < MIN_COLS) {
            throw new IllegalArgumentException("This matrix must have at least " + MIN_ROWS + " row(s) and " + MIN_COLS + " column(s)!");
        }

        this.matrix = vals;
        this.rows = vals.length;
        this.cols = vals[0].length;
    }

    /**
     * Creates a matrix object from the data provided in the form of a vector. This method was added after I submitted
     * my application to your institution.
     *
     * @param vals A 1d array of doubles representing the matrix entries.
     * @throws IllegalArgumentException If the matrix does not meet the minimum required dimensions
     * ({@value MIN_COLS} column(s)).
     */
    public Matrix(double[] vals) {
        if (vals.length < MIN_ROWS) {
            throw new IllegalArgumentException("This matrix must have at least " + MIN_COLS + " column(s)!");
        }

        this.matrix = new double[vals.length][1];

        for (int i = 0; i < vals.length; i++) {
            matrix[i][0] = vals[i];
        }

        this.rows = vals.length;
        this.cols = 1;
    }

    /**
     * Creates a blank matrix object of the size provided.
     *
     * @param rows The number of rows the matrix will have.
     * @param cols The number of columns the matrix will have.
     * @throws IllegalArgumentException If the matrix does not meet the minimum required dimensions
     * ({@value MIN_ROWS} row(s) and {@value MIN_COLS} column(s)).
     */
    public Matrix(int rows, int cols) {
        if (rows < MIN_ROWS || cols < MIN_COLS) {
            throw new IllegalArgumentException("This matrix must have at least " + MIN_ROWS + " row(s) and " + MIN_COLS + " column(s)!");
        }

        matrix = new double[rows][cols];
        this.rows = rows;
        this.cols = cols;
    }

    /**
     * Returns the entries of the matrix.
     *
     * @return The 2d array of doubles representing the matrix's entries.
     */
    public double[][] getMatrix() {
        return matrix;
    }

    /**
     * Returns to total number of rows in the matrix.
     *
     * @return The number of rows the matrix has.
     */
    public int getRows() {
        return rows;
    }

    /**
     * Returns to total number of columns in the matrix.
     *
     * @return The number of columns the matrix has.
     */
    public int getCols() {
        return cols;
    }

    /**
     * Returns a boolean representing if the matrix is square.
     *
     * @return True if the matrix is square, false if otherwise.
     */
    public boolean isSquare() {
        return rows == cols;
    }

    /**
     * Returns a boolean representing if the matrix is augmented.
     *
     * @return True if the matrix is augmented, false if otherwise.
     */
    public boolean getAugmentation () {
        return augmented;
    }

    /**
     * Sets the matrix parameter governing if the matrix is augmented.
     *
     * @param augmented The boolean that the internal augmented variable will be set to.
     */
    public void setAugmentation (boolean augmented) {
        this.augmented = augmented;
    }

    /**
     * Sets the value of the matrix at the specified position using 1-based indexing.
     *
     * @param row The 1-based row index where the value will be collected.
     * @param col The 1-based column index where the value will be collected.
     * @return The value at the specified position in the matrix.
     * @throws IllegalArgumentException If the specified row or column index is out of the matrix's bounds
     * (1 ≤ row ≤ {@code rows}, 1 ≤ col ≤ {@code cols}).
     */
    public double getValue(int row, int col) {
        if (row <= 0 || row > rows || col <= 0 || col > cols) {
            throw new IllegalArgumentException("Row: " + row + ", Col: " + col + " is out of bounds!");
        }

        return matrix[row - 1][col - 1];
    }

    /**
     * Sets the value of the matrix at the specified position using 1-based indexing.
     *
     * @param row The 1-based row index where the value will be set.
     * @param col The 1-based column index where the value will be set.
     * @param value The value to set at the specified position in the matrix.
     * @throws IllegalArgumentException If the specified row or column index is out of the matrix's bounds
     * (1 ≤ row ≤ {@code rows}, 1 ≤ col ≤ {@code cols}).
     */
    public void setValue (int row, int col, double value) {
        if (row <= 0 || row > rows || col <= 0 || col > cols) {
            throw new IllegalArgumentException("Row: " + row + ", Col: " + col + " is out of bounds!");
        }

        matrix[row-1][col-1] = value;
    }

    /**
     * Adds a row to the top of the matrix with the values provided. This method was added after I submitted my
     * application to your institution.
     *
     * @param vals The 1d array of values to be added to the matrix as a row.
     * @throws IllegalArgumentException If the input row's length does not match the number of columns in the matrix.
     */
    public void addRowTop(double[] vals) {
        if (vals.length != cols) {
            throw new IllegalArgumentException("This row is not the same size as the others in this matrix!");
        }

        rows++;

        double[][] tempMatrix = new double[rows][cols];

        for (int i = 1; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                tempMatrix[i][j] = matrix[i - 1][j];
            }
        }

        for (int i = 0; i < cols; i++) {
            tempMatrix[0][i] = vals[i];
        }

        matrix = tempMatrix;
    }

    /**
     * Adds a row to the top of the matrix with the values provided. This method was added after I submitted my
     * application to your institution.
     *
     * @param vals The 2d array of values to be added to the matrix as a row.
     * @throws IllegalArgumentException If the input row's length does not match the number of columns in the matrix.
     */
    public void addRowTop(double[][] vals) {
        if (vals.length != cols) {
            throw new IllegalArgumentException("This row is not the same size as the others in this matrix!");
        }

        rows++;

        double[][] tempMatrix = new double[rows][cols];

        for (int i = 1; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                tempMatrix[i][j] = matrix[i - 1][j];
            }
        }

        for (int i = 0; i < cols; i++) {
            tempMatrix[0][i] = vals[i][0];
        }

        matrix = tempMatrix;
    }

    /**
     * Adds a row to the bottom of the matrix with the values provided.
     *
     * @param vals The 1d array of values to be added to the matrix as a row.
     * @throws IllegalArgumentException If the input row's length does not match the number of columns in the matrix.
     */
    public void addRowBottom(double[] vals) {
        if (vals.length != cols) {
            throw new IllegalArgumentException("This row is not the same size as the others in this matrix!");
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

    /**
     * Adds a row to the bottom of the matrix with the values provided. This method was added after I submitted my
     * application to your institution.
     *
     * @param vals The 2d array of values to be added to the matrix as a row.
     * @throws IllegalArgumentException If the input row's length does not match the number of columns in the matrix.
     */
    public void addRowBottom(double[][] vals) {
        if (vals.length != cols) {
            throw new IllegalArgumentException("This row is not the same size as the others in this matrix!");
        }

        rows++;

        double[][] tempMatrix = new double[rows][cols];

        for (int i = 0; i < rows - 1; i++) {
            for (int j = 0; j < cols; j++) {
                tempMatrix[i][j] = matrix[i][j];
            }
        }

        for (int i = 0; i < cols; i++) {
            tempMatrix[rows - 1][i] = vals[i][0];
        }

        matrix = tempMatrix;
    }


    /**
     * Removes the matrix's row at the specified position using 1-based indexing.
     *
     * @param row The 1-based index where the row will be removed.
     * @throws IllegalArgumentException If the specified row index is outside the matrix's bounds
     * (1 ≤ row ≤ {@code rows}).
     */
    public void removeRow(int row) {
        if (rows < row || row <= 0) {
            throw new IllegalArgumentException("This row is out of bounds!");
        }

        row--;
        rows--;

        double[][] tempMatrix = new double[rows][cols];

        for (int i = 0; i < row; i++) {
            tempMatrix[i] = matrix[i];
        }

        for (int i = row; i < rows; i++) {
            tempMatrix[i] = matrix[i+1];
        }

        matrix = tempMatrix;
    }

    /**
     * Adds a column to the left side of the matrix with the values provided. This method was added after I submitted my
     * application to your institution.
     *
     * @param vals The 1d array of values to be added to the matrix as a column.
     * @throws IllegalArgumentException If the input column's length does not match the number of rows in the matrix.
     */
    public void addColLeft(double[] vals) {
        if (vals.length != rows) {
            throw new IllegalArgumentException("This column is not the same size as the others in this matrix!");
        }

        cols++;
        double[][] tempMatrix = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 1; j < cols; j++) {
                tempMatrix[i][j] = matrix[i][j- 1];
            }
        }

        for (int i = 0; i < rows; i++) {
            tempMatrix[i][0] = vals[i];
        }

        matrix = tempMatrix;
    }

    /**
     * Adds a column to the left side of the matrix with the values provided. This method was added after I submitted my
     * application to your institution.
     *
     * @param vals The 2d array of values to be added to the matrix as a column.
     * @throws IllegalArgumentException If the input column's length does not match the number of rows in the matrix.
     */
    public void addColLeft(double[][] vals) {
        if (vals.length != rows) {
            throw new IllegalArgumentException("This column is not the same size as the others in this matrix!");
        }

        cols++;
        double[][] tempMatrix = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 1; j < cols; j++) {
                tempMatrix[i][j] = matrix[i][j- 1];
            }
        }

        for (int i = 0; i < rows; i++) {
            tempMatrix[i][0] = vals[i][0];
        }

        matrix = tempMatrix;
    }

    /**
     * Adds a column to the right side of the matrix with the values provided.
     *
     * @param vals The 1d array of values to be added to the matrix as a column.
     * @throws IllegalArgumentException If the input column's length does not match the number of rows in the matrix.
     */
    public void addColRight(double[] vals) {
        if (vals.length != rows) {
            throw new IllegalArgumentException("This column is not the same size as the others in this matrix!");
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

    /**
     * Adds a column to the right side of the matrix with the values provided. This method was added after I submitted
     * my application to your institution.
     *
     * @param vals The 2d array of values to be added to the matrix as a column.
     * @throws IllegalArgumentException If the input column's length does not match the number of rows in the matrix.
     */
    public void addColRight(double[][] vals) {
        if (vals.length != rows) {
            throw new IllegalArgumentException("This column is not the same size as the others in this matrix!");
        }

        cols++;
        double[][] tempMatrix = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols - 1; j++) {
                tempMatrix[i][j] = matrix[i][j];
            }
        }

        for (int i = 0; i < rows; i++) {
            tempMatrix[i][cols - 1] = vals[i][0];
        }

        matrix = tempMatrix;
    }

    /**
     * Removes the matrix's column at the specified position using 1-based indexing.
     *
     * @param col The 1-based index where the column will be removed.
     * @throws IllegalArgumentException If the specified column index is outside the matrix's bounds
     * (1 ≤ col ≤ {@code cols}).
     */
    public void removeCol(int col) {
        if (cols < col || col <= 0) {
            throw new IllegalArgumentException("This column is out of bounds");
        }

        col--;
        cols--;

        double[][] tempMatrix = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < col; j++) {
                tempMatrix[i][j] = matrix[i][j];
            }
        }

        for (int i = 0; i < rows; i++) {
            for (int j = col; j < cols; j++) {
                tempMatrix[i][j] = matrix[i][j+1];
            }
        }

        matrix = tempMatrix;
    }

    /**
     * Sets the column in the matrix at the specified position to the one provided using 1-based indexing.
     *
     * @param col The 1-based column index where the values will be set.
     * @param vals A 1d array of values that will be inserted into the matrix at the given location.
     * @throws IllegalArgumentException If the specified column index is out of the matrix's bounds
     * (1 ≤ col ≤ {@code cols}).
     * @throws IllegalArgumentException If the input column's length does not match the number of rows in the matrix.
     */
    public void setCol(int col, double[] vals) {
        col--;

        if (col < 0 || col >= cols) {
            throw new IllegalArgumentException("The column out of bounds!");
        }

        if (vals.length != rows) {
            throw new IllegalArgumentException("The given column is not the correct size!");
        }

        for (int i = 0; i < rows; i++) {
            matrix[i][col] = vals[i];
        }
    }

    /**
     * Sets the column in the matrix at the specified position to the one provided using 1-based indexing.
     *
     * @param col The 1-based column index where the values will be set.
     * @param vals A 2d array of values that will be inserted into the matrix at the given location.
     * @throws IllegalArgumentException If the specified column index is out of the matrix's bounds
     * (1 ≤ col ≤ {@code cols}).
     * @throws IllegalArgumentException If the input column's length does not match the number of rows in the matrix.
     */
    public void setCol(int col, double[][] vals) {
        col--;

        if (col < 0 || col >= cols) {
            throw new IllegalArgumentException("The column out of bounds!");
        }

        if (vals.length != rows || vals[0].length != 1) {
            throw new IllegalArgumentException("The given column is not the correct size!");
        }

        for (int i = 0; i < rows; i++) {
            matrix[i][col] = vals[i][0];
        }
    }

    /**
     * Sets the row in the matrix at the specified position to the one provided using 1-based indexing.
     *
     * @param row The 1-based row index where the values will be set.
     * @param vals A 1d array of values that will be inserted into the matrix at the given location.
     * @throws IllegalArgumentException If the specified row index is out of the matrix's bounds
     * (1 ≤ row ≤ {@code rows}).
     * @throws IllegalArgumentException If the input column's length does not match the number of rows in the matrix.
     */
    public void setRow(int row, double[] vals) {
        row--;

        if (row < 0 || row >= rows) {
            throw new IllegalArgumentException("The row out of bounds!");
        }

        if (vals.length != cols) {
            throw new IllegalArgumentException("The given row is not the correct size!");
        }

        for (int i = 0; i < cols; i++) {
            matrix[row][i] = vals[i];
        }
    }

    /**
     * Sets the row in the matrix at the specified position to the one provided using 1-based indexing.
     *
     * @param row The 1-based row index where the values will be set.
     * @param vals A 2d array of values that will be inserted into the matrix at the given location.
     * @throws IllegalArgumentException If the specified row index is out of the matrix's bounds
     * (1 ≤ row ≤ {@code rows}).
     * @throws IllegalArgumentException If the input column's length does not match the number of rows in the matrix.
     */
    public void setRow(int row, double[][] vals) {
        row--;

        if (row < 0 || row >= rows) {
            throw new IllegalArgumentException("The row out of bounds!");
        }

        if (vals[0].length != cols || vals.length != 1) {
            throw new IllegalArgumentException("The given row is not the correct size!");
        }

        for (int i = 0; i < cols; i++) {
            matrix[row][i] = vals[0][i];
        }
    }

    /**
     * The elementary matrix operation of adding one row to another within the matrix. Row R1 = row R1 + scalar * row R2.
     *
     * @param R1 The 1-based index where the row will be replaced by the sum of it and another row.
     * @param R2 The 1-based index of the row that will be multiplied by a scalar and added to the other.
     * @param scalar The value the second row (R2) will be multiplied by in the row addition.
     * @throws IllegalArgumentException * @throws IllegalArgumentException If the provided row indices are out of bounds
     * (1 ≤ R1, R2 ≤ {@code rows}).
     */
    public void addRows(int R1, int R2, double scalar) {
        if (R1 <= 0 || R1 > rows || R2 <= 0 || R2 > rows) {
            throw new IllegalArgumentException("One or more of the rows are out of bounds!");
        }

        for (int i = 0; i < cols; i++) {
            matrix[R1 - 1][i] += matrix[R2 - 1][i] * scalar;
        }
    }

    /**
     * The elementary matrix operation of swapping two rows within the matrix. Row R1 = row R2, row R2 = row R1.
     *
     * @param R1 The 1-based index of the first row to be swapped.
     * @param R2 The 1-based index of the second row to be swapped.
     * @throws IllegalArgumentException * @throws IllegalArgumentException If the provided row indices are out of bounds
     * (1 ≤ R1, R2 ≤ {@code rows}).
     */
    public void swapRows(int R1, int R2) {
        if (R1 <= 0 || R1 > rows || R2 <= 0 || R2 > rows) {
            throw new IllegalArgumentException("One or more of the rows are out of bounds!");
        }

        double[] temRow = matrix[R1 - 1];
        matrix[R1 - 1] = matrix[R2 - 1];
        matrix[R2 - 1] = temRow;
    }

    /**
     * The elementary matrix operation of scaling a row within the matrix by a scalar. Row R = Row R * scalar.
     *
     * @param R The 1-based index of the row to be scaled.
     * @param scalar The value the row will be multiplied by.
     * @throws IllegalArgumentException If the specified row index is outside the matrix's bounds (1 ≤ R ≤ {@code rows}).
     */
    public void scaleRow(int R, double scalar) {
        if (R <= 0 || R > rows) {
            throw new IllegalArgumentException("This row is out of bounds!");
        }

        for (int i = 0; i < cols; i++) {
            matrix[R - 1][i] = matrix[R - 1][i] * scalar;
        }
    }

    /**
     * Compares two rows of the matrix. Determines which row has the largest leading value.
     *
     * @param R1 The 1-based index of the first row to compare.
     * @param R2 The 1-based index of the second row to compare.
     * @param col The column in which row values will be compared.
     * @return An integer indicating which row is larger:
     *         <ul>
     *             <li> 1 if the first row is larger</li>
     *             <li> -1 if the second row is larger</li>
     *             <li> 0 if the rows are equal</li>
     *         </ul>
     * @throws IllegalArgumentException * @throws IllegalArgumentException If the provided row indices are out of bounds
     * (1 ≤ R1, R2 ≤ {@code rows}).
     */
    public int compareRows(int R1, int R2, int col) {
        int result = 0;
        double val1 = Math.abs(matrix[R1 - 1][col - 1]);
        double val2 = Math.abs(matrix[R2 - 1][col - 1]);

        if (isZero(val1 - val2) && (col != cols)) {
            result = compareRows(R1, R2, col + 1);
        }
        else if (!isZero(val1 - val2)) {
            result = 1;
        }
        else if (val1 - val2 < -1.0 * tol) {
            result = -1;
        }

        return result;
    }

    /**
     * Compares two rows of the matrix after each is scaled by dividing all elements by the largest absolute value
     * in the row. Determines which row has the largest leading value.
     *
     * @param R1 The 1-based index of the first row to compare.
     * @param R2 The 1-based index of the second row to compare.
     * @return An integer indicating which row is larger:
     *         <ul>
     *             <li> 1 if the first row is larger after scaling</li>
     *             <li> -1 if the second row is larger after scaling</li>
     *             <li> 0 if the rows are equal after scaling</li>
     *         </ul>
     * @throws IllegalArgumentException If the provided row indices are out of bounds (1 ≤ R1, R2 ≤ {@code rows}).
     */
    public int compareScaledRows(int R1, int R2) {
        if (R1 <= 0 || R1 > rows || R2 <= 0 || R2 > rows) {
            throw new IllegalArgumentException("One or both row indices are invalid!");
        }

        double[] r1 = matrix[R1 - 1].clone();
        double[] r2 = matrix[R2 - 1].clone();

        r1 = scaleRowToOne(r1);
        r2 = scaleRowToOne(r2);

        return compareScaledRowsRec(r1, r2, 1);
    }

    //Scales an input row to one; used in scaled row comparison.
    private double[] scaleRowToOne(double [] row) {
        double max = 0.0;
        double augFactor = 0.0;

        if (augmented) {
            augFactor = 1.0;
        }

        for (int i = 0; i < row.length - augFactor; i++) {
            if (Math.abs(row[i]) > Math.abs(max)) {
                max = row[i];
            }
        }

        for (int i = 0; i < row.length; i++) {
            row[i] = row[i] / max;
        }

        return row;
    }

    //Recursively compares two scaled arrays.
    private int compareScaledRowsRec(double[] r1, double[] r2, int col) {
        int result = 0;
        double val1 = Math.abs(r1[col - 1]);
        double val2 = Math.abs(r2[col - 1]);

        if (isZero(val1  - val2) && col != cols){
            result = compareScaledRowsRec(r1, r2, col + 1);
        }
        else if (val1 - val2 > tol){
            result = 1;
        }
        else if (val1 - val2 < -1.0 * tol){
            result = -1;
        }
        return result;
    }

    /**
     * Returns the matrix as a String in NumPy formating. This method was added after I submitted my application to your
     * institution.
     *
     * @return The matrix as a multi-line String, with each print line containing a row. The matrix's entries are printed
     * to four decimal places.
     */
    public String npString() {
        String s = "";
        if (rows > 1) {
            s += "[";
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (j == 0) {
                    s += "[";
                }
                s += String.format("%.4f", matrix[i][j]);
                if (j == cols - 1 && i != rows - 1) {
                    s += "], ";
                }
                else if (j == cols - 1 && i == rows - 1) {
                    s += "]";
                }
                else {
                    s += ", ";
                }
            }
        }

        if (rows > 1) {
            s += "]";
        }

        return s;
    }

    /**
     * Returns the matrix as a String.
     *
     * @return The matrix as a multi-line String, with each print line containing a row. The matrix's entries are printed
     * to four decimal places.
     */
    @Override
    public String toString() {
        String s = "";
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (j == 0) {
                    s += "[";
                }
                s += String.format("%.4f", matrix[i][j]);
                if (j == cols - 1 && i != rows - 1) {
                    s += "]\n";
                }
                else if (j == cols - 1 && i == rows - 1) {
                    s += "]";
                }
                else {
                    s += " ";
                }
            }
        }

        return s;
    }

    /**
     * Clones the matrix as an object.
     *
     * @return A clone of the matrix as an object.
     */
    @Override
    protected Object clone() throws CloneNotSupportedException { //Copies the array and then makes a new matrix with the copy
        double[][] vals = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                vals[i][j] = matrix[i][j];
            }
        }

        Matrix newMatrix = new Matrix(vals.clone());
        newMatrix.setAugmentation(augmented);

        return newMatrix;
    }

    /**
     * Evaluates if this matrix and another object are the same based on class; row and column count; and matrix entries.
     *
     * @param obj The object the matrix will be compared to.
     * @return A boolean representing if the matrix and the provided object are the same.
     */
    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }

        if (obj.getClass() != this.getClass()) {
            return false;
        }

        Matrix matrix2 = (Matrix) obj;
        if (rows != matrix2.rows || cols != matrix2.cols) {
            return false;
        }

        for (int i = 1; i <= rows; i++) {
            for (int j = 1; j <= cols; j++) {
                if (!isZero((getValue(i, j)) - (matrix2.getValue(i, j)))) {
                    return false;
                }
            }
        }

        return true;
    }
}