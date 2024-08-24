import java.util.Scanner;

public class eigenApproximator{
    //Print Matrix Method
    private static void printMatrix(double [][] matrix) {
        int i;
        int j;
        int size;
        size = matrix.length;

        for (i = 0; i < size; i++) {
            for (j = 0; j < size; j++) {
                if (j == 0) {
                    System.out.print("[");
                    System.out.printf("%.2f", matrix[i][j]);
                    System.out.print(" ");
                }
                else if (j == (size -1)) {
                    System.out.printf("%.2f", matrix[i][j]);
                    System.out.print("]");
                }
                else {
                    System.out.printf("%.2f", matrix[i][j]);
                    System.out.print(" ");
                }
            }
            System.out.println();
        }
    }

    private static void printVector(double [][] vector) {
        int i;
        int size;
        size = vector.length;

        for (i = 0; i < size; i++) {
            System.out.print("[");
            System.out.printf("%.4f", vector[i][0]);
            System.out.println("]");
        }
    }

    //Multiple Matrixes Together
    private  static double[][] multiplyMatrix(double[][] matrix, double[][]vector) {
        double[][] newVector = new double [matrix.length][1];
        int i;
        int j;
        double value = 0.0;
        for (i = 0; i < matrix.length; ++i) {
            for (j = 0; j < matrix.length; ++j) {
                value = value + (matrix[i][j] * vector[j][0]);
            }
            newVector[i][0] = value;
            value = 0.0;
        }

        return newVector;
    }

    //Multiply Matrix by a Scalar
    private static double[][] multiplyScalar(double[][] matrix, double scalar) {
        int i;
        int j;

        for (i = 0; i < matrix.length; ++i) {
            for (j = 0; j < matrix[0].length; ++j) {
                matrix[i][j] = scalar * matrix[i][j];
            }
        }
        return matrix;
    }

    //Find Greatest Vector Value
    private static double greatestValue(double[][] vector) {
        double greatestValue = 0.0;
        int i;
        for (i = 0; i < vector.length; ++i) {
            if (Math.abs(vector[i][0]) > Math.abs(greatestValue)) {
                greatestValue = vector[i][0];
            }
        }

        return greatestValue;
    }

    //Approximator
    private static void approximator(double[][] matrix, double[][] vector, int k) {
        double mew = 0.0;
        int i;
        double[][] tempVector = new double [matrix.length][1];

        for (i = 0; i < k; ++i) {
            tempVector = multiplyMatrix(matrix, vector);
            mew = greatestValue(tempVector);
            vector = multiplyScalar(tempVector, (1.0 / mew));
        }
        System.out.println("Your Dominant Eigenvalue is:");
        System.out.printf("%.4f", mew);
        System.out.println();
        System.out.println();
        System.out.println("Your Dominant Eigenvector is:");
        printVector(vector);
    }


    public static void main (String[] args) {
        Scanner scnr = new Scanner(System.in);
        int size = 3;
        int iterations = 500;
        int i;
        int j;
        int tempRow;
        int tempCol;
        double newValue;
        int choice1 = 0;
        double maxVal;

        //Gets the size of the user's matrix and makes an nxn matrix.
        System.out.println("What is the size of your matrix?");
        size = scnr.nextInt();
        System.out.println();
        double[][] matrix = new double [size][size];
        double[][] initialVector = new double [size][1];

        //Recives user doubles and puts them into the matrix.
        System.out.print("Please enter the values in your matrix.");
        for (i = 0; i < size; ++i) {
            for (j = 0; j < size; ++j ) {
                matrix[i][j] = scnr.nextDouble();
            }
            if (i < size - 1) {
                System.out.println("You are now on row " + (i+2) + ".");
            }

        }
        System.out.println();
        while (choice1 == 0) {
            System.out.println("Is this the matrix you want to use?");
            printMatrix(matrix);
            System.out.println();
            System.out.println("Type 1 if the matrix is correct. Type 0 if it is incorrect and you would like to replace an element in the matrix.");
            choice1 = scnr.nextInt();
            //Allows the user to change a matrix element
            if (choice1 == 0) {
                System.out.println("What row is the element that you want to change in?");
                tempRow = scnr.nextInt();
                System.out.println("What column is the element that you want to change in?");
                tempCol = scnr.nextInt();
                System.out.println("Please enter the new value for " + tempRow +", " + tempCol +".");
                newValue = scnr.nextDouble();
                matrix[tempRow-1][tempCol-1] = newValue;
            }
        }
        choice1 = 0;

        System.out.println();
        System.out.println("Please enter the values in your initial vector.");
        for (i = 0; i < size; ++i) {
            initialVector[i][0] = scnr.nextDouble();
        }
        System.out.println();
        while (choice1 == 0) {
            System.out.println("Is this the vector you want to use?");
            printVector(initialVector);
            System.out.println();
            System.out.println("Type 1 if the vector is correct. Type 0 if it is incorrect and you would like to replace an element in the matrix.");
            choice1 = scnr.nextInt();
            //Allows the user to change a matrix element
            if (choice1 == 0) {
                System.out.println("What row is the element that you want to change in?");
                tempRow = scnr.nextInt();
                System.out.println("Please enter the new value for " + tempRow +", " + "1" +".");
                newValue = scnr.nextDouble();
                initialVector[tempRow-1][0] = newValue;
            }
        }
        System.out.println();
        System.out.println("How many iterations do you want to do?");
        iterations = scnr.nextInt();
        //Makes the initial vector have the greatest value of 1.
        maxVal = (1.0 / greatestValue(initialVector));
        initialVector = multiplyScalar(initialVector, maxVal);

        System.out.println();
        approximator(matrix, initialVector, iterations);



    }
}