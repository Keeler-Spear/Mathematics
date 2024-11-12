import javax.swing.*;

public class UI {
    private static int curLetter = 0;
    private static final String[] names = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"};
    private static Matrix[] matrices = new Matrix[26];
    private static String[] options = {"Create Matrix", "Perform Matrix Operations", "View Saved Matrices", "TO BE ADDED", "Exit"};
    private static String[] opOptions = {"Determinant", "Eigenvalues", "LU Decomposition", "Matrix Inverse", "Matrix Multiplication", "QR Decomposition", "RREF"};

    private static final String SAVE_FILE = "SavedVars.csv";

    private static Matrix findMatrix(String name) {
        int j = 0;
        while (!name.equals(names[j])) {
            j++; //J will have the index for the correct matrix in the matrices array
        }
        return matrices[j];
    }

    public static void createMatrix () {
        int correct = 1;
        int curRow;
        int curCol;
        double tempval;

        if (curLetter == 26) {
            throw new ArrayStoreException("There is not enough memory to create a new matrix!"); //Is this the right type
        }
        JOptionPane.showMessageDialog(null, "You are creating matrix " + names[curLetter]);
        int rows = Integer.parseInt(JOptionPane.showInputDialog("How many rows does your matrix have?"));
        int cols = Integer.parseInt(JOptionPane.showInputDialog("How many columns does your matrix have?"));
        double[][] vals = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                vals[i][j] = Double.parseDouble(JOptionPane.showInputDialog("Enter value at " + (i + 1) + ", " + (j + 1)));
            }
        }

        Matrix matrixObj = new Matrix(vals);
        String matrix = "Is this correct? \n" + matrixObj.toString();
        correct = JOptionPane.showConfirmDialog(null, matrix);
        while (correct == 1) {
            curRow = Integer.parseInt(JOptionPane.showInputDialog("What row needs to be changed?"));
            curCol = Integer.parseInt(JOptionPane.showInputDialog("What column needs to be changed?"));
            tempval = Double.parseDouble(JOptionPane.showInputDialog("Enter value at " + curRow + ", " + curCol));
            matrixObj.setValue(curRow, curCol, tempval);
            matrix = "Is this correct? \n" + matrixObj.toString();
            correct = JOptionPane.showConfirmDialog(null, matrix);
        }

        matrices[curLetter] = matrixObj;
        curLetter++;
    }

    public static void savedMatrices() {
        String curMatrices = "";
        int j = 0;

        for (int i = 0; i < curLetter; i++) {
            curMatrices += names[i];
            if (i != curLetter - 1) {
                curMatrices += ", ";
            }
        }
        String name = JOptionPane.showInputDialog("Current Matrices: \n" + curMatrices + "\nWhich matrix do you want to view?");
        Matrix A = findMatrix(name);
        String matrix = A.toString();

        JOptionPane.showMessageDialog(null, "Matrix:\n" + matrix);
    }

    public static void matrixOperations () {
        //Todo: make this box vertical
        int choice = JOptionPane.showOptionDialog(null, "Calculate:", "Select one:", 0, 3, null, opOptions, opOptions[6]);
        String name = JOptionPane.showInputDialog("Which matrix do you want to use?");
        Matrix A = findMatrix(name);

        if (choice == 0) { //Det
            JOptionPane.showMessageDialog(null, "The determinant of your matrix is: " + LinearAlgebra.determinant(A));

        }
        else if (choice == 1) { //Eig
            double[] eVals =  LinearAlgebra.eig(A);
            String eigs = "";
            for (int i = 0; i < eVals.length; i++) {
                eigs += String.format("%.4f", eVals[i]);
                if (i != eVals.length - 1) {
                    System.out.println(",");
                }
                eigs += " ";
            }
            JOptionPane.showMessageDialog(null, "The eigenvalues of your matrix are: " + eigs);
        }
        else if (choice == 2) { //LU
            Matrix[] LU = LinearAlgebra.LUDecomp(A);
            String L = LU[0].toString();
            String U = LU[1].toString();
            JOptionPane.showMessageDialog(null, "Matrix L:\n" + L);
            JOptionPane.showMessageDialog(null, "Matrix U:\n" + U);
        }
        else if (choice == 3) { //Inv
            String Inv = LinearAlgebra.matrixInverse(A).toString();
            JOptionPane.showMessageDialog(null, "Inverse of your matrix:\n" + Inv);
        }
        else if (choice == 4) { //AB
            name = JOptionPane.showInputDialog("What other matrix do you want to use?");
            Matrix B = findMatrix(name);
            String mul = LinearAlgebra.multiplyMatrices(A, B).toString();
            JOptionPane.showMessageDialog(null, "Product of your matrices:\n" + mul);
        }
        else if (choice == 5) { //QR
            Matrix[] QR = LinearAlgebra.QRFactorization(A);
            String Q = QR[0].toString();
            String R = QR[1].toString();
            JOptionPane.showMessageDialog(null, "Matrix Q:\n" + Q);
            JOptionPane.showMessageDialog(null, "Matrix R:\n" + R);
        }
        else if (choice == 6) { //RREF
            String RREF = LinearAlgebra.RREF(A).toString();
            JOptionPane.showMessageDialog(null, "RREF of your matrix:\n" + RREF);
        }
    }

    public static void main(String[] args) {
        int choice;
        boolean exit = false;

        while (!exit) {
            //Todo: Make this box vertical
            choice = JOptionPane.showOptionDialog(null, "Options:", "Select one:", 0, 3, null, options, options[4]);
            if (choice == 0) {
                createMatrix();
            }
            else if (choice == 1) {
                matrixOperations();
            }
            else if (choice == 2) {
                savedMatrices();
            }
            else if (choice == 3) {

            }
            else {
                exit = true;
            }
        }
    }

}