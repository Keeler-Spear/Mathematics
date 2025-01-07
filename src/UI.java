import javax.swing.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * A main class that functions as an interface between a user and the LinearAlgebra and Matrix classes.
 * <p>
 *     This class supports matrix operations such as addition, multiplication, RREF, and more.
 * </p>
 *
 * @author Keeler Spear
 * @version %I%, %G%
 * @since 1.0
 */
public class UI {
    private static final String[] NAMES = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"};
    private static final String[] OPTIONS = {"Create Matrix", "Perform Matrix Operations", "View Saved Matrices", "Load Data", "Save Data", "Reset Data", "Exit"};
    private static final String[] OPOPTIONS = {"Determinant", "Eigenvalues", "LU Decomposition", "Matrix Inverse", "Matrix Multiplication", "QR Decomposition", "RREF"};
    private static final String SAVE_FILE = "SavedVars.csv";

    private static Matrix[] matrices = new Matrix[26];
    private static ArrayList<String> curNames = new ArrayList<>();

    //Lets a user select a saved matrix.
    private static Matrix selectMatrix() {
        int choice = JOptionPane.showOptionDialog(null, "Which matrix do you want to use?", "Select one:", 0, 3, null, curNames.toArray(), curNames.getFirst());

        return matrices[choice];
    }

    //Lets a user create and save a matrix.
    private static void createMatrix () {
        int correct;
        int curRow;
        int curCol;
        double tempval;

        if (curNames.size() >= 27) {
            throw new ArrayStoreException("There is not enough memory to create a new matrix!"); //Is this the right type
        }

        JOptionPane.showMessageDialog(null, "You are creating matrix " + NAMES[curNames.size()]);
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

        matrices[curNames.size()] = matrixObj;
        curNames.add(NAMES[curNames.size()]);
    }

    //Displays a list of saved matrices and lets the user view one.
    private static void viewSavedMatrices() {
        String curMatrices = "";

        for (int i = 0; i < curNames.size(); i++) {
            curMatrices += curNames.get(i);
            if (i != curNames.size() - 1) {
                curMatrices += ", ";
            }
        }

        Matrix A = selectMatrix();
        String matrix = A.toString();

        JOptionPane.showMessageDialog(null, "Matrix:\n" + matrix);
    }

    //Lets the user perform matrix operations on a saved matrix.
    private static void matrixOperations () {
        int choice = JOptionPane.showOptionDialog(null, "Calculate:", "Select one:", 0, 3, null, OPOPTIONS, OPOPTIONS[6]);
        Matrix A = selectMatrix();

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
            String P = LU[2].toString();
            JOptionPane.showMessageDialog(null, "Matrix L:\n" + L);
            JOptionPane.showMessageDialog(null, "Matrix U:\n" + U);
            JOptionPane.showMessageDialog(null, "Matrix P:\n" + P);

        }
        else if (choice == 3) { //Inv
            String Inv = LinearAlgebra.matrixInverse(A).toString();
            JOptionPane.showMessageDialog(null, "Inverse of your matrix:\n" + Inv);
        }
        else if (choice == 4) { //AB
            Matrix B = selectMatrix();
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

    //Loads the matrices from SAVE_FILE into the list of saved matrices, with their designations assigned in alphabetical
    //order starting from the first matrix in SAVE_FILE.
    private static void loadData() {
        try {
            Scanner in = new Scanner(new FileInputStream(SAVE_FILE));

            if (in.hasNextLine()) {
                String line = in.nextLine();
                int rows;
                int cols;
                double[][] temp = new double[100][100];
                double[][] vals;

                while (in.hasNextLine()) {
                    rows = 0;
                    cols = 0;
                    while(!line.equals("NEWMATRIX")) {
                        String data[] = line.split(" ");
                        for (int i = 0; i < data.length; i++) {
                            data[i] = data[i].replace("[", "");
                            data[i] = data[i].replace("]", "");
                            temp[rows][i] = Double.parseDouble(data[i]);
                        }
                        cols = data.length;
                        rows++;
                        line = in.nextLine();
                    }
                    vals = new double[rows][cols]; //Creates the matrix from the saved data
                    for (int i = 0; i < rows; i++) {
                        for (int j = 0; j < cols; j++) {
                            vals[i][j] = temp[i][j];
                        }
                    }
                    matrices[curNames.size()] = new Matrix(vals);;
                    curNames.add(NAMES[curNames.size()]);
                    if (in.hasNextLine()) {
                        line = in.nextLine();
                    }
                }
                in.close();
                JOptionPane.showMessageDialog(null, "Data Loaded!");
            }
        }
        catch (IOException ex) {
            System.out.println("Error: " + ex);
        }
    }

    //Saves the current matrices in alphabetical order, starting with A, to SAVE_FILE.
    private static void saveData() {
        try {
            PrintStream out = new PrintStream(new FileOutputStream(SAVE_FILE));

            for (int i = 0; i < curNames.size(); i++) {
                out.println(matrices[i]);
                out.println("NEWMATRIX");
            }

            out.close();
            JOptionPane.showMessageDialog(null, "Data Saved!");
        } catch (FileNotFoundException ex) {
            System.out.println("Error: " + ex);
        }
    }

    //Deleted all the data is SAVE_FILE.
    private static void resetData() {
        try {
            new FileOutputStream(SAVE_FILE).close();
            JOptionPane.showMessageDialog(null, "Data Reset!");
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * An interactive interface that allows the user to create matrices and perform linear algebra operations.
     * The following interface operations are supported:
     *         <ul>
     *             <li> Create Matrix</li>
     *             <li> Perform Matrix Operations</li>
     *             <li> View Saved Matrices</li>
     *             <li> Load Data</li>
     *             <li> Save Data</li>
     *             <li> Reset Data</li>
     *             <li> Exit</li>
     *         </ul>
     * The following mathematical calculations and operations are supported:
     *         <ul>
     *             <li> Determinant</li>
     *             <li> Eigenvalues</li>
     *             <li> LU Decomposition</li>
     *             <li> Matrix Inverse</li>
     *             <li> Matrix Multiplication</li>
     *             <li> QR Decomposition</li>
     *             <li> RREF</li>
     *         </ul>
     */
    public static void main(String[] args) {
        int choice;
        boolean exit = false;

        while (!exit) {
            choice = JOptionPane.showOptionDialog(null, "Options:", "Select one:", 0, 3, null, OPTIONS, OPTIONS[6]);

            if (choice == 0) {
                createMatrix();
            }
            else if (choice == 1) {
                if (curNames.size() > 0) {
                    matrixOperations();
                }
                else {
                    JOptionPane.showMessageDialog(null, "There are no matrices to use!");
                }
            }
            else if (choice == 2) {
                if (curNames.size() > 0) {
                    viewSavedMatrices();
                }
                else {
                    JOptionPane.showMessageDialog(null, "There are no saved matrices!");
                }
            }
            else if (choice == 3) {
                loadData();
            }
            else if (choice == 4) {
                saveData();
            }
            else if (choice == 5) {
                resetData();
            }
            else {
                exit = true;
            }
        }
    }

}