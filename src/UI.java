import javax.swing.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Scanner;

public class UI {
    private static final String[] NAMES = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"};
    private static final String[] OPTIONS = {"Create Matrix", "Perform Matrix Operations", "View Saved Matrices", "Load Data", "Save Data", "Reset Data", "Exit"};
    private static final String[] OPOPTIONS = {"Determinant", "Eigenvalues", "LU Decomposition", "Matrix Inverse", "Matrix Multiplication", "QR Decomposition", "RREF"};
    private static final String SAVE_FILE = "SavedVars.csv";

    private static Matrix[] matrices = new Matrix[26];
    private static ArrayList<String> curNames = new ArrayList<>();


    private static Matrix selectMatrix() {
        int choice = JOptionPane.showOptionDialog(null, "Which matrix do you want to use?", "Select one:", 0, 3, null, curNames.toArray(), curNames.getFirst());
        return matrices[choice];
    }

    public static void createMatrix () {
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

    public static void viewSavedMatrices() {
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

    public static void matrixOperations () {
        //Todo: make this box vertical
        int choice = JOptionPane.showOptionDialog(null, "Calculate:", "Select one:", 0, 3, null, OPOPTIONS, OPOPTIONS[6]);
        Matrix A = selectMatrix();
        //ToDo: Make an option to save these matrices
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

    private static void loadData() { //Loads the matrices from SAVE_FILE
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
                    //Creates the matrix from the saved data
                    vals = new double[rows][cols];
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

    private static void saveData() { //Saves the current matrices in alphabetical order to SAVE_FILE
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

    private static void resetData() {
        try {
            new FileOutputStream(SAVE_FILE).close();
            JOptionPane.showMessageDialog(null, "Data Reset!");
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static void main(String[] args) {
        int choice;
        boolean exit = false;
        while (!exit) {
            //Todo: Make this box vertical
            choice = JOptionPane.showOptionDialog(null, "Options:", "Select one:", 0, 3, null, OPTIONS, OPTIONS[6]);
            if (choice == 0) {
                createMatrix();
            }
            else if (choice == 1) {
                if (curNames.size() > 0) {
                    matrixOperations();
                }
                else {
                    JOptionPane.showMessageDialog(null, "There are matrices to use!");
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