import javax.swing.*;

public class UI {
    private static int letter = 0;
    private static final String[] names = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"};
    private static Matrix[] matrices = new Matrix[26];
    private static String options = "Options:\n" + "(1) Create Matrix\n" + "(2) Perform Matrix Operations\n" + "(3) View Saved Matrices\n" + "(4) TO BE ADDED\n" + "(5) Exit\n";

    private static final String SAVE_FILE = "SavedVars.csv";

    public static void createMatrix () {
        int correct = 1;
        int curRow;
        int curCol;
        double tempval;

        if (letter == 26) {
            throw new ArrayStoreException("There is not enough memory to create a new matrix!"); //Is this the right type
        }
        JOptionPane.showMessageDialog(null, "You are creating matrix " + names[letter]);
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

        matrices[letter] = matrixObj;
        letter++;
    }

    public static void savedMatrices() {
        String curMatrices = "";
        int j = 0;

        for (int i = 0; i < letter; i++) {
            curMatrices += names[i] + ", ";
        }
        String name = JOptionPane.showInputDialog("Current Matrices: \n" + curMatrices + "\nWhich Matrix do you want to view?");
        while (!name.equals(names[j])) {
            j++; //J will have the index for the correct matrix in the matrices array
        }
        Matrix A = matrices[j];
        String matrix = A.toString();

        JOptionPane.showMessageDialog(null, "Matrix " + names[j] + ": \n" + matrix);
    }

    public static void matrixOperations () {}

    public static void main(String[] args) {

        System.out.println(options);
        int choice = 1;
        boolean exit = false;

        while (!exit) {
            choice = Integer.parseInt(
                    JOptionPane.showInputDialog(options));
            if (choice == 1) {
                createMatrix();
            }
            else if (choice == 2) {
                matrixOperations();
            }
            else if (choice == 3) {
                savedMatrices();
            }
            else if (choice == 4) {

            }
            else if (choice == 5) {
                exit = true;
            }
        }
    }

}
