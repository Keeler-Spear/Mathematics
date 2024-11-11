import java.util.Scanner;

public class UI {
    private int letter = 0;
    private final String[] names = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"};
    private String[] matrices = new String[26];

    private static final String SAVE_FILE = "SavedVars.csv";

    public static void createMatrix () {}

    public static void saveMatrix () {}

    public static void savedMatrices() {}

    public static void matrixOperations () {}

    public static void main(String[] args) {
        Scanner scnr = new Scanner(System.in);

        boolean exit = false;
        while (!exit) {
            System.out.println("Options:\n" + "(1) Create Matrix\n" + "(2) Save Matrix\n" + "(3) View Saved Matrices\n" + "(4) Perform Matrix Operations\n" + "(5) Exit");
            String option = scnr.nextLine();
            if (option.equals("1")) {
                createMatrix();
            }
            else if (option.equals("2")) {
                saveMatrix();
            }
            else if (option.equals("3")) {
                savedMatrices();
            }
            else if (option.equals("4")) {
                matrixOperations();
            }
            else if (option.equals("5")) {
                exit = true;
            }
        }
    }
}
