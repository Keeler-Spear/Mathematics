public class Main {

    public static void main(String[] args) {
        double[][] vals = new double[3][4];
        vals[0][0] = 7;
        vals[0][1] = 0;
        vals[0][2] = 1;
        vals[0][3] = 4;
        vals[1][0] = 1;
        vals[1][1] = 9;
        vals[1][2] = 4;
        vals[1][3] = 2;
        vals[2][0] = 6;
        vals[2][1] = 2;
        vals[2][2] = 10;
        vals[2][3] = 3;


        Matrix matrix = new Matrix(vals);
        System.out.println(matrix);
        matrix.removeCol(1);
        matrix.removeRow(1);
        System.out.println(matrix);

    }
}
