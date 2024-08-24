public class Main {

    public static void main(String[] args) {
        double[][] vals = new double[3][2];
        vals[0][0] = 7;
        vals[0][1] = 0;
//        vals[0][2] = 1;
        vals[1][0] = 1;
        vals[1][1] = 9;
//        vals[1][2] = 4;
        vals[2][0] = 6;
        vals[2][1] = 2;
//        vals[2][2] = 10;
        Matrix matrix = new Matrix(vals);

        double[][] vals2 = new double[3][4];
        vals2[0][0] = 8;
        vals2[0][1] = 9;
        vals2[0][2] = 8;
        vals2[0][3] = 7;
        vals2[1][0] = 2;
        vals2[1][1] = 9;
        vals2[1][2] = 1;
        vals2[1][3] = 4;
        vals2[2][0] = 8;
        vals2[2][1] = 9;
        vals2[2][2] = 5;
        vals2[2][3] = 4;
        Matrix matrix2 = new Matrix(vals2);

        System.out.println(matrix);
        matrix = LinearAlgebraCalculator.RREF(matrix);
        System.out.println(matrix);

    }
}
