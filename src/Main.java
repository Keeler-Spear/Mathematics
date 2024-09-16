public class Main {

    public static void main(String[] args) {
        double[][] vals = new double[3][3];
        vals[0][0] = 6;
        vals[0][1] = 2;
        vals[0][2] = 1;
        vals[1][0] = 2;
        vals[1][1] = 4;
        vals[1][2] = 4;
        vals[2][0] = 3;
        vals[2][1] = 6;
        vals[2][2] = 10;
        Matrix matrix = new Matrix(vals);

        double[][] vals2 = new double[3][4];
        vals2[0][0] = 1;
        vals2[0][1] = 2;
        vals2[0][2] = 1;
        vals2[0][3] = 3;
        vals2[1][0] = 0;
        vals2[1][1] = 1;
        vals2[1][2] = 0;
        vals2[1][3] = 4;
        vals2[2][0] = 0;
        vals2[2][1] = 0;
        vals2[2][2] = 0;
        vals2[2][3] = 0;
        Matrix matrix2 = new Matrix(vals2);

        double[][] vals3 = new double[2][2];
        vals3[0][0] = 3;
        vals3[0][1] = 4;
        vals3[1][0] = 5;
        vals3[1][1] = 6;
        Matrix matrix3 = new Matrix(vals3);

        System.out.println(matrix);
        System.out.println(LinearAlgebraCalculator.matrixInverse(matrix));

    }
}
