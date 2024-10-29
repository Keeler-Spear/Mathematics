//8th Programming Assignment for Numerical Analysis

public class PA8 {

    public static void RREFandLU(Matrix matrix, Matrix b) {
        Matrix[] LU = LinearAlgebraCalculator.LUDecomp(matrix);
        Matrix L = LU[0];
        Matrix U = LU[1];
        Matrix P = LU[2];
        long start = System.nanoTime();
        Matrix x1 = LinearAlgebraCalculator.LUSolve(matrix, b);
        long end = System.nanoTime();
        long time = end - start;
        System.out.println("LU Decomposition:");
        System.out.println("\nP:");
        System.out.println(P);
        System.out.println("L:");
        System.out.println(L);
        System.out.println("U:");
        System.out.println(U);
        System.out.println("LU Solution:");
        System.out.println(x1);
        System.out.println("LU Time (nanoseconds): " + time +"\n");
        matrix = LinearAlgebraCalculator.augmentMatrix(matrix, b);
        start = System.nanoTime();
        Matrix x2 = LinearAlgebraCalculator.RREFSolve(matrix);
        end = System.nanoTime();
        time = end - start;
        System.out.println("RREF Solution:");
        System.out.println(x2);
        System.out.println("RREF Time (nanoseconds): " + time +"\n");

    }

    public static void P1 () {
        double[][] vals = new double[3][3];
        vals[0][0] = 2;
        vals[0][1] = -1;
        vals[0][2] = 1;
        vals[1][0] = 3;
        vals[1][1] = 3;
        vals[1][2] = 9;
        vals[2][0] = 3;
        vals[2][1] = 3;
        vals[2][2] = 5;
        Matrix matrix1 = new Matrix(vals);

        double[][] v1 = new double[3][1];
        v1[0][0] = -1;
        v1[1][0] = 0;
        v1[2][0] = 4;
        Matrix b1 = new Matrix(v1);
        RREFandLU(matrix1, b1);
    }

    public static void P2 () {
        double[][] vals = new double[3][3];
        vals[0][0] = 1;
        vals[0][1] = -1;
        vals[0][2] = 0;
        vals[1][0] = 2;
        vals[1][1] = 2;
        vals[1][2] = 3;
        vals[2][0] = -1;
        vals[2][1] = 3;
        vals[2][2] = 2;
        Matrix matrix1 = new Matrix(vals);
        double[][] v1 = new double[3][1];
        v1[0][0] = 2;
        v1[1][0] = -1;
        v1[2][0] = 4;
        Matrix b1 = new Matrix(v1);
        RREFandLU(matrix1, b1);
    }

    public static void P3 () {
        double[][] vals = new double[4][4];
        vals[0][0] = 2.1756;
        vals[0][1] = 4.0231;
        vals[0][2] = -2.1732;
        vals[0][3] = 5.1967;
        vals[1][0] = -4.0231;
        vals[1][1] = 6;
        vals[1][2] = 0;
        vals[1][3] = 1.1973;
        vals[2][0] = -1;
        vals[2][1] = -5.2107;
        vals[2][2] = 1.1111;
        vals[2][3] = 0;
        vals[3][0] = 6.0235;
        vals[3][1] = 7;
        vals[3][2] = 0;
        vals[3][3] = -4.1561;
        Matrix matrix1 = new Matrix(vals);
//
        double[][] v1 = new double[4][1];
        v1[0][0] = 17.102;
        v1[1][0] = -6.1593;
        v1[2][0] = 3.0004;
        v1[3][0] = 0;
        Matrix b1 = new Matrix(v1);
        RREFandLU(matrix1, b1);
    }

    public static void P4 () {
        double[][] vals = new double[4][4];
        vals[0][0] = 2;
        vals[0][1] = 1;
        vals[0][2] = 0;
        vals[0][3] = 0;
        vals[1][0] = -1;
        vals[1][1] = 3;
        vals[1][2] = 3;
        vals[1][3] = 0;
        vals[2][0] = 2;
        vals[2][1] = -2;
        vals[2][2] = 1;
        vals[2][3] = 4;
        vals[3][0] = -2;
        vals[3][1] = 2;
        vals[3][2] = 2;
        vals[3][3] = 5;
        Matrix matrix1 = new Matrix(vals);

        double[][] v1 = new double[4][1];
        v1[0][0] = 0;
        v1[1][0] = 5;
        v1[2][0] = -2;
        v1[3][0] = 6;
        Matrix b1 = new Matrix(v1);
        RREFandLU(matrix1, b1);
    }

    public static void main(String[] args) {
        P4();
    }
}
