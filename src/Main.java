import javax.sound.sampled.Line;
import java.sql.SQLOutput;

public class Main {

    public static void main(String[] args) {
        double[][] vals = new double[3][3];
        vals[0][0] = 1;
        vals[0][1] = 2;
        vals[0][2] = 3;
        vals[1][0] = 4;
        vals[1][1] = 5;
        vals[1][2] = 6;
        vals[2][0] = 7;
        vals[2][1] = 8;
        vals[2][2] = 12;
        Matrix matrix1 = new Matrix(vals);

        double[][] vals2 = new double[4][4];
        vals2[0][0] = 3;
        vals2[0][1] = -7;
        vals2[0][2] = -2;
        vals2[0][3] = 2;
        //vals2[0][4] = -2;
        vals2[1][0] = -3;
        vals2[1][1] = 5;
        vals2[1][2] = 1;
        vals2[1][3] = 0;
        //vals2[1][4] = 1;
        vals2[2][0] = 6;
        vals2[2][1] = -4;
        vals2[2][2] = 0;
        vals2[2][3] = -5;
        //vals2[2][4] = 8;
        vals2[3][0] = -9;
        vals2[3][1] = 5;
        vals2[3][2] = -5;
        vals2[3][3] = 12;
        //vals2[3][4] = 1;
        Matrix matrix2 = new Matrix(vals2);

        double[][] vals3 = new double[3][3];
        vals3[0][0] = 3;
        vals3[0][1] = 4;
        vals3[0][2] = -1;
        vals3[1][0] = 5;
        vals3[1][1] = 6;
        vals3[1][2] = -2;
        Matrix matrix3 = new Matrix(vals3);

        double[][] vals4 = new double[4][1];
        vals4[0][0] = -9;
        vals4[1][0] = 5;
        vals4[2][0] = 7;
        vals4[3][0] = 11;
        Matrix matrix4 = new Matrix(vals4);

        //System.out.println(matrix2);
        Matrix x1 = LinearAlgebraCalculator.LUSolve(matrix2, matrix4);
        matrix2 = LinearAlgebraCalculator.augmentMatrix(matrix2, matrix4);
        //System.out.println(matrix2);
        Matrix x2 = LinearAlgebraCalculator.RREF(matrix2);
        System.out.println(x1);
        System.out.println(x2);



    }
}
