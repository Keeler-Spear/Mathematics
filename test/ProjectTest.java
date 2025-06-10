
import Mathematics.Calculus;
import Mathematics.LinearAlgebra;
import Mathematics.Matrix;
import org.junit.jupiter.api.Test;

import java.util.function.Function;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class ProjectTest {

    private static final double tol = 0.0001;

    //All matrices used in these tests are from Linear Algebra and Its Applications, Sixth Edition by Lay, Lay, and McDonald.

    @Test //Section 2.1.2 Example 1
    public void additionTest() {
        double[][] valsA = {
                {1, 2, 1},
                {2, 1, 2},
                {1, 0, 1}
        };
        Matrix A = new Matrix(valsA);

        double[][] valsB = {
                {1, -5, 2},
                {2, -1, 0},
                {1, 0, 2}
        };
        Matrix B = new Matrix(valsB);

        Matrix C = LinearAlgebra.addMatrices(A, B, 1);

        double[][] correctVals = {
                {2, -3, 3},
                {4, 0, 2},
                {2, 0, 3}
        };
        Matrix correctMatrix = new Matrix(correctVals);

        assertEquals(correctMatrix, C);
    }

    @Test //Section 2.1.4 Example 7
    public void multiplicationTest() {
        double[][] valsA = {
                {1, 2, 3},
                {4, 5,6},
        };
        Matrix A = new Matrix(valsA);

        double[][] valsB = {
                {7, 8, 9, 10},
                {11, 12, 13, 14},
                {15, 16, 17, 18}
        };
        Matrix B = new Matrix(valsB);

        Matrix C = LinearAlgebra.multiplyMatrices(A, B);

        double[][] correctVals = {
                {74, 80, 86, 92},
                {173, 188, 203, 218},
        };
        Matrix correctMatrix = new Matrix(correctVals);

        assertEquals(correctMatrix, C);
    }

    @Test //Section 1.2.7 Reasonable Answers
    public void RREFTest() {
        double[][] valsA = {
                {1, -2, 1, 2},
                {1, -1, 2, 5},
                {0, 1, 1, 3}
        };
        Matrix A = new Matrix(valsA);

        Matrix C = LinearAlgebra.RREF(A);

        double[][] correctVals = {
                {1, 0, 3, 8},
                {0, 1, 1, 3},
                {0, 0, 0, 0}
        };
        Matrix correctMatrix = new Matrix(correctVals);

        assertEquals(correctMatrix, C);
    }

    @Test //Section 2.5.2 Example 1
    public void LUTest() {
        double[][] valsA = {
                {3, -7, -2, 2},
                {-3, 5, 1, 0},
                {6, -4, 0, -5},
                {-9, 5, -5, 12}
        };
        Matrix A = new Matrix(valsA);

        double[][] bVals = {
                {-9},
                {5},
                {7},
                {11}
        };
        Matrix b = new Matrix(bVals);

        Matrix[] LU = LinearAlgebra.LUDecomp(A);

        Matrix L = LU[0];
        Matrix U = LU[1];
        Matrix P = LU[2];

        Matrix x = LinearAlgebra.LUSolve(A, b);

        double[][] correctxVals = {
                {3},
                {4},
                {-6},
                {-1}
        };
        Matrix correctx = new Matrix(correctxVals);

        assertEquals(LinearAlgebra.multiplyMatrices(P, A), LinearAlgebra.multiplyMatrices(L, U)); //Is PA = LU

        assertEquals(correctx, x);
    }

    @Test //Section 2.2.4 Example 7
    public void InverseTest() {
        double[][] valsA = {
                {0, 1, 2},
                {1, 0, 3},
                {4, -3, 8}
        };
        Matrix A = new Matrix(valsA);

        Matrix C = LinearAlgebra.matrixInverse(A);

        double[][] correctVals = {
                {-9/2.0, 7, -3/2.0},
                {-2, 4, -1},
                {3/2.0, -2.0, 1/2.0}
        };
        Matrix correctMatrix = new Matrix(correctVals);

        assertEquals(correctMatrix, C);
    }

    @Test //Section 3.1.3 Practice Problem 1
    public void DeterminantTest() {
        double[][] valsA = {
                {5, -7, 2, 2},
                {0, 3, 0, -4},
                {-5, -8, 0, 3},
                {0, 5, 0, -6}
        };
        Matrix A = new Matrix(valsA);

        double det = LinearAlgebra.determinant(A);

        double trueDet = 20;

        assertEquals(trueDet, det, tol);
    }

    @Test //Section 6.4.3 Example 4
    public void QRTest() {
        double[][] valsA = {
                {1, -3, 0},
                {1, 1, -2/3.0},
                {1, 1, 1/3.0},
                {1, 1, 1/3.0}
        };
        Matrix A = new Matrix(valsA);

        Matrix[] QR = LinearAlgebra.QRFactorization(A);

        Matrix Q = QR[0];
        Matrix R = QR[1];

        assertEquals(A, LinearAlgebra.multiplyMatrices(Q, R)); //Is A = QR
    }

    @Test //Section 5.6.7 Practice Problem 1
    public void EigTest() {
        double[][] valsA = {
                {7, -2, 0},
                {-2, 6, 2},
                {0, 2, 5}
        };
        Matrix A = new Matrix(valsA);
        A = LinearAlgebra.scaleMatrix(A, 1/9.0);

        double[] eigs = LinearAlgebra.eig(A);

        double[] trueEigs = {1, 2/3.0, 1/3.0};

        assertArrayEquals(trueEigs, eigs, tol);

    }

    //The functions use in the following texts is e^x^2, which is set in Function.java. The true values were calculated
    //with WolframAlpha

    @Test
    public void diffTest() {
        Function<Double, Double> fnc = x -> Math.pow(Math.E, Math.pow(x, 2.0));

        double der = Calculus.differentiate(fnc, 2);

        double trueDer = 218.3934009;

        assertEquals(trueDer, der, tol);
    }

    @Test
    public void intTest() {
        Function<Double, Double> fnc = x -> Math.pow(Math.E, Math.pow(x, 2.0));

       double area = Calculus.integrate(fnc, 0, 2);

       double trueArea = 16.4526;

        assertEquals(trueArea, area, tol);
    }
}
