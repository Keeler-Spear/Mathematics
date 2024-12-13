import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

//ToDo: Add more tests
class LinearAlgebraTest {

    final double h = 0.01;

    @Test
    void testRREF0() {
        double[][] vals = new double[4][5];
        vals[0][0] = 3;
        vals[0][1] = -13;
        vals[0][2] = 9;
        vals[0][3] = 3;
        vals[0][4] = -19;

        vals[1][0] = -6;
        vals[1][1] = 4;
        vals[1][2] = 1;
        vals[1][3] = -18;
        vals[1][4] = -34;

        vals[2][0] = 6;
        vals[2][1] = -2;
        vals[2][2] = 2;
        vals[2][3] = 4;
        vals[2][4] = 16;

        vals[3][0] = 12;
        vals[3][1] = -8;
        vals[3][2] = 6;
        vals[3][3] = 10;
        vals[3][4] = 26;
        Matrix matrix = new Matrix(vals);
        matrix = LinearAlgebra.RREF(matrix);

        Matrix correctMatrix = LinearAlgebra.identityMatrix(4);
        double[][] xVals = new double[4][1];
        xVals[0][0] = 3;
        xVals[1][0] = 1;
        xVals[2][0] = -2;
        xVals[3][0] = 1;
        Matrix x = new Matrix(xVals);

        correctMatrix = LinearAlgebra.augmentMatrix(correctMatrix, x);
        assertEquals(correctMatrix, matrix);
    }

    //The following tests were written by ChatGPT-4o and verified by me
    @Test
    public void testRREF1() {
        double[][] vals = {
                {2, 1, -1, 8},
                {-3, -1, 2, -11},
                {-2, 1, 2, -3}
        };
        Matrix matrix = new Matrix(vals);
        matrix = LinearAlgebra.RREF(matrix);

        double[][] correctVals = {
                {1, 0, 0, 2},
                {0, 1, 0, 3},
                {0, 0, 1, -1}
        };
        Matrix correctMatrix = new Matrix(correctVals);
        assertEquals(correctMatrix, matrix);
    }

    @Test
    public void testRREF2() {
        double[][] vals = {
                {1, 2, 3, 4},
                {2, 4, 6, 8},
                {3, 6, 9, 12}
        };
        Matrix matrix = new Matrix(vals);
        matrix = LinearAlgebra.RREF(matrix);

        double[][] correctVals = {
                {1, 2, 3, 4},
                {0, 0, 0, 0},
                {0, 0, 0, 0}
        };
        Matrix correctMatrix = new Matrix(correctVals);
        assertEquals(correctMatrix, matrix);
    }

    @Test
    public void testRREF3() {
        double[][] vals = {
                {1, -2, 3},
                {2, -4, 6},
                {3, -6, 9}
        };
        Matrix matrix = new Matrix(vals);
        matrix = LinearAlgebra.RREF(matrix);

        double[][] correctVals = {
                {1, -2, 3},
                {0, 0, 0},
                {0, 0, 0}
        };
        Matrix correctMatrix = new Matrix(correctVals);
        assertEquals(correctMatrix, matrix);
    }

    @Test
    public void testRREF4() {
        double[][] vals = {
                {1, 3, 1},
                {1, 1, 1},
                {3, 7, 3}
        };
        Matrix matrix = new Matrix(vals);
        matrix = LinearAlgebra.RREF(matrix);

        double[][] correctVals = {
                {1, 0, 1},
                {0, 1, 0},
                {0, 0, 0}
        };
        Matrix correctMatrix = new Matrix(correctVals);
        assertEquals(correctMatrix, matrix);
    }

    @Test
    public void testRREF5() {
        double[][] vals = {
                {1, 2, 1, -1},
                {2, 4, -1, -2},
                {-1, -2, 3, 4}
        };
        Matrix matrix = new Matrix(vals);
        matrix = LinearAlgebra.RREF(matrix);

        double[][] correctVals = {
                {1, 2, 0, 0},
                {0, 0, 1, 0},
                {0, 0, 0, 1}
        };
        Matrix correctMatrix = new Matrix(correctVals);
        assertEquals(correctMatrix, matrix);
    }

    @Test
    public void testRREF6() {
        double[][] vals = {
                {2, 4, -2, 2},
                {-4, -8, 6, -4},
                {1, 2, -1, 1}
        };
        Matrix matrix = new Matrix(vals);
        matrix = LinearAlgebra.RREF(matrix);

        double[][] correctVals = {
                {1, 2, 0, 1},
                {0, 0, 1, 0},
                {0, 0, 0, 0}
        };
        Matrix correctMatrix = new Matrix(correctVals);
        assertEquals(correctMatrix, matrix);
    }

    @Test
    public void testRREF7() {
        double[][] vals = {
                {1, 1, 1, 1},
                {0, 1, -1, -1},
                {2, 3, 0, 1}
        };
        Matrix matrix = new Matrix(vals);
        matrix = LinearAlgebra.RREF(matrix);

        double[][] correctVals = {
                {1, 0, 0, 2},
                {0, 1, 0, -1},
                {0, 0, 1, 0}
        };
        Matrix correctMatrix = new Matrix(correctVals);
        assertEquals(correctMatrix, matrix);
    }

    @Test
    public void testRREF8() {
        double[][] vals = {
                {1, 2, 3, 4, 5},
                {0, 1, 2, 3, 4},
                {0, 0, 1, 2, 3}
        };
        Matrix matrix = new Matrix(vals);
        matrix = LinearAlgebra.RREF(matrix);

        double[][] correctVals = {
                {1, 0, 0, 0, 0},
                {0, 1, 0, -1, -2},
                {0, 0, 1, 2, 3}
        };
        Matrix correctMatrix = new Matrix(correctVals);
        assertEquals(correctMatrix, matrix);
    }

    @Test
    public void testRREF9() {
        double[][] vals = {
                {1, -1, 2, -3},
                {3, -3, 6, -9},
                {1, -2, 3, -4}
        };
        Matrix matrix = new Matrix(vals);
        matrix = LinearAlgebra.RREF(matrix);

        double[][] correctVals = {
                {1, -1, 2, -3},
                {0, 0, 0, 0},
                {0, 0, 0, 0}
        };
        Matrix correctMatrix = new Matrix(correctVals);
        assertEquals(correctMatrix, matrix);
    }

    @Test
    public void testRREF10() {
        double[][] vals = {
                {1, 2, 3},
                {4, 5, 6},
                {7, 8, 10}
        };
        Matrix matrix = new Matrix(vals);
        matrix = LinearAlgebra.RREF(matrix);

        double[][] correctVals = {
                {1, 0, 0},
                {0, 1, 0},
                {0, 0, 1}
        };
        Matrix correctMatrix = new Matrix(correctVals);
        assertEquals(correctMatrix, matrix);
    }

    @Test
    public void testLU0() {
        double[][] vals = {
                {1, 2, 3},
                {4, 5, 6},
                {7, 8, 12}
        };
        Matrix matrix = new Matrix(vals);
        Matrix[] LU = LinearAlgebra.LUDecomp(matrix);
        Matrix L = LU[0];
        Matrix U = LU[1];

        double[][] LVALS = {
                {1, 0, 0},
                {4, 1, 0},
                {7, 2, 1}
        };
        Matrix correctL = new Matrix(LVALS);

        double[][] UVALS = {
                {1, 2, 3},
                {0, -3, -6},
                {0, 0, 3}
        };
        Matrix correctU = new Matrix(LVALS);

        assertEquals(correctL, L);

        assertEquals(correctU, U);
    }
}