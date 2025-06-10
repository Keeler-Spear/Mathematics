import Mathematics.LinearAlgebra;
import Mathematics.Matrix;
import Mathematics.ODE;
import Mathematics.TriFunction;
import org.junit.jupiter.api.Test;

import java.util.function.Function;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class ODETest {

    private static final double TOL = 0.1;

    @Test
    public void BVPTest1() {

        TriFunction<Double, Double, Double, Double> ode = (t, y, yp) -> 8 * yp - 16 * y;
        Matrix approx = LinearAlgebra.vectorFromColumn(ODE.solveBVP(ode, 0, 1, 1, 0, 0.001), 1);
        Function<Double, Double> f = (t) -> Math.exp(4 * t) - t * Math.exp(4 * t);
        Matrix tLin = LinearAlgebra.linSpace(0, 1, 0.001);
        Matrix sol = LinearAlgebra.applyFunction(tLin, f);
        assertEquals(sol, approx);
    }

}
