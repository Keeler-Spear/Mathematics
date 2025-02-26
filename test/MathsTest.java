import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

//ToDo: Do relative error comparisons
class MathsTest {

    final double h = 0.000001;

    @org.junit.Test
    void powTest() {
        assertEquals(1, Maths.pow(143, 0), h);
        assertEquals(9, Maths.pow(3, 2), h);
        assertEquals(512, Maths.pow(2, 9), h);
        assertEquals(0.1, Maths.pow(10, -1), h);
        assertEquals(0.000540833, Maths.pow(43, -2), h);
    }
    @Test
    void powDecTest() {
        assertEquals(Math.pow(3, 1.5), Maths.pow(3, 1.5), h);
        assertEquals(Math.pow(4, 3.2), Maths.pow(4, 3.2), h);
        assertEquals(Math.pow(43, 17.2), Maths.pow(43, 17.2), h);

    }
    @Test
    void absTest() {
        assertEquals(0, Maths.abs(0), h);
        assertEquals(4, Maths.abs(4), h);
        assertEquals(0.37, Maths.abs(-0.37), h);
        assertEquals(4723, Maths.abs(4723), h);
        assertEquals(10274, Maths.abs(-10274), h);
    }
    @Test
    void rootTest() {
        assertEquals(0, Maths.sqrt(0), h);
        assertEquals(4, Maths.sqrt(16), h);
        assertEquals(244.810130509, Maths.sqrt(59932), h);
        assertEquals(0, Maths.nthRoot(0, 9), h);
        assertEquals(2, Maths.nthRoot(8, 3), h);
        assertEquals(5.22443184738, Maths.nthRoot(745, 4), h);
        Exception exception = assertThrows(IllegalArgumentException.class, () -> {Maths.nthRoot(3, -1);});
        assertEquals("n must be positive!", exception.getMessage());
    }
}