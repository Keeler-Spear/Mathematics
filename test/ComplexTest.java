import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

//Todo: Add more tests
class ComplexTest {

    final double h = 0.1;

    @Test
    void testAdd() {
        ComplexNumber z0 = new ComplexNumber(3, 4);
        ComplexNumber z1 = new ComplexNumber(-2, 7);
        ComplexNumber z = new ComplexNumber(1, 11);
        assertEquals(z, Complex.add(z0, z1));
    }

    @Test
    void testSub() {
        ComplexNumber z0 = new ComplexNumber(3, 4);
        ComplexNumber z1 = new ComplexNumber(-2, 7);
        ComplexNumber z = new ComplexNumber(5, -3);
        assertEquals(z, Complex.sub(z0, z1));
    }

    @Test
    void testMul() {
        ComplexNumber z0 = new ComplexNumber(3, 4);
        ComplexNumber z1 = new ComplexNumber(-2, 7);
        ComplexNumber z = new ComplexNumber(-34, 13);
        assertEquals(z, Complex.mul(z0, z1));
    }

    @Test
    void testDiv() {
        ComplexNumber z0 = new ComplexNumber(3, 4);
        ComplexNumber z1 = new ComplexNumber(-2, 7);
        ComplexNumber z = new ComplexNumber(22.0/53.0, -29.0/53.0);
        assertEquals(z, Complex.div(z0, z1));
    }

    @Test
    void testConjugate0() {
        ComplexNumber z0 = new ComplexNumber(3, 4);
        ComplexNumber z = new ComplexNumber(3, -4);
        assertEquals(z, Complex.conjugate(z0));
    }

    @Test
    void testConjugate1() {
        ComplexNumber z0 = new ComplexNumber(-2, 7);
        ComplexNumber z = new ComplexNumber(-2, -7);
        assertEquals(z, Complex.conjugate(z0));
    }

    @Test
    void testModulus0() {
        ComplexNumber z0 = new ComplexNumber(3, -4);
        assertEquals(5, Complex.modulus(z0));
    }

    @Test
    void testModulus1() {
        ComplexNumber z0 = new ComplexNumber(-2, 7);
        assertEquals(Math.sqrt(53), Complex.modulus(z0));
    }


}