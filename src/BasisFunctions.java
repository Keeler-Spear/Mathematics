import java.util.function.Function;

/**
 * A class that creates basis functions.
 * <p>
 *     This class and all functions within were added after I submitted my application to your institution.
 * </p>
 *
 * @author Keeler Spear
 * @version %I%, %G%
 * @since 1.0
 */
public class BasisFunctions {

    /**
     * Creates a set of polynomial basis functions of the size specified.
     *
     * @param n The highest polynomial order in the basis function set.
     * @throws IllegalArgumentException If the set's order is negative.
     */
    public static Function[] polynomials(int n) {
        if (n < 0) {
            throw new IllegalArgumentException("The set's order must be at least zero!");
        }

        Function<Double, Double>[] fncs = new Function[n + 1];

        for (int i = 0; i <= n; i++) {
            double I = i;
            fncs[i] = x -> Math.pow(x, I);
        }

        return fncs;
    }

    /**
     * Creates a set of trigonometric basis functions of the size specified.
     *
     * @param n The highest polynomial order in the basis function set.
     * @throws IllegalArgumentException If the set's order is negative.
     */
    public static Function[] trig(int n) { //Orthogonal on [-PI, PI]
        if (n < 0) {
            throw new IllegalArgumentException("The set's order must be at least zero!");
        }

        Function<Double, Double>[] fncs = new Function[n + 1];

        fncs[0] = x -> 1.0;

        if (n >= 1) {
            fncs[1] = x -> Math.sin(x);
        }

        if (n >= 2) {
            fncs[2] = x -> Math.cos(x);
        }

        for (int i = 3; i <= n; i++) {
            int I = i;
            if (i % 2 == 1) { //i is odd
                fncs[i] = x -> Math.sin(((I + 1.0) / 2.0) * x);
            }
            else { //i is even
                fncs[i] = x -> Math.cos((I / 2.0) * x);
            }
        }

        return fncs;
    }

    public static Function[] legendre(int n) {
        if (n < 0) {
            throw new IllegalArgumentException("The set's order must be at least zero!");
        }

        Function<Double, Double>[] fncs = new Function[n + 1];

        fncs[0] = x -> 1.0;

        if (n >= 1) {
            fncs[1] = x -> x;
        }

        for (int i = 2; i <= n; i++) {
            int I = i;
            fncs[i] = x -> x * ((2 * I + 1) / (I + 1)) * fncs[I - 1].apply(x) - (I / (I + 1)) * fncs[I - 2].apply(x);
        }

        return fncs;
    }

}
