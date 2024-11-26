//My own implementation of the Java Math package
public class Maths {

    public final static double PI = 3.14159265358979323846;
    public final static double e = 2.718281828459045;
    private static final double tol = 0.00000001;
    private static final int SIGFIG = 4;

    //Todo Implement myself; currently this method is written by ChatGPT 4.0.
    private static double round(double x) {
        // Step 1: Multiply by 1000 to shift the decimal point 3 places to the right
        x = x * pow(10, SIGFIG); //EDITED FROM CHATGPT to use SIGFIG

        // Step 2: Add 0.5 to implement rounding (positive or negative)
        if (x >= 0) {
            x = (long)(x + 0.5);
        } else {
            x = (long)(x - 0.5);
        }

        // Step 3: Divide by 1000 to shift the decimal point back to the left
        return x / 1000.0;
    }

    //Todo make it so the method can handle non-integer powers
    public static double pow(double x, int y) { //X^Y
        double result = 1.0;
        if (y >= 0) {
            for (int i = 0; i < y; i++) {
                result *= x;
            }
        }
        else { //Exponent must be negative
            for (int i = 0; i < abs(y); i++) {
                result /= x;
            }
        }
        return result;
    }

    public static double abs(double x) {
        double result = x;
        if (x < 0) {
            result = -x;
        }
        return result;
    }

    public static double sqrt(double x) {
        double xk = x/10; //Guess for the sqrt
        while (abs(xk * xk - x) > tol) {
            xk = xk - (xk*xk - x) / (2*xk); //Newton's Method. Find the roots of y = x^2 - x
        }
        return xk;
    }

    //TODO Finish the method
    //X ^ Y = Z --> Log(x)Z = Y
    public static double log (double x, double y) { //log(x)(y) = result
        double result = 0.0;
        System.out.println(result);
        return result;
    }

    public static double log (double x) {
        return log(10, x);
    }

    public static double ln(double x) {
        return log(e, x);
    }


}