package Mathematics;//WORK IN PROGRESS
//WORK IN PROGRESS
//WORK IN PROGRESS
//WORK IN PROGRESS
//WORK IN PROGRESS

//My own implementation of the Java Math package. It is currently not in use.
public class Maths {

    public final static double PI = 3.14159265358979323846;
    public final static double e = 2.718281828459045;
    public static final double tol = 0.00000001;
    public static final int SIGFIG = 4;

    private static double round(double x) {
        x = x * pow(10, SIGFIG);

        if (x >= 0) {
            x = (long)(x + 0.5);
        } else {
            x = (long)(x - 0.5);
        }

        return x / 1000.0;
    }

    //Todo make it so the method can handle non-integer powers
    //Break the number into an integer and a decimal (x/1000 or something) to put into faction form, then recursivly call this function
    public static double pow(double x, double n) { //x^n
        double result = 1.0;
        if (n - (int) n < tol) { //If n is an integer
            if (n >= 0) {
                for (int i = 0; i < n; i++) {
                    result *= x;
                }
            }
            else { //Exponent must be negative
                for (int i = 0; i < abs(n); i++) {
                    result /= x;
                }
            }
            return result;
        }
        //Converting n into a fraction\
        String num = Double.toString(abs(n));
        int length = num.indexOf('.');
        length = num.length() - length - 1;

        int denominator = (int) pow(10, length);
        int numerator = (int) (n * denominator);
        result = nthRoot(x, denominator);
        result = pow(result, numerator);
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
        return nthRoot(x, 2);
    }

    public static double nthRoot(double x, int n) {
        double xk = x/10; //Guess for the root
        if (n < 0) {
            throw new IllegalArgumentException("n must be positive!");
        }
        while (abs(pow(xk, n) - x) > tol) {
            xk = xk - (pow(xk, n) - x) / (n*pow(xk, n-1)); //Newton's Method. Find the roots of y = x^n - root^n
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
