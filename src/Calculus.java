//This class will solve single-variable calculus problems NUMERICALLY

public class Calculus {
    private static final double TOL = 0.00000001;
    private static final double hDEF = 0.001;

    private static double centeredDifference (double f1, double fm1, double h) { //O(h^2). f1 = f(x0+h), fm1 = f(x0-h)
        return (f1 - fm1) / (2 * h);
    }

    private static double leftSidedThreePoint (double f0, double f1, double f2, double h) {//O(h^2).  f0 = f(x0), f1 = f(x0+h), f2 = f(x0+2h).
        return (4*f1 - f2 - 3*f0) / (2 * h);
    }

    private static double rightSidedThreePoint (double f0, double fm1, double fm2, double h) {//O(h^2).  f0 = f(x0), fm1 = f(x0-h), fm2 = f(x0-2h).
        return (4*fm1 - fm2 - 3*f0) / (-2 * h);
    }

    public static double differentiate(Function function, double x0) { //Returns the derivative of a known function at x0 with O(hDEF^2) error
        return centeredDifference(function.eval(x0 + hDEF), function.eval(x0 - hDEF), hDEF);
    }

    public static double[] differentiate(double[] x, double[] y) { //Returns derivatives at each data point with O(h^2) error. This method assumes equally spaced data
        if (x.length != y.length) {
            throw new IllegalArgumentException("You must have the same number of data points in x and f(x)!");
        }
        if (x.length < 3) {
            throw new IllegalArgumentException("There are not enough data points to calculate the derivative!");
        }
        double[] derivative = new double[x.length];
        double h = x[1] - x[0];

        derivative[0] = leftSidedThreePoint(y[0], y[1], y[2], h);
        for (int i = 1; i < x.length - 1; i++) {
            derivative[i] = centeredDifference(y[i+1], y[i-1], h);
        }
        derivative[derivative.length - 1] = rightSidedThreePoint(y[y.length-1], y[y.length-2], y[y.length-3], h);
        return derivative;
    }

    public static double integrate(double[] y, double a, double b) { //Returns the integral with O(h^4) error. This method assumes equally spaced data
        //Simpson's Composite Quadrature Method
        if (y.length < 3) {
            throw new IllegalArgumentException("There are not enough data points to calculate the integral!");
        }
        double area = 0.0;
        double h = (b - a) / (y.length - 1);
         area += -y[0] + y[y.length - 1]; //Correcting before the loop
         for (int i = 1; i < y.length; i += 2) {
             area += 4 * y[i];
             area += 2 * y[i - 1];
         }
         return (h/3)*area;
    }

    public static double integrate(Function function, double a, double b) { //Returns the integral with O(h^4) error.
        //Simpson's Composite Quadrature Method
        double area = 0.0;
        area += -function.eval(a) + function.eval(b); //Correcting before the loop
        for (double i = a + hDEF; i < b; i += 2 * hDEF) {
            area += 4 * function.eval(i);
            area += 2 * function.eval(i-hDEF);
        }
        return (hDEF/3)*area;
    }
}
