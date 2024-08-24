import java.util.Scanner;

public class rungeKuttaMethod {

    //Eulers Improved formula :Y(n+1)=Y(n)+h/2 (f(xn,yn)+f(xn+1,y*n+1)
    // y*n+1=yn+hf(xn,yn)

    //function
    public static double function(double x, double y) {
        y = x + y * y;
        // Change to your function

        return y;

    }

    public static double k1(double x, double y, double h) {
        y = function(x,y);
        return y;

    }

    public static double k2(double x, double y, double h) {
        y = function(x + 0.5 * h,y + 0.5 * h * k1(x, y, h));
        return y;
    }

    public static double k3(double x, double y, double h) {
        y = function(x + 0.5 * h,y + 0.5 * h * k2(x, y, h));
        return y;

    }

    public static double k4(double x, double y, double h) {
        y = function(x + h,y + h * k3(x, y, h));
        return y;

    }


    public static void main (String[] args) {
        Scanner scnr = new Scanner(System.in);
        //Declared Variable
        int i;
        double x=0.0;
        double y=0.0;
        int n;
        double h;
        double ys;

        //Asks for variables
        System.out.println("How many iterations do you want to do?");
        n = scnr.nextInt();
        System.out.println("What is h?");
        h = scnr.nextDouble();
        System.out.println("What is x?");
        x = scnr.nextDouble();
        System.out.println("What is y?");
        y = scnr.nextDouble();
        System.out.println();

        for (i = 0; i < n; ++i) {
            System.out.print("F(");
            System.out.printf("%.4f", (x + h));
            System.out.print(") = ");
            y = y + (h/6.0) * (k1(x, y, h) + 2 * k2(x, y, h) + 2 * k3(x, y, h) + k4(x, y, h));
            System.out.printf("%.4f", y);
            System.out.println();
            x = x + h;
        }






    }

}
