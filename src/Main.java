//WORK IN PROGRESS
//WORK IN PROGRESS
//WORK IN PROGRESS
//WORK IN PROGRESS
//WORK IN PROGRESS

//This class is used for informal testing. Once the entire project is complete, the UI will be updated and implemented here.
public class Main {

    public static void main(String[] args) {

        //y1' = y2
        //y2' = 2y2 - y1

        TriFunction<Double, Double, Double, Double> fnc1 = (t, y1, y2) -> y2;

        TriFunction<Double, Double, Double, Double> fnc2 = (t, y1, y2) -> 2 * y2 - y1;

        double[] y0 = {1, 1};

        double[][] v = ODE.rk4System(fnc1, fnc2, 0, y0, 0.2, 0.1);

        Matrix vals = new Matrix(v);

        System.out.println(vals);
    }
}
