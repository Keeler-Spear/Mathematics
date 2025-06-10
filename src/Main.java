import Mathematics.*;

import static Mathematics.LinearAlgebra.addMatrices;
import static Mathematics.LinearAlgebra.randMatrix;


//This class is used for informal testing. Once the entire project is complete, the UI will be updated and implemented here.
public class Main {

    public static void main(String[] args) {

        Matrix A = randMatrix(3, 3, -10.0, 10.0);
        Matrix B = randMatrix(3, 3, -10.0, 10.0);

        System.out.println(addMatrices(A, B));


    }
}
