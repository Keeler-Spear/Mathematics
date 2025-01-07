//WORK IN PROGRESS
//WORK IN PROGRESS
//WORK IN PROGRESS
//WORK IN PROGRESS
//WORK IN PROGRESS
public class ComplexNumber {
    private double real;
    private double imaginary;
    private double modulus;
    private double theta;
    private static final double h = 0.00000000001;

    public ComplexNumber(double real, double imaginary) {
        this.real = real;
        this.imaginary = imaginary;
        updateModulus();
        updateTheta();
    }

    public ComplexNumber() {
        real = 0;
        imaginary = 0;
        updateModulus();
        updateTheta();
    }

    public double getReal() {
        return real;
    }

    public double getImaginary() {
        return imaginary;
    }

    public double getModulus() {
        return modulus;
    }

    public double getTheta () {
        return theta;
    }

    public void setReal(double real) {
        this.real = real;
        updateModulus();
        updateTheta();
    }

    public void setImaginary(double imaginary) {
        this.imaginary = imaginary;
        updateModulus();
        updateTheta();
    }

    private void updateModulus () {
        modulus = Math.sqrt(Math.pow(real, 2) + Math.pow(imaginary, 2));
    }

    private void updateTheta() {
        theta = Math.atan2(imaginary, real);
    }

    @Override
    public String toString() {
        String s = String.format("%.4f", real);
        if (imaginary < 0) {
            s += " - ";
        }
        else {
            s += " + ";
        }
        s += String.format("%.4f", Math.abs(imaginary));
        return s + "i";
    }

    @Override
    protected Object clone() throws CloneNotSupportedException {
        return new ComplexNumber(real, imaginary);
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }

        if (obj.getClass() != this.getClass()) {
            return false;
        }

        ComplexNumber z = (ComplexNumber) obj;
        if ((Math.abs(z.real - this.real) > h) || (Math.abs(z.imaginary - this.imaginary) > h)) {
            return false;
        }
        return true;
    }

}
