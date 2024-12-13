public class ComplexNumber {
    private double real;
    private double imaginary;
    private static final double h = 0.00000000001;

    public ComplexNumber(double real, double imaginary) {
        this.real = real;
        this.imaginary = imaginary;
    }

    public ComplexNumber() {
        real = 0;
        imaginary = 0;
    }

    public double getReal() {
        return real;
    }

    public double getImaginary() {
        return imaginary;
    }

    public void setReal(double real) {
        this.real = real;
    }

    public void setImaginary(double imaginary) {
        this.imaginary = imaginary;
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
    protected Object clone() throws CloneNotSupportedException { //Copies the array and then makes a new matrix with the copy
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
