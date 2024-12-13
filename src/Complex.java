public class Complex {

    public static ComplexNumber add(ComplexNumber z0, ComplexNumber z1) {
        ComplexNumber z = new ComplexNumber();
        z.setReal(z0.getReal() + z1.getReal());
        z.setImaginary(z0.getImaginary() + z1.getImaginary());
        return z;
    }

    public static ComplexNumber sub(ComplexNumber z0, ComplexNumber z1) { //z0 - z1
        ComplexNumber z = new ComplexNumber();
        z.setReal(z0.getReal() - z1.getReal());
        z.setImaginary(z0.getImaginary() - z1.getImaginary());
        return z;
    }

    public static ComplexNumber mul(ComplexNumber z0, ComplexNumber z1) {
        ComplexNumber z = new ComplexNumber();
        z.setReal(z0.getReal() * z1.getReal() - z0.getImaginary() * z1.getImaginary());
        z.setImaginary(z0.getReal() * z1.getImaginary() + z0.getImaginary() * z1.getReal());
        return z;
    }

    public static ComplexNumber div(ComplexNumber z0, ComplexNumber z1) { //z0/z1
        ComplexNumber z;
        ComplexNumber z1Con = conjugate(z1);
        ComplexNumber denominator =  mul(z1, z1Con);
        z = mul(z0, z1Con);
        z.setReal(z.getReal() / denominator.getReal());
        z.setImaginary(z.getImaginary() / denominator.getReal());
        return z;
    }

    public static ComplexNumber conjugate(ComplexNumber z0) {
        ComplexNumber z = new ComplexNumber();
        z.setReal(z0.getReal());
        z.setImaginary(-1*z0.getImaginary());
        return z;
    }

}
