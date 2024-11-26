//Todo

public class Function {

    private String function;

    public Function (String function) {
        this.function = function;
    }

    public String getFunction() {
        return function;
    }

    public void setFunction(String function) {
        this.function = function;
    }

    public double evaluate(double x) {
        return 2*x*x;
    }

}
