//WORK IN PROGRESS
//WORK IN PROGRESS
//WORK IN PROGRESS
//WORK IN PROGRESS
//WORK IN PROGRESS

//ToDo: Implement a way for a String to be converted into a mathematical function.
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

    //For now, just type your function here to use the Calculus.java operations.
    public double eval(double x) {
        return Math.exp(x*x);
    }

}
