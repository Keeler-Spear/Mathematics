# Project Summary

  This project is a coding project I worked on in the latter half of 2024. It can perform the fundamental operations of linear algebra and single-variable calculus via analytical and/or numerical methods. All code is in Java, as it is the language I feel the most comfortable with.
  
# File summary

-Mathematics.BasisFunctions.java: Static class that generates basis functions.\
-Mathematics.Calculus.java: Static class with calculus operations.\
-Mathematics.Complex.java: Static class with complex arithmetic, WIP.\
-Mathematics.ComplexNumber.java: Non-static class with the properties of a complex number, WIP.\
-Mathematics.Error.java: Static class with error evaluation operations.\
-LinReg.Java: Static class with numerical linear regression operations.\
-Mathematics.LinearAlgebra.java: Static class with linear algebra operations.\
-Main.java: Main class for user input, WIP.\
-Mathematics.Maths.java: Static class that can perform some of the operations found in the Java Math library, WIP.\
-Mathematics.Matrix.java: Non-static class with the properties of a mathematical matrix and one’s based indexing.\
-Mathematics.ODE.java: Static class with numerical ordinary differential equation operations, WIP.\
-Mathematics.TriFunction.java: Interface that acts as a mathematical function with three independent variables.\
-UI.java: Main class that provides a user interface for the most important linear algebra methods.

# Main Operations and Method Type

-Mathematics.Matrix Addition (Analytically)\
-Mathematics.Matrix Multiplication (Analytically)\
-Mathematics.Matrix Row Reduction (Analytically)\
-LU Decomposition (Analytically)\
-Mathematics.Matrix Inversion (Analytically)\
-Determinant Calculation (Analytically)\
-Gram-Schmidt Orthogonalization (Analytically)\
-Eigenvalue Calculation (Numerically via the QR Method)\
-Differentiation (Numerically via O(h^2) Methods)\
-Integration (Numerically via O(h^4) Methods)\
-Mathematics.Regression (Numerically via gradient descent)

# User Guide

  There are two ways to use the code. First, if you plan to perform a common linear algebra operation, use UI.java for a smooth experience. Second, if the operation is not offered in UI.java or you prefer not to use UI.java, use the public methods provided to execute your operation in Main.java. To use any of the calculus methods with a function (as opposed to a data set), you must use java.util.function.Function interface.
  
# Motivations

  My desire to create this project was driven by three factors. First, I did so to deepen my knowledge of numerical analysis and to apply the theory I learned. In parallel with my work on this project, I took Numerical Analysis I, in which I learned to derive numerical methods and evaluate their accuracy/error, efficiency/convergence, and stability; however, it did not heavily emphasize the coding implementation of the methods taught. Therefore, I explored this aspect of numerical analysis on my own time by creating this project. Second, I made this project as a tool. Whenever I mess around with circuits and attempt to analyze them, linear algebra operations consume much of my time. By using my code, I can analyze circuits instantly (once they have been translated into a math problem). While other software exists, such as MATLAB, it is often expensive, cumbersome, and slow (at least on my machine), therefore I created my own tool. Third, I made this project for fun. I find engineering and design enjoyable, so making this project provided me the opportunity to learn while having a blast. For these three reasons, I made this project.
  
# Constraints

  The only notable constraint for this project was time. As I used it as an educational opportunity, my work corresponded to my learning of certain numerical methods in class, resulting in a work period of roughly one week per method(s). Due to this, much of my code is unoptimized (while it is insanely fast for the problems I usually solve, it is inefficient for larger problems), and not every method I wanted to code was implemented, such as stationary matrix methods.
  
# Contributions

  Outside of the Java libraries used (Math, Swing, JUnit Jupiter, IO, and Util) and the MeersenneTwister class, the code is entirely my own. The theory and formulas behind the methods implemented were taught to me by my Numerical Analysis I professor and my Linear Algebra professor. However, Mathematics.Calculus.java’s leftSidedThreePoint (double f0, double f1, double f2, double h) and rightSidedThreePoint (double f0, double fm1, double fm2, double h) formulas were derived on my own.

  # Favorite Algorithm

  Of the algorithms I have coded, gradient descent (linearReg in Mathematics.Regression.java) is my favorite (made in mid-February). It can calculate the least squares fit of any data set, including those of multiple variables, and model it with the functions provided. For simplicity, I have created a couple of function sets under Mathematics.BasisFunctions.
