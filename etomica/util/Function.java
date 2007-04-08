package etomica.util;

/**
 * Interface for the basic features of a function that maps a double onto
 * another double.
 */
public interface Function {

    public double f(double x);

    /**
     * The function f(x) = constant
     */
    public static class Constant implements FunctionInvertible, FunctionDifferentiable, java.io.Serializable {

        private double c;

        public Constant() {
            this(0.0);
        }

        public Constant(double c) {
            this.c = c;
        }

        public double f(double x) {
            return c;
        }

        public double inverse(double f) {
            throw new RuntimeException("function not invertible");
        }

        public double dfdx(double x) {
            return 0.0;
        }
    }

    /**
     * The function f(x) = x
     */
    public static class Identity implements FunctionInvertible, FunctionDifferentiable, java.io.Serializable {

        public double f(double x) {
            return x;
        }

        public double inverse(double f) {
            return f;
        }

        public double dfdx(double x) {
            return 1.0;
        }

        public final static Identity INSTANCE = new Identity();
    }

    /**
     * The function f(x) = 1/x
     */
    public static class Reciprocal implements FunctionInvertible, FunctionDifferentiable, java.io.Serializable {

        public double f(double x) {
            return 1.0 / x;
        }

        public double dfdx(double x) {
            return -1.0 / (x * x);
        }

        public double inverse(double x) {
            return 1.0 / x;
        }

        public final static Reciprocal INSTANCE = new Reciprocal();
    }

    /**
     * The function f(x) = exp(x)
     */
    public static class Exp implements FunctionInvertible, FunctionDifferentiable, java.io.Serializable {

        public double f(double x) {
            return Math.exp(x);
        }

        public double dfdx(double x) {
            return Math.exp(x);
        }

        public double inverse(double x) {
            return Math.log(x);
        }

        public final static Exp INSTANCE = new Exp();
    }

    /**
     * The function f(x) = ln(x)
     */
    public static class Log implements FunctionInvertible, FunctionDifferentiable, java.io.Serializable {

        public double f(double x) {
            return Math.log(x);
        }

        public double dfdx(double x) {
            return 1.0 / x;
        }

        public double inverse(double x) {
            return Math.exp(x);
        }

        public final static Log INSTANCE = new Log();
    }

    /**
     * The function f(x) = sqrt(x)
     */
    public static class Sqrt implements FunctionInvertible, FunctionDifferentiable, java.io.Serializable {

        public double f(double x) {
            return Math.sqrt(x);
        }

        public double dfdx(double x) {
            return 0.5 / Math.sqrt(x);
        }

        public double inverse(double x) {
            return x * x;
        }

        public final static Sqrt INSTANCE = new Sqrt();
    }

    /**
     * The function f(x) = abs(x)
     */
    public static class Abs implements FunctionInvertible, FunctionDifferentiable, java.io.Serializable {

        public double f(double x) {
            return Math.abs(x);
        }

        public double dfdx(double x) {
            return x > 0 ? 1 : -1;
        }

        public double inverse(double x) {
            throw new RuntimeException("function not invertible");
        }

        public final static Abs INSTANCE = new Abs();
    }

    /**
     * The function f(x) = a*x + b
     */
    public static class Linear implements FunctionInvertible, FunctionDifferentiable, java.io.Serializable {

        private final double a, b, ra;

        public Linear(double slope, double intercept) {
            this.a = slope;
            this.b = intercept;
            ra = 1.0 / a;
        }

        public double f(double x) {
            return a * x + b;
        }

        public double inverse(double f) {
            return ra * (f - b);
        }

        public double dfdx(double x) {
            return a;
        }
    }
    
    /**
     * Returns the input value if its magnitude is greater than value specified at
     * construction; otherwise returns zero.  Used to zero values that are nonzero
     * only because of roundoff errors.
     */
    public static class Chop implements Function, java.io.Serializable {
        
        private final double eps;
        
        /**
         * Constructs with default cutoff value of 1e-10.
         */
        public Chop() {
            this(1.e-10);
        }
        
        public Chop(double eps) {
            if(eps < 0) eps = -eps;
            this.eps = eps;
        }
        
        public double f(double x) {
            return (Math.abs(x) > eps) ? x : 0.0;
        }
    }

}
