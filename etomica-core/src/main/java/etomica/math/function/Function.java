/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.function;

/**
 * Interface for the basic features of a function that maps a double onto
 * another double.
 */
public interface Function extends IFunction {

    /**
     * The function f(x) = constant
     */
    public static class Constant implements FunctionDifferentiable, java.io.Serializable {

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

        public double df(int n, double x) {
            if(n < 0) {
                throw new IllegalArgumentException("Order of derivative must be non-negative");
            }
            switch(n) {
            case 0:
                return f(x);
            default:
                return 0.0;
            }
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

        public double df(int n, double x) {
            if(n < 0) {
                throw new IllegalArgumentException("Order of derivative must be non-negative");
            }
            switch(n) {
            case 1:
                return 1.0;
            case 0:
                return f(x);
            default:
                return 0.0;
            }
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

        public double df(int n, double x) {
            if(n < 0) {
                throw new IllegalArgumentException("Order of derivative must be non-negative");
            }
            double prod = 1.0/x;
            for(int i=0; i<n; i++) {
                prod *= -(i+1)/x;
            }
            return prod;
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

        public double df(int n, double x) {
            if(n < 0) {
                throw new IllegalArgumentException("Order of derivative must be non-negative");
            }
            return f(x);
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

        public double df(int n, double x) {
            if(n < 0) {
                throw new IllegalArgumentException("Order of derivative must be non-negative");
            }
            switch(n) {
            case 1:
                return 1.0/x;
            case 0:
                return f(x);
            default:
                double prod = 1.0/x;
                for(int i=1; i<n; i++) {
                    prod *= -i/x;
                }
                return prod;
            }
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

        public double df(int n, double x) {
            if(n < 0) {
                throw new IllegalArgumentException("Order of derivative must be non-negative");
            }
            switch(n) {
            case 1:
                return 0.5/Math.sqrt(x);
            case 0:
                return f(x);
            default:
                double prod = 0.5/Math.sqrt(x);
                for(int i=1; i<n; i++) {
                    prod *= -(i-0.5)/x;
                }
                return prod;
            }
        }

        public double inverse(double x) {
            return x * x;
        }

        public final static Sqrt INSTANCE = new Sqrt();
    }

    /**
     * The function f(x) = abs(x)
     */
    public static class Abs implements FunctionDifferentiable, java.io.Serializable {

        public double f(double x) {
            return Math.abs(x);
        }

        public double df(int n, double x) {
            if(n < 0) {
                throw new IllegalArgumentException("Order of derivative must be non-negative");
            }
            switch(n) {
            case 1:
                return x > 0 ? 1 : -1;
            case 0:
                return f(x);
            default:
                return 0.0;
            }
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

        public double df(int n, double x) {
            if(n < 0) {
                throw new IllegalArgumentException("Order of derivative must be non-negative");
            }
            switch(n) {
            case 1:
                return a;
            case 0:
                return f(x);
            default:
                return 0.0;
            }
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
         * Constructs with default cutoff value of 2 Double.MIN_VALUE.
         */
        public Chop() {
            this(2.0*Double.MIN_VALUE);
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
