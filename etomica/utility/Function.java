package etomica.utility;

/**
 * Interface for the basic features of a function that maps a double onto another double.
 */

public interface Function {
    
    public double f(double x);
    public double inverse(double f);
    public double dfdx(double x);
 
    /**
     * The function f(x) = x
     */
    public static class Identity implements Function {
        public double f(double x) {return x;}
        public double inverse(double f) {return f;}
        public double dfdx(double x) {return 1.0;}
    }
    
    /**
     * The function f(x) = 1/x
     */
    public static class Reciprocal implements Function {
            
        public double f(double x){
                return 1.0/x;
            }
            
        public double dfdx(double x){
                return -1.0/(x*x);
            }
            
        public double inverse(double x){
                return 1.0/x;
            }
    }
    
    /**
     * The function f(x) = ln(x)
     */
    public static class Log implements Function{
            
        public double f(double x){
                return Math.log(x);
            }
            
        public double dfdx(double x){
                return 1.0/x;
            }
            
        public double inverse(double x){
                return Math.exp(x);
            }
         
        }

    /**
     * The function f(x) = a*x + b
     */
    public static class Linear implements Function {
        private final double a, b, ra;
        public Linear(double slope, double intercept) {
            this.a = slope;
            this.b = intercept;
            ra = 1.0/a;
        }
        public double f(double x) {return a*x + b;}
        public double inverse(double f) {return ra*(f-b);}
        public double dfdx(double x) {return a;}
    }
    

}