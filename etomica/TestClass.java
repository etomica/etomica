package etomica;

public class TestClass {
    
    public static void main(String arg[]) {
        
        double x0 = -Double.MAX_VALUE;
        double x1 = +Double.MAX_VALUE;
        double x2 = -x1;
        
        System.out.println("x0 = -Double.MAX_VALUE: "+x0);
        System.out.println("x1 = +Double.MAX_VALUE: "+x1);
        System.out.println("x2 = -x1: "+x2);
        System.out.println("x0 <= -Double.MAX_VALUE: "+(x0 <= -Double.MAX_VALUE));
        System.out.println("x1 <= -Double.MAX_VALUE: "+(x1 <= -Double.MAX_VALUE));
        System.out.println("x2 <= -Double.MAX_VALUE: "+(x2 <= -Double.MAX_VALUE));
        
    }
}
        