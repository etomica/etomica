package etomica;

import etomica.utility.Arrays;

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
        
        Object[] object1 = new Object[2];
        Object[] object2 = new Object[2];
        object1[0] = new Object();
        object1[1] = new Object();
        object2[0] = object1[0];
        object2[1] = object1[1];
 //       System.out.println(object2.equals(object1));
 //       System.out.println(java.util.Arrays.equals(object1,object2));
        object2 = (Object[])object1.clone();
        System.out.println(object2.equals(object1));
        System.out.println(object2[0].equals(object1[0]));
        System.out.println(object2[1].equals(object1[1]));
        
        System.out.println(Arrays.toString(new int[] {}));
        System.out.println(Arrays.toString(new int[] {3}));
        System.out.println(Arrays.toString(new int[] {3,-1,5,20}));
        System.out.println(Arrays.toString((int[])null));
        System.out.println(Arrays.toString(new double[] {}));
        System.out.println(Arrays.toString(new double[] {3.0}));
        System.out.println(Arrays.toString(new double[] {3.2,4./3.,-5,Double.MAX_VALUE}));
        
    }
}
        