package etomica;

import java.beans.BeanInfo;
import java.beans.IntrospectionException;
import java.beans.Introspector;

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
        
        InnerA innerA = new InnerA();
        InnerA innerAB = new InnerB();
        InnerB innerB = new InnerB();
        
        System.out.println(innerA.name());
        System.out.println(innerAB.name());
        System.out.println(innerB.name());
        System.out.println(((InnerA)innerB).name());
        
        //Introspection to get array of all properties
        java.beans.PropertyDescriptor[] properties = null;
        BeanInfo bi = null;
        try {
            bi = Introspector.getBeanInfo(innerB.getClass());
            properties = bi.getPropertyDescriptors();
        } 
        catch (IntrospectionException ex) {
            ex.printStackTrace();
        }
        for(int i=0; i<properties.length; i++) {
            System.out.println(properties[i].getName());
        }
    }
    
    private static class InnerA {
        public String name() {return "InnerA";}
        public void setProperty1(double prop1) {}
        public double getProperty1() {return 0.0;}
        public void setMyProperty2(int prop2) {}
        public int getMyProperty2() {return 0;}
        public void setOn(boolean prop3) {}
        public boolean isOn() {return true;}
    }
    
    private static class InnerB extends InnerA {
        public String name() {return "InnerB";}
        public void setPropB(int b) {}
        public int getPropB() {return 0;}
    }
}
        