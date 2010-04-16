package etomica.models.oneDHardRods;

import etomica.math.SpecialFunctions;

public class MakeData {
    
    double n1, n2;
    
    public MakeData(int n, double density){
        n1 = 0.0;
        n2 = 0.0;
        for(int i = 1; i < n; i++){
            n1=n2;
            n2 = -(i)*Math.log((i+1)/density-(i+1)) + SpecialFunctions.lnFactorial(i) ;
             System.out.println(i + " delta " + (n2-n1));
//            System.out.println(i + " " + n2);
        }
    }
    
    public static void main(String[] args) {
        int rods = 101;
        double density = 0.5;
        MakeData sim = new MakeData(rods, density);
    }
}
