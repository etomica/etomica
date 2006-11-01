package etomica.models.hexane;

import etomica.eigenstuff.MyEigenvalueDecomposition;

public class NormalModeEigenGetter {
    
    public static void main (String[] args){
        NormalModeReader reader = new NormalModeReader();
        String filename = "normal_modes";
        if (args.length > 0) {
            filename = args[0];
        }
        reader.readFromFile(filename);
        double[][] q = reader.q;
        double[][][] S = reader.S;
        for (int i=0; i<q.length; i++) {
            MyEigenvalueDecomposition evaler = new MyEigenvalueDecomposition(q[0].length, S[i]);
            double[] d = evaler.getRealEigenvalues();
            double[] e = evaler.getImagEigenvalues();
            double[][] eigenvectors = evaler.getV();

            System.out.print("q ");
            for (int k=0; k<q[i].length; k++) {
                System.out.print(q[i][k]+" ");
            }
            System.out.println("");

            System.out.println("eigenvectors");
            for (int j=0; j<eigenvectors.length; j++) {
                for (int k=0; k<eigenvectors[j].length; k++) {
                    System.out.print(eigenvectors[j][k]+" ");
                }
                System.out.println("");
            }
            System.out.println("real eigenvalues:");
            for (int j=0; j<d.length; j++) {
                System.out.print(d[j]+" ");
            }
            System.out.println("");
            
//            System.out.println("imaginary eigenvalues:");
//            for (int j=0; j<d.length; j++) {
//                System.out.print(e[j]+" ");
//            }
//            System.out.println("");
        }
    }
    
}
