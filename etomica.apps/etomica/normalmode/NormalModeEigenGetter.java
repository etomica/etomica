package etomica.normalmode;

import java.io.FileWriter;
import java.io.IOException;

import etomica.eigenstuff.MyEigenvalueDecomposition;

public class NormalModeEigenGetter {
    
    public static void main (String[] args){
        String filename = "normal_modes4000";
        if (args.length > 0) {
            filename = args[0];
        }
        double[][][] S = ArrayReader2D.getFromFile(filename+".S");
        try {
            FileWriter fileWriterVal = new FileWriter(filename+".val");
            FileWriter fileWriterVec = new FileWriter(filename+".vec");
            for (int i=0; i<S.length; i++) {
                MyEigenvalueDecomposition evaler = new MyEigenvalueDecomposition(S[0].length, S[i]);
                double[] d = evaler.getRealEigenvalues();
    //            double[] e = evaler.getImagEigenvalues();
                double[][] eigenvectors = evaler.getV();
    
                for (int j=0; j<eigenvectors.length; j++) {
                    for (int k=0; k<eigenvectors[j].length; k++) {
                        fileWriterVec.write(eigenvectors[j][k]+" ");
                    }
                    fileWriterVec.write("\n");
                }
                for (int j=0; j<d.length; j++) {
                    fileWriterVal.write(d[j]+" ");
                }
                fileWriterVal.write("\n");
                
    //            System.out.println("imaginary eigenvalues:");
    //            for (int j=0; j<d.length; j++) {
    //                System.out.print(e[j]+" ");
    //            }
    //            System.out.println("");
            }
            fileWriterVal.close();
            fileWriterVec.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
    
}
