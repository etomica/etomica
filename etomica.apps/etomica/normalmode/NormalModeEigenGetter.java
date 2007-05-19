package etomica.normalmode;

import java.io.FileWriter;
import java.io.IOException;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

public class NormalModeEigenGetter {
    
    public static void main (String[] args){
        String filename = "normal_modes_LJ_3D_32";
        if (args.length > 0) {
            filename = args[0];
        }
        doit(filename);
    }
    
    public static void doit(String filename) {
        //first index of S indicates the wave vector
        //the remaining two indices describe a square matrix of dimension equal to coordinateDim
        //e.g., coordinateDim = D for monatomic spherical molecules in D dimensions
        double[][][] S = ArrayReader2D.getFromFile(filename+".S");
        double maxIgnoredVal = 0;
        try {
            FileWriter fileWriterVal = new FileWriter(filename+".val");
            FileWriter fileWriterVec = new FileWriter(filename+".vec");
            for (int i=0; i<S.length; i++) {
                EigenvalueDecomposition evaler = new EigenvalueDecomposition(new Matrix(S[i]));
                double[] d = evaler.getRealEigenvalues();
    //            double[] e = evaler.getImagEigenvalues();
                Matrix eigenvectors = evaler.getV();
    
                for (int k=0; k<eigenvectors.getColumnDimension(); k++) {
                    for (int j=0; j<eigenvectors.getRowDimension(); j++) {
                        fileWriterVec.write(eigenvectors.get(j,k)+" ");
                    }
                    fileWriterVec.write("\n");
                }
                for (int j=0; j<d.length; j++) {
                    if (i==0 && j<3) {
                        if (Math.abs(d[j]) > maxIgnoredVal) {
                            maxIgnoredVal = d[j];
                        }
                        fileWriterVal.write("0 ");
                    }
                    else {
                        if (maxIgnoredVal > d[j] * 1E-10) {
                            System.err.println("Ignored "+maxIgnoredVal+" and now I found "+d[j]);
                        }
                        fileWriterVal.write(d[j]+" ");
                    }
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
