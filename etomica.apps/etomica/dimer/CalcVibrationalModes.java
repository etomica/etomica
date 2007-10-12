package etomica.dimer;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

public class CalcVibrationalModes {

    double[][] forceConstant;
    double [] modes;
    double [] frequencies;
    int[] modeSigns;
    Matrix fC;
    EigenvalueDecomposition eigenDecomp;
    
    public CalcVibrationalModes(double [][] aForceConstantArray){
        
        modeSigns = new int[3];
        
        this.forceConstant = aForceConstantArray;
        fC = new Matrix(aForceConstantArray);
        
        // Finds Eigenvalues of Matrix fC
        eigenDecomp = new EigenvalueDecomposition(fC);
        
        }
    
    public double[] getLambdas(){
        return eigenDecomp.getRealEigenvalues();
    }
    
    public int[] getModeSigns(){
        
        modes = getLambdas();
        
        for(int i=0; i<modes.length; i++){
            
            if(modes[i]>0.0){modeSigns[0]++;}
            else {modeSigns[1]++;}
        }
        
        modeSigns[2] = modes.length;
        
        return modeSigns;
    }
    
    public double[] getFrequencies(){
        // where lambda = 4 * pi^2 * omega^2;
        modes = getLambdas();
        
        frequencies = new double[modes.length];
        
        for(int i=0; i<frequencies.length; i++){
            
            //Negative mode catch
            if(modes[i]<0.0){
                frequencies[i] = 0.0;
                continue;
            }
            
            frequencies[i] = Math.sqrt(modes[i]) / (2*Math.PI);
        }
        
        return frequencies;
    }
}
