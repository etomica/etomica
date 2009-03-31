package etomica.normalmode;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.api.IBox;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.space.ISpace;

/**
 *
 *
 * @author Tai Boon Tan
 */
public class NormalModes3D implements NormalModes {

    /**
     */
    public NormalModes3D(ISpace _space, Primitive primitive, Basis basis) {
    	
    	this.space = _space;
    	this.primitive = primitive;
    	this.basis = basis;
    	waveVectorFactory = new WaveVectorFactorySimple(primitive, space);
        harmonicFudge = 1;
        eigenvalues = ArrayReader1D.getFromFile("DB_FCC_n12_N5.val");
        eigenvectors = ArrayReader2D.getFromFile("DB_FCC_n12_N5.vec");
    }

    public double[][] getOmegaSquared(IBox box) {
   
    	double[][] omega2 = new double[eigenvalues.length][eigenvalues[0].length];
    	for (int i=0; i<omega2.length; i++) {
    		for (int j=0; j<omega2[i].length; j++) {
    			omega2[i][j] = 1.0/eigenvalues[i][j];
               }
    	}
    	
    	/*
    	 *  set all the eigenvalues to 0.001
    	 *  omega2 = 1/eigenvalue
    	 */
    	/*
        for (int i=0; i<omega2.length; i++) {
            for (int j=0; j<omega2[i].length; j++) {
            	if (i==0 && j<12){
            		omega2[i][j] = Double.POSITIVE_INFINITY;
            		
            	} else {
            		omega2[i][j] = 1000*temperature;
            		
            	}
            }
        }
        */
    	
    	
        return omega2;
    }

    public double[][][] getEigenvectors(IBox box) {
  
    	/*
    	for(int i=0; i<numWV; i++){
    		for(int j=0; j<numEigen; j++){
    			for(int k=0; k<numEigen; k++){
    				if(j==k){
    					eigenvectors[i][j][k] = 0.8;
    					
    				} else {
    					eigenvectors[i][j][k] = 0.0;
    					
    				}
    			}
    		}
    	}
    	*/
    	
        return eigenvectors;
    }
                      
                      
    public WaveVectorFactory getWaveVectorFactory() {
        return waveVectorFactory;
    }
    
    public void setHarmonicFudge(double newHarmonicFudge) {
        harmonicFudge = newHarmonicFudge;
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }


    
    public static void main(String[] args){
    	double[][] array = {{1,-1,0},{-1,2,-1},{0,-1,1}};
    	Matrix matrix = new Matrix(array);
    	for (int i=0; i<3; i++){
    		for (int j=0; j<3; j++){
    			System.out.print(matrix.get(i, j)+ " ");
    		}
    		System.out.println("");
    	}
    	
    	EigenvalueDecomposition ed = matrix.eig();
    	double[] eVal = ed.getRealEigenvalues();
    	double[][] eVec = ed.getV().getArray();
    	
    	System.out.print("\nEigenvalues: ");
    	for (int i=0; i<3; i++){

    		System.out.print(eVal[i]+ " ");
    		
    	}
    	
    	System.out.println("\n\nEigenvectors: ");
    	for (int i=0; i<3; i++){
    		for (int j=0; j<3; j++){
    			System.out.print(eVec[i][j]+ " ");
    		}
    		System.out.println("");
    	}
   		
   		
    }
    
    protected WaveVectorFactory waveVectorFactory;
    protected double harmonicFudge;
    protected double temperature;
    protected Primitive primitive;
    protected ISpace space;
    protected Basis basis;
    protected double[][][] eigenvectors;
    protected double[][] eigenvalues;

}
