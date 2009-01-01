package etomica.normalmode;

import etomica.api.IBox;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.space.ISpace;

/**
 *
 * 
 */
public class NormalModes2D implements NormalModes {

    /**
     */
    public NormalModes2D(ISpace _space, Primitive primitive, Basis basis) {
    	
    	this.space = _space;
    	this.primitive = primitive;
    	this.basis = basis;
    	waveVectorFactory = new WaveVectorFactory2D(primitive, space);
        harmonicFudge = 1;
        
    }

    public double[][] getOmegaSquared(IBox box) {
    	
    	double[] d = primitive.getSize();
    	int[] numCells = new int[space.D()];
    	
    	for(int i=0; i<space.D(); i++){
    		numCells[i] = (int)Math.round(box.getBoundary().getDimensions().x(i)/(d[i]));
    		
    	}
    	
    	int numWV = waveVectorFactory.getWaveVectors().length;
    	int numEigen = basis.getScaledCoordinates().length*2;
    	
    	double[][] omega2 = new double[numWV][numEigen];
    	
    	/*
    	 *  set all the eigenvalues to 0.001
    	 *  omega2 = 1/eigenvalue
    	 */
    	
        for (int i=0; i<omega2.length; i++) {
            for (int j=0; j<omega2[i].length; j++) {
            	if (i==0 && j<1){
            		omega2[i][j] = Double.POSITIVE_INFINITY;
            		
            	} else {
            		omega2[i][j] = 1000*temperature;
            		
            	}
            }
        }
        return omega2;
    }

    
    
    public double[][][] getEigenvectors(IBox box) {
    	
    	double[] d = primitive.getSize();
    	int[] numCells = new int[space.D()];
    	
    	for(int i=0; i<space.D(); i++){
    		numCells[i] = (int)Math.round(box.getBoundary().getDimensions().x(i)/(d[i]));
    		
    	}
    	
    	int numWV = waveVectorFactory.getWaveVectors().length;
    	int numEigen = basis.getScaledCoordinates().length*2;
    	
    	eigenvectors = new double[numWV][numEigen][numEigen];
    	
    	for(int i=0; i<numWV; i++){
    		for(int j=0; j<numEigen; j++){
    			for(int k=0; k<numEigen; k++){
    				if(j==k){
    					eigenvectors[i][j][k] = 1.0;
    					
    				} else {
    					eigenvectors[i][j][k] = 0.0;
    					
    				}
    			}
    		}
    	}
    	
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


    protected double[][][] eigenvectors;
    protected WaveVectorFactory waveVectorFactory;
    protected double harmonicFudge;
    protected double temperature;
    protected Primitive primitive;
    protected ISpace space;
    protected Basis basis;
    

}
