package etomica.models.oneDHardRods;

import etomica.api.IPotential;
import etomica.api.IPotentialMaster;
import etomica.api.IVector;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.units.Null;


/**
 * Uses brute force to calculate energy of a system with a set number of 
 * normal modes, and the rest of the degrees of freedom taken by Gaussians.
 * Assumes 1D system - otherwise, choose a mode and eliminate i loops.
 * 
 * @author cribbin
 *
 */
public class MeterConvertModeBrute extends DataSourceScalar {
    int numTrials, numAccept;
    IPotential potentialTarget, potentialHarmonic;
    MeterPotentialEnergy meterPE;
    
    private double eigenvectors[][][];
    private IVector[] wavevectors;
    int convertedWV;
    protected double temperature;
    private double[] wavevectorCoefficients;
    private double wvc;
    private CoordinateDefinition coordinateDefinition;
    private double[] realT, imagT;
    private double[][] uOld, omegaSquared;
    private double[] uNow, deltaU;
    int coordinateDim;
    private double energyNM, energyOP;
    
    private static final long serialVersionUID = 1L;
    
    
    public MeterConvertModeBrute(IPotentialMaster potentialMaster){
        super("meterConvertMode", Null.DIMENSION);
        realT = new double[coordinateDim];
        imagT = new double[coordinateDim];
        deltaU = new double[coordinateDim];
        meterPE = new MeterPotentialEnergy(potentialMaster);
        
    }    
    
    public double getDataAsScalar() {
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        BasisCell cell = cells[0];
        uOld = new double[cells.length][coordinateDim];
        energyNM = 0.0;
        energyOP = 0.0;
        
        //get normal mode coordinate of "last" wavevector
        coordinateDefinition.calcT(wavevectors[convertedWV], realT, imagT);
        double realCoord = 0.0, imagCoord = 0.0;
        for(int i = 0; i < coordinateDim; i++){  //Loop would go away
            for(int j = 0; j < coordinateDim; j++){
                realCoord += eigenvectors[convertedWV][i][j] * realT[j];
                imagCoord += eigenvectors[convertedWV][i][j] * imagT[j];
            }
        }
        for(int iCell = 0; iCell < cells.length; iCell++){
            //store original positions
            uNow = coordinateDefinition.calcU(cells[iCell].molecules);
            System.arraycopy(uNow, 0, uOld[iCell], 0, coordinateDim);
            cell = cells[iCell];
            for(int j = 0; j < coordinateDim; j++){
                deltaU[j] = 0.0;
            }
            
            //Calculate the contributions to the current position of the 
            //zeroed mode, and subtract it from the overall position.
            double kR = wavevectors[convertedWV].dot(cell.cellPosition);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);
            for(int i = 0; i < coordinateDim; i++){  //Loop would go away
                //Calculate the current coordinates.
                for(int j = 0; j < coordinateDim; j++){
                    deltaU[j] -= wvc*eigenvectors[convertedWV][i][j] *
                        2.0 * (realCoord*coskR - imagCoord*sinkR);
                }
            }
            
            for(int i = 0; i < coordinateDim; i++) {
                uNow[i] += deltaU[i];
            }
            coordinateDefinition.setToU(cells[iCell].molecules, uNow);
        }
        energyNM = meterPE.getDataAsScalar();
        
        //Calculate the energy due to the Gaussian modes.
        for(int i = 0; i < coordinateDim; i++){  //Loop would go away
            if(Double.isInfinite(omegaSquared[convertedWV][i])){
                continue;
            }
            double normalCoord = realCoord*realCoord + imagCoord * imagCoord;
            energyOP += wvc * normalCoord * omegaSquared[convertedWV][i];
        }
        
        // Set all the atoms back to the old values of u
        for (int iCell = 0; iCell<cells.length; iCell++) {
            cell = cells[iCell];
            coordinateDefinition.setToU(cell.molecules, uOld[iCell]);
        }
        
        return energyNM + energyOP;
    }

    
    public void setEigenvectors(double[][][] eigenvectors) {
        this.eigenvectors = eigenvectors;
    }
    public void setWavevectors(IVector[] wavevectors) {
        this.wavevectors = wavevectors;
    }
    public void setConvertedWV(int convertedWV) {
        this.convertedWV = convertedWV;
        wvc = wavevectorCoefficients[convertedWV];
    }
    public void setTermperature(double temperature) {
        this.temperature = temperature;
    }
    public void setWavevectorCoefficients(double[] wavevectorCoefficients) {
        this.wavevectorCoefficients = wavevectorCoefficients;
    }
    public void setWvc(double wvc) {
        this.wvc = wvc;
    }
    public void setCoordinateDefinition(CoordinateDefinition cd){
        coordinateDefinition = cd;
        coordinateDim = coordinateDefinition.getCoordinateDim();
    }
    
    public void setSpringConstants(double[][] sc){
        omegaSquared = sc;
    }
    
    public void setOmegaSquared(double[][] sc){
        omegaSquared = sc;
    }
    
    

}
