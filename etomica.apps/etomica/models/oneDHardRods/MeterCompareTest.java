package etomica.models.oneDHardRods;

import etomica.api.IBox;
import etomica.api.IPotential;
import etomica.api.IPotentialMaster;
import etomica.api.IVectorMutable;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.units.Null;


/**
 * Ugly test code for hard rods in 1D.  making sure that our code is doing
 * what it is supposed to do.
 * 
 * 
 * @author cribbin
 *
 */
public class MeterCompareTest extends DataSourceScalar {
    int numTrials, numAccept;
    IPotential potentialTarget, potentialHarmonic;
    MeterPotentialEnergy meterPE;
    
    private double eigenVectors[][][];
    private IVectorMutable[] waveVectors;
    int comparedWV;
    protected double temperature;
    private double[] waveVectorCoefficients;
    private double wvc;
    private CoordinateDefinition coordinateDefinition;
    private double[] realT, imagT;
    private double[][] uOld, omegaSquared;
    private double[] uNow, deltaU;
    int coordinateDim;
    private double energyHardRod, energyHarmonic, energyOld;
    public double[] energyHardRodArray, energyHarmonicArray, energyOldArray;
    
    private static final long serialVersionUID = 1L;
    
    public MeterCompareTest(IPotentialMaster potentialMaster, CoordinateDefinition cd, IBox box){
        this("meterCompareTEST", potentialMaster, cd, box);
    }
    
    public MeterCompareTest(String string, IPotentialMaster potentialMaster, CoordinateDefinition cd, IBox box){
        super(string, Null.DIMENSION);
        setCoordinateDefinition(cd);
        realT = new double[coordinateDim];
        imagT = new double[coordinateDim];
        deltaU = new double[coordinateDim];
        meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
    }
        
    
    public double getDataAsScalar() {
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        BasisCell cell = cells[0];
        uOld = new double[cells.length][coordinateDim];
        double normalization = 1/Math.sqrt(cells.length);
        
        int wvlength = waveVectors.length;
        
        energyHardRodArray = new double[wvlength];
        energyHarmonicArray = new double[wvlength];
        energyOldArray = new double[wvlength];
        
        
        
        
        for(int countWV = wvlength-1; countWV > -1; countWV--){
            setComparedWV(countWV);
        
            energyHardRod = 0.0;
            energyHarmonic = 0.0;
            energyOld = meterPE.getDataAsScalar();
            energyOldArray[countWV] = energyOld;
            
            //get normal mode coordinate of "last" waveVector
            coordinateDefinition.calcT(waveVectors[comparedWV], realT, imagT);
            double realCoord = 0.0, imagCoord = 0.0;
            for(int i = 0; i < coordinateDim; i++){  //Loop would go away
                for(int j = 0; j < coordinateDim; j++){
                    realCoord += eigenVectors[comparedWV][i][j] * realT[j];
                    imagCoord += eigenVectors[comparedWV][i][j] * imagT[j];
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
                double kR = waveVectors[comparedWV].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                for(int i = 0; i < coordinateDim; i++){  //Loop would go away
                    //Calculate the current coordinates.
                    for(int j = 0; j < coordinateDim; j++){
                        deltaU[j] -= wvc*eigenVectors[comparedWV][i][j] *
                            2.0 * (realCoord*coskR - imagCoord*sinkR);
                    }
                }
    
                for(int i = 0; i < coordinateDim; i++){
                    deltaU[i] *= normalization;
                }
                
                for(int i = 0; i < coordinateDim; i++) {
                    uNow[i] += deltaU[i];
                }
                coordinateDefinition.setToU(cells[iCell].molecules, uNow);
            }
            energyHardRod = meterPE.getDataAsScalar();
            energyHardRodArray[countWV] += energyHardRod;
            
            //Calculate the energy due to the Gaussian modes.
            for(int i = 0; i < coordinateDim; i++){  //Loop would go away
                if(Double.isInfinite(omegaSquared[comparedWV][i])){
                    continue;
                }
                double normalCoord = realCoord*realCoord + imagCoord * imagCoord;
                energyHarmonic += wvc * normalCoord * omegaSquared[comparedWV][i];
            }
        energyHarmonicArray[countWV] += energyHarmonic;
        }
        
        // Set all the atoms back to the old values of u
        for (int iCell = 0; iCell<cells.length; iCell++) {
            cell = cells[iCell];
            coordinateDefinition.setToU(cell.molecules, uOld[iCell]);
        }
        
        System.out.println("HardRodenergies Harmonicenergies OldEnergies ");
        for(int i = 0; i < 17; i++){
            System.out.println(energyHardRodArray[i]+ " "+energyHarmonicArray[i]+ " "+energyOldArray[i]);
        }
        
        
//        if(getDataInfo().getLabel() == "meterBinA" && energyHardRod != 0.0 ){
//            System.out.println("energyOld  " + energyOld);
//            System.out.println("energyNM  " + energyHardRod);
//            System.out.println("energyOP  " + energyHarmonic);
//        }
        
        return energyHardRod + energyHarmonic;
    }

    
    public void setEigenVectors(double[][][] eigenVectors) {
        this.eigenVectors = eigenVectors;
    }
    public void setWaveVectors(IVectorMutable[] waveVectors) {
        this.waveVectors = waveVectors;
    }
    public void setComparedWV(int comparedWV) {
        this.comparedWV = comparedWV;
        wvc = waveVectorCoefficients[comparedWV];
    }
    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }
    public void setWaveVectorCoefficients(double[] waveVectorCoefficients) {
        this.waveVectorCoefficients = waveVectorCoefficients;
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
