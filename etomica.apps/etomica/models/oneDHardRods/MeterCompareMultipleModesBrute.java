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

public class MeterCompareMultipleModesBrute extends DataSourceScalar {
    int numTrials, numAccept;
    IPotential potentialTarget, potentialHarmonic;
    MeterPotentialEnergy meterPE;
    
    private double eigenVectors[][][];
    private IVectorMutable[] waveVectors;
    int[] comparedWVs;
    protected double temperature;
    private double[] waveVectorCoefficients;
    private CoordinateDefinition coordinateDefinition;
    private double[] realT, imagT;
    private double[][] uOld, omegaSquared;
    private double[] uNow, deltaU;
    int coordinateDim;
    private double energyHardRod, energyHarmonic;
    
    private static final long serialVersionUID = 1L;
    
    public MeterCompareSingleModeBrute single;
    
    public boolean isA;
    
    
    public MeterCompareMultipleModesBrute(IPotentialMaster potentialMaster, 
            CoordinateDefinition cd, IBox box){
        this("meterCompareMultipleModes", potentialMaster, cd, box);
    }
    
    public MeterCompareMultipleModesBrute(String string, IPotentialMaster
            potentialMaster, CoordinateDefinition cd, IBox box){
        super(string, Null.DIMENSION);
        setCoordinateDefinition(cd);
        realT = new double[coordinateDim];
        imagT = new double[coordinateDim];
        deltaU = new double[coordinateDim];
        meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
        
        single = new MeterCompareSingleModeBrute("single", potentialMaster, cd, box);
    }
    
    
    
    public double getDataAsScalar(){
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        BasisCell cell = cells[0];
        uOld = new double[cells.length][coordinateDim];
        double normalization = 1/Math.sqrt(cells.length);
        int numWV = comparedWVs.length;
        double[] realCoord = new double[numWV];
        double[] imagCoord = new double[numWV];
        energyHardRod = 0.0;
        energyHarmonic = 0.0;
        
        

        System.out.println("single: " + single.getDataAsScalar());
        
        
        //Get the normal mode coordinates of the compared waveVectors, and
        // store them in realCoord and imagCoord for further use.
        for(int wvcount = 0; wvcount < numWV; wvcount++){
            coordinateDefinition.calcT(waveVectors[wvcount], realT, imagT);
//            System.out.println("real " +realT[0]);
//            System.out.println("imag " +imagT[0]);

            
            imagCoord[wvcount] = 0.0;
            for(int i = 0; i < coordinateDim; i++){  //Loop would go away
                for(int j = 0; j < coordinateDim; j++){
                    realCoord[wvcount] += eigenVectors[wvcount][i][j] * realT[j];
                    imagCoord[wvcount] += eigenVectors[wvcount][i][j] * imagT[j];
                }
            }
        }
        
        //Store the original positions
        for(int i = 0; i < cells.length; i++){
            uNow = coordinateDefinition.calcU(cells[i].molecules);
            System.arraycopy(uNow, 0, uOld[i], 0, coordinateDim);
        }
        
        //Remove the effects of the compared modes.
        for(int wvcount = 0; wvcount < numWV; wvcount++){
            for(int iCell = 0; iCell < cells.length; iCell++){
                cell = cells[iCell];
                for(int j = 0; j < coordinateDim; j++){
                    deltaU[j] = 0.0;
                }
                
                //Calculate the contributions to the current position of the 
                //zeroed mode, and subtract it from the overall position.
                double kR = waveVectors[wvcount].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                for(int i = 0; i < coordinateDim; i++){  //Loop would go away
                    //Calculate the current coordinates.
                    for(int j = 0; j < coordinateDim; j++){
                        deltaU[j] -= waveVectorCoefficients[wvcount] * 
                            eigenVectors[wvcount][i][j] * 2.0 * (realCoord[wvcount] *
                            coskR - imagCoord[wvcount]*sinkR);
                    }
                }

                for(int i = 0; i < coordinateDim; i++){
                    deltaU[i] *= normalization;
                }
                
                for(int i = 0; i < coordinateDim; i++) {
                    uNow[i] += deltaU[i];
                }
                coordinateDefinition.setToU(cells[iCell].molecules, uNow);
            }//end of cell loop
        }//end of wvcount loop
        energyHardRod = meterPE.getDataAsScalar();
        
        //Calculate the energy due to the compared modes
        for(int wvcount = 0; wvcount < numWV; wvcount++){
            for(int i = 0; i < coordinateDim; i++){  //Loop would go away
                if(Double.isInfinite(omegaSquared[wvcount][i])){
                    continue;
                }
                double normalCoord = realCoord[wvcount]*realCoord[wvcount] + 
                    imagCoord[wvcount] * imagCoord[wvcount];
                energyHarmonic += waveVectorCoefficients[wvcount] * normalCoord 
                    * omegaSquared[wvcount][i];
            }
        }
        
     // Set all the atoms back to the old values of u
        for (int iCell = 0; iCell<cells.length; iCell++) {
            cell = cells[iCell];
            coordinateDefinition.setToU(cell.molecules, uOld[iCell]);
        }
        
//        if(isA){
//            double dork = energyHardRod + energyHarmonic;
//            System.out.println("single: " + single.getDataAsScalar());
//            System.out.println("HR: " + energyHardRod);
//            System.out.println("Harm: " + energyHarmonic);
//        }
        
        return energyHardRod + energyHarmonic;
    }
    
    
    public void setEigenVectors(double[][][] eigenVectors) {
        this.eigenVectors = eigenVectors;
    }
    public void setWaveVectors(IVectorMutable[] waveVectors) {
        this.waveVectors = waveVectors;
    }
    public void setComparedWVs(int[] cwvs) {
        this.comparedWVs = cwvs;
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

    public MeterCompareSingleModeBrute getSingle() {
        return single;
    }

    public void setSingle(MeterCompareSingleModeBrute single) {
        this.single = single;
    }

    public boolean isA() {
        return isA;
    }

    public void setA(boolean isA) {
        this.isA = isA;
    }
}
