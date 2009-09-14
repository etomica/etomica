package etomica.models.oneDHardRods;

import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;

/**
 * A Monte Carlo move which compares one hard rod normal mode to a harmonic 
 * normal mode and explores the phase space defined by the remaining normal modes.
 * 
 * comparedWV is the wavevector that will be compared.  It cannot be selected
 * by the "change a random mode" section of the doTrial() method
 * 
 * @author cribbin
 *
 */
public class MCMoveCompareWVLoopLEFT extends MCMoveBoxStep{

    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    private final AtomIteratorLeafAtoms iterator;
    protected double[][] uOld;
    protected double[] deltaU;
    protected final IRandom random;
    protected double energyOld, energyNew, energyEvenLater; /*, latticeEnergy*/;
    protected final MeterPotentialEnergy energyMeter;
    private double[][][] eigenVectors;
    private IVectorMutable[] waveVectors;
    private double[][] omega2;
    int comparedWV;
    private double[] gaussian;
    protected double temperature;
    private double[][] stdDev;
    private double[] rRand, iRand, realT, imagT;
    private double[] waveVectorCoefficients;
    private double wvc;
    double[] uNow;
    
    
    public MCMoveCompareWVLoopLEFT(IPotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster);
        
        this.random = random;
        iterator = new AtomIteratorLeafAtoms();
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        gaussian = new double[2];
    }


    public boolean doTrial() {
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        rRand = new double[coordinateDim];
        iRand = new double[coordinateDim];
        realT = new double[coordinateDim];
        imagT = new double[coordinateDim];
        
        //nan These lines make it a single atom-per-molecule class.
        BasisCell cell = cells[0];
        uOld = new double[cells.length][coordinateDim];
        double normalization = 1/Math.sqrt(cells.length);

        //Get normal mode coordinate information
        coordinateDefinition.calcT(waveVectors[comparedWV], realT, imagT);
        
        
//ZERO OUT A NORMAL MODE.
        for(int iCell = 0; iCell < cells.length; iCell++){
            //store old positions.
            uNow = coordinateDefinition.calcU(cells[iCell].molecules);
            System.arraycopy(uNow, 0, uOld[iCell], 0, coordinateDim);
            cell = cells[iCell];
            //rezero deltaU
            for(int j = 0; j < coordinateDim; j++){
                deltaU[j] = 0.0;
            }
            
            //Calculate the contributions to the current position of the zeroed
            //mode, and subtract it from the overall position
            double kR = waveVectors[comparedWV].dot(cell.cellPosition);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);
            for(int i = 0; i < coordinateDim; i++){
                //Calculate the current coordinate:
                double realCoord = 0, imagCoord = 0;
                for (int j=0; j<coordinateDim; j++) {
                    realCoord += eigenVectors[comparedWV][i][j] * realT[j];
                    imagCoord += eigenVectors[comparedWV][i][j] * imagT[j];
                }
                for(int j = 0; j < coordinateDim; j++){
                    deltaU[j] -= wvc*eigenVectors[comparedWV][i][j] * 2.0 *
                        (realCoord*coskR - imagCoord*sinkR);
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
        energyOld = energyMeter.getDataAsScalar();
        if(Double.isInfinite(energyOld)){
            throw new IllegalStateException("Overlap after the removal of a mode!");
        }
        
//MOVE A RANDOM (N-1) MODE, AND MEASURE energyNew
        //equivalent to MCMoveChangeMode
        if(comparedWV != (waveVectorCoefficients.length)) {
            //Select the wave vector whose eigenvectors will be changed.
            //The the compared wavevector, and any wavevector
            //numbered higher than it, are.rejected as possibilities.
            int changedWV = random.nextInt(waveVectorCoefficients.length - comparedWV);
            changedWV += comparedWV;
//            System.out.println(changedWV);
            
            //calculate the new positions of the atoms.
            //loop over cells
            double[] delta = new double[coordinateDim*2];
            for ( int i = 0; i < coordinateDim*2; i++) {
                delta[i] = (2*random.nextDouble()-1) * stepSize;
            }
            for(int iCell = 0; iCell < cells.length; iCell++){
                uNow = coordinateDefinition.calcU(cells[iCell].molecules);
                cell = cells[iCell];
                //rezero deltaU
                for(int j = 0; j < coordinateDim; j++){
                    deltaU[j] = 0.0;
                }
                //loop over the wavevectors, and sum contribution of each to the
                //generalized coordinates.  Change the selected wavevectors eigen
                //-vectors at the same time!
                double kR = waveVectors[changedWV].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                for(int i = 0; i < coordinateDim; i++){
                    if( !(Double.isInfinite(omega2[changedWV][i])) ){
                        for(int j = 0; j < coordinateDim; j++){
                             deltaU[j] += waveVectorCoefficients[changedWV] * 
                                 eigenVectors[changedWV][i][j] * 2.0 * (delta[j]*coskR
                                 - delta[j+coordinateDim]*sinkR);
                        }
                    }
                }
                for(int i = 0; i < coordinateDim; i++) {
                    uNow[i] += deltaU[i];
                }
                coordinateDefinition.setToU(cells[iCell].molecules, uNow);
            }
        }
        energyNew = energyMeter.getDataAsScalar();
        
//MOVE THE NORMAL MODE THAT WAS ZEROED OUT.
        //set up the gaussian values
        double sqrtT = Math.sqrt(temperature);
        for (int j=0; j<coordinateDim; j++) {
            if (stdDev[comparedWV][j] == 0) continue;
            //generate real and imaginary parts of random normal-mode coordinate Q
            double realGauss = random.nextGaussian() * sqrtT;
            double imagGauss = random.nextGaussian() * sqrtT;
            
            //XXX we know that if c(k) = 0.5, one of the gaussians will be ignored, but
            // it's hard to know which.  So long as we don't put an atom at the origin
            // (which is true for 1D if c(k)=0.5), it's the real part that will be ignored.
            if (wvc == 0.5) imagGauss = 0;
            rRand[j] = realGauss * stdDev[comparedWV][j];
            iRand[j] = imagGauss * stdDev[comparedWV][j];
            gaussian[0] = realGauss;
            gaussian[1] = imagGauss;
            
        }

        //calculate the new positions of the atoms.
        for(int iCell = 0; iCell < cells.length; iCell++){
            uNow = coordinateDefinition.calcU(cells[iCell].molecules);
            cell = cells[iCell];
            //rezero deltaU
            for(int j = 0; j < coordinateDim; j++){
                deltaU[j] = 0.0;
            }
            // Calculate the change in position due to the substitution of a 
            //  Gaussian.
            double kR = waveVectors[comparedWV].dot(cell.cellPosition);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);
            for(int i = 0; i < coordinateDim; i++){
                for(int j = 0; j < coordinateDim; j++){
                    deltaU[j] += wvc*eigenVectors[comparedWV][i][j] * 2.0 *
                        (rRand[i]*coskR - iRand[i]*sinkR);
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
        
        return true;
    }
    
    public double getA() {
        return 1;
    }

    public double getB() {
        return -(energyNew - energyOld);
    }
    
    public void acceptNotify() {
    }

    public double energyChange() {
        return energyNew - energyOld;
    }

    public void rejectNotify() {
//        System.out.println("reject " + energyNew);
        // Set all the atoms back to the old values of u
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        for (int iCell = 0; iCell<cells.length; iCell++) {
            BasisCell cell = cells[iCell];
            coordinateDefinition.setToU(cell.molecules, uOld[iCell]);
        }
    }

    public void setBox(IBox newBox) {
        super.setBox(newBox);
        iterator.setBox(newBox);
        energyMeter.setBox(newBox);
    }

    public AtomIterator affectedAtoms() {return iterator;}
    
    public void setCoordinateDefinition(CoordinateDefinition newCD) {
        coordinateDefinition = newCD;
        deltaU = new double[coordinateDefinition.getCoordinateDim()];
        uOld = null;
        realT = new double[coordinateDefinition.getCoordinateDim()];
        imagT = new double[coordinateDefinition.getCoordinateDim()];
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    /**
     * Set the wave vectors accessible to the move.
     * 
     * @param wv
     */
    public void setWaveVectors(IVectorMutable[] wv){
        waveVectors = new IVectorMutable[wv.length];
        waveVectors = wv;
    }    
    public void setWaveVectorCoefficients(double[] newWaveVectorCoefficients) {
        waveVectorCoefficients = newWaveVectorCoefficients;
    }
    /**
     * Informs the move of the eigenvectors 
     */
    public void setEigenVectors(double[][][] newEigenVectors) {
        eigenVectors = newEigenVectors;
    }
    
    public void setOmegaSquared(double[][] o2, double[] coeff) {
        this.omega2 = o2;
        stdDev = new double[o2.length][o2[0].length];
        for (int i=0; i<stdDev.length; i++) {
            for (int j=0; j<stdDev[i].length; j++) {
                stdDev[i][j] = Math.sqrt(1.0/(2.0*o2[i][j]*coeff[i]));
            }
        }
    }
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    public double[] getGaussian(){
        return gaussian;
    }
    public void setComparedWV(int wv){
        if(wv == waveVectorCoefficients.length) {System.out.println("3D system is now entirely Gaussian!");};
        comparedWV = wv;
        wvc = waveVectorCoefficients[wv];
    }
    public int getComparedWV(){
        return comparedWV;
    }
    public double[] getLastGaussian(){
        return gaussian;
    }

}