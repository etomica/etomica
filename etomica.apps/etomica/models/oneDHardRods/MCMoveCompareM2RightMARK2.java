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
public class MCMoveCompareM2RightMARK2 extends MCMoveBoxStep{

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
    int comparedWV1, comparedWV2;
    private double[] gaussian;
    protected double temperature;
    private double[][] stdDev;
    private double[] rRand1, iRand1, rRand2, iRand2, realT1, imagT1, realT2, imagT2;
    private double[] waveVectorCoefficients;
    double[] uNow;
    
;
    
    public MCMoveCompareM2RightMARK2(IPotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster);
        
        this.random = random;
        iterator = new AtomIteratorLeafAtoms();
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        gaussian = new double[4];
//        count = 0;
    }


    public boolean doTrial() {
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        rRand1 = new double[coordinateDim];
        iRand1 = new double[coordinateDim];
        rRand2 = new double[coordinateDim];
        iRand2 = new double[coordinateDim];
        realT1 = new double[coordinateDim];
        imagT1 = new double[coordinateDim];
        realT2 = new double[coordinateDim];
        imagT2 = new double[coordinateDim];
        
        //nan These lines make it a single atom-per-molecule class.
        BasisCell cell = cells[0];
        uOld = new double[cells.length][coordinateDim];
        double normalization = 1/Math.sqrt(cells.length);

        //Get normal mode coordinate information
        coordinateDefinition.calcT(waveVectors[comparedWV1], realT1, imagT1);
        
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
            double kR = waveVectors[comparedWV1].dot(cell.cellPosition);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);
            for(int i = 0; i < coordinateDim; i++){
                //Calculate the current coordinate:
                double realCoord = 0, imagCoord = 0;
                for (int j=0; j<coordinateDim; j++) {
                    realCoord += eigenVectors[comparedWV1][i][j] * realT1[j];
                    imagCoord += eigenVectors[comparedWV1][i][j] * imagT1[j];
                }
                for(int j = 0; j < coordinateDim; j++){
                    deltaU[j] -= waveVectorCoefficients[comparedWV1]*eigenVectors[comparedWV1][i][j] * 2.0 *
                        (realCoord*coskR - imagCoord*sinkR);
                }
            }
            
            kR = waveVectors[comparedWV2].dot(cell.cellPosition);
            coskR = Math.cos(kR);
            sinkR = Math.sin(kR);
            for(int i = 0; i < coordinateDim; i++){
                //Calculate the current coordinate:
                double realCoord = 0, imagCoord = 0;
                for (int j=0; j<coordinateDim; j++) {
                    realCoord += eigenVectors[comparedWV2][i][j] * realT2[j];
                    imagCoord += eigenVectors[comparedWV2][i][j] * imagT2[j];
                }
                for(int j = 0; j < coordinateDim; j++){
                    deltaU[j] -= waveVectorCoefficients[comparedWV2]*eigenVectors[comparedWV2][i][j] * 2.0 *
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
        if(energyOld != 0.0){
            for(int k = 0; k < waveVectors.length; k++){
                System.out.println(k + " " +((IAtomPositioned)coordinateDefinition.getBox().getLeafList().getAtom(k)).getPosition());
            }
            throw new IllegalStateException("Overlap after the removal of a mode!");
        }
        
//MOVE A RANDOM (N-1) MODE, AND MEASURE energyNew
        //equivalent to MCMoveChangeMode
        if(comparedWV1 != 1) {
            //Select the wave vector whose eigenvectors will be changed.
            //The zero wavevector is center of mass motion, and is rejected as a 
            //possibility, as is the compared wavevector and any wavevector
            //number higher than it.
            int changedWV = random.nextInt(comparedWV1-1);
            changedWV += 1;
            
//            System.out.println(changedWV);
            
            //calculate the new positions of the atoms.
            //loop over cells
            double delta1 = (2*random.nextDouble()-1) * stepSize;
            double delta2 = (2*random.nextDouble()-1) * stepSize;
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
                    for(int j = 0; j < coordinateDim; j++){
                         deltaU[j] += waveVectorCoefficients[changedWV] * 
                             eigenVectors[changedWV][i][j] * 2.0 * (delta1*coskR
                             - delta2*sinkR);
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
        }
        energyNew = energyMeter.getDataAsScalar();
        
//MOVE THE NORMAL MODE THAT WAS ZEROED OUT.
        //set up the gaussian values
        double sqrtT = Math.sqrt(temperature);
        for (int j=0; j<coordinateDim; j++) {
            if (stdDev[comparedWV1][j] == 0) continue;
            //generate real and imaginary parts of random normal-mode coordinate Q
            double realGauss = random.nextGaussian() * sqrtT;
            double imagGauss = random.nextGaussian() * sqrtT;
            
            //XXX we know that if c(k) = 0.5, one of the gaussians will be ignored, but
            // it's hard to know which.  So long as we don't put an atom at the origin
            // (which is true for 1D if c(k)=0.5), it's the real part that will be ignored.
            if (waveVectorCoefficients[comparedWV1] == 0.5) imagGauss = 0;
            rRand1[j] = realGauss * stdDev[comparedWV1][j];
            iRand1[j] = imagGauss * stdDev[comparedWV1][j];
            gaussian[0] = realGauss;
            gaussian[1] = imagGauss;
            
            
            if (stdDev[comparedWV2][j] == 0) continue;
            //generate real and imaginary parts of random normal-mode coordinate Q
            realGauss = random.nextGaussian() * sqrtT;
            imagGauss = random.nextGaussian() * sqrtT;
                        
            //XXX we know that if c(k) = 0.5, one of the gaussians will be ignored, but
            // it's hard to know which.  So long as we don't put an atom at the origin
            // (which is true for 1D if c(k)=0.5), it's the real part that will be ignored.
            if (waveVectorCoefficients[comparedWV2] == 0.5) imagGauss = 0;
            rRand2[j] = realGauss * stdDev[comparedWV2][j];
            iRand2[j] = imagGauss * stdDev[comparedWV2][j];
            gaussian[2] = realGauss;
            gaussian[3] = imagGauss;
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
            double kR = waveVectors[comparedWV1].dot(cell.cellPosition);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);
            for(int i = 0; i < coordinateDim; i++){
                for(int j = 0; j < coordinateDim; j++){
                    deltaU[j] += waveVectorCoefficients[comparedWV1] *
                        eigenVectors[comparedWV1][i][j] * 2.0 *
                        (rRand1[i]*coskR - iRand1[i]*sinkR);
                }
            }
            
            kR = waveVectors[comparedWV2].dot(cell.cellPosition);
            coskR = Math.cos(kR);
            sinkR = Math.sin(kR);
            for(int i = 0; i < coordinateDim; i++){
                for(int j = 0; j < coordinateDim; j++){
                    deltaU[j] += waveVectorCoefficients[comparedWV2] *
                        eigenVectors[comparedWV2][i][j] * 2.0 *
                        (rRand2[i]*coskR - iRand2[i]*sinkR);
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
        System.out.println("accept MCMoveCompareM2RightMARK2");
//        count ++;
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
        realT1 = new double[coordinateDefinition.getCoordinateDim()];
        imagT1 = new double[coordinateDefinition.getCoordinateDim()];
        realT2 = new double[coordinateDefinition.getCoordinateDim()];
        imagT2 = new double[coordinateDefinition.getCoordinateDim()];
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
    
    public void setOmegaSquared(double[][] omega2, double[] coeff) {
        stdDev = new double[omega2.length][omega2[0].length];
        for (int i=0; i<stdDev.length; i++) {
            for (int j=0; j<stdDev[i].length; j++) {
                stdDev[i][j] = Math.sqrt(1.0/(2.0*omega2[i][j]*coeff[i]));
            }
        }
    }
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    public double[] getGaussian(){
        return gaussian;
    }
    public void setComparedWV(int wv, int wv2){
        if(wv == 1 && wv2 == 2) {System.out.println("System is now entirely Gaussian!");};
        comparedWV1 = wv;
        comparedWV2 = wv2;
    }

    public double[] getLastGaussian(){
        return gaussian;
    }

}