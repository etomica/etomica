package etomica.models.oneDHardRods;

import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;

/**
 * A Monte Carlo move which selects a wave vector, and an eigenvector allowed 
 * by that wave vector.
 * 
 * @author cribbin
 *
 */
public class MCMoveConvertMode extends MCMoveBoxStep{

    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    private final AtomIteratorAllMolecules iterator;
    protected double[][] uOld;
    protected double[] deltaU;
    protected final IRandom random;
    protected double energyOld, energyNew /*, latticeEnergy*/;
    protected final MeterPotentialEnergy energyMeter;
    private double[][][] eigenVectors;
    private IVector[] waveVectors;
    int convertedWV, changedWV;
    private double[] gaussian;
    protected double temperature;
    private double[][] stdDev;
    private double[] rRand, iRand, realT, imagT;
    private double[] waveVectorCoefficients;
    private double wvc;
    public long count;
    double[] uNow;

    public MCMoveConvertMode(IPotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster);
        
        this.random = random;
        iterator = new AtomIteratorAllMolecules();
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        gaussian = new double[2];
        count = 0;
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
        
        //Get normal mode coordinate information
        coordinateDefinition.calcT(waveVectors[convertedWV], realT, imagT);
        
//ZERO OUT A NORMAL MODE, AND MEASURE energyOld
        for(int iCell = 0; iCell < cells.length; iCell++){
            //store old positions.
            uNow = coordinateDefinition.calcU(cells[iCell].molecules);
            System.arraycopy(uNow, 0, uOld[iCell], 0, coordinateDim);
            cell = cells[iCell];
            
            //Calculate the contributions to the current position of the zeroed
            //mode, and subtract it from the overall position
            double kR = waveVectors[convertedWV].dot(cell.cellPosition);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);
            for(int i = 0; i < coordinateDim; i++){
                //Calculate the current coordinate:
                double realCoord = 0, imagCoord = 0;
                double[] coordToChange = new double[coordinateDim];
                
                for (int j=0; j<coordinateDim; j++) {
                    realCoord += eigenVectors[convertedWV][i][j] * realT[j];
                    System.out.println("realcoord " + realCoord);
                    imagCoord += eigenVectors[convertedWV][i][j] * imagT[j];
                    System.out.println("imagcoord " + imagCoord);
                    coordToChange[j] = 0.0;
                }
                for(int j = 0; j < coordinateDim; j++){
//                    coordToChange[j] += wvc*eigenVectors[convertedWV][i][j] * 2.0 *
//                        (realCoord*coskR - imagCoord*sinkR);
//                    deltaU[j] = -coordToChange[j];
                    
                    deltaU[j] = -wvc*eigenVectors[convertedWV][i][j] * 2.0 *
                         (realCoord*coskR - imagCoord*sinkR);
                    System.out.println(deltaU[j]);
                }
            }
            
            for(int i = 0; i < coordinateDim; i++) {
                uNow[i] += deltaU[i];
                System.out.println("i: " + i + " "+uNow[i]);
            }
            coordinateDefinition.setToU(cells[iCell].molecules, uNow);
            
        }
        energyOld = energyMeter.getDataAsScalar();
        System.out.println("NRGOLD: "+ energyOld);

//MOVE A RANDOM GAUSSIAN MODE, AND MEASURE energyNew
        if(convertedWV != 1) {
            //Select the wave vector whose eigenvectors will be changed.
            //The zero wavevector is center of mass motion, and is rejected as a 
            //possibility, as is the converted wavevector.
            changedWV = random.nextInt(convertedWV-2);
            changedWV +=1;
            
            //calculate the new positions of the atoms.
            //loop over cells
            double delta1 = (2*random.nextDouble()-1) * stepSize;
            double delta2 = (2*random.nextDouble()-1) * stepSize;
            delta1 = 0.25; delta2 = 0.5;
            for(int iCell = 0; iCell < cells.length; iCell++){
                cell = cells[iCell];
                
                //loop over the wavevectors, and sum contribution of each to the
                //generalized coordinates.  Change the selected wavevectors eigen
                //-vectors at the same time!
                double kR = waveVectors[changedWV].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                for(int i = 0; i < coordinateDim; i++){
                    for(int j = 0; j < coordinateDim; j++){
                         deltaU[j] += eigenVectors[changedWV][i][j]*2.0*
                             (delta1*coskR - delta2*sinkR);
                    }
                }
                double normalization = 1/Math.sqrt(cells.length);
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
            if (stdDev[convertedWV][j] == 0) continue;
            //generate real and imaginary parts of random normal-mode coordinate Q
            double realGauss = random.nextGaussian() * sqrtT;
            double imaginaryGauss = random.nextGaussian() * sqrtT;
            rRand[j] = realGauss * stdDev[convertedWV][j];
            iRand[j] = imaginaryGauss * stdDev[convertedWV][j];
            gaussian[0] = realGauss;
            gaussian[1]= imaginaryGauss;
            //XXX we know that if c(k) = 0.5, one of the gaussians will be ignored, but
            // it's hard to know which.  So long as we don't put an atom at the origin
            // (which is true for 1D if c(k)=0.5), it's the real part that will be ignored.
            if (wvc == 0.5) imaginaryGauss = 0;
        }
        

        //calculate the new positions of the atoms.
        for(int iCell = 0; iCell < cells.length; iCell++){
            cell = cells[iCell];
            
            // Calculate the change in position due to the substitution of a 
            //  Gaussian.
            double kR = waveVectors[convertedWV].dot(cell.cellPosition);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);
            for(int i = 0; i < coordinateDim; i++){
                //Calculate the current coordinate:
                double realCoord = 0, imaginaryCoord = 0;
                double[] coordToChange = new double[coordinateDim];
                
                for (int j=0; j<coordinateDim; j++) {
                    realCoord += eigenVectors[convertedWV][i][j] * realT[j];
                    imaginaryCoord += eigenVectors[convertedWV][i][j] * imagT[j];
                    coordToChange[j] = 0.0;
                }
                for(int j = 0; j < coordinateDim; j++){
                    coordToChange[j] += wvc*eigenVectors[convertedWV][i][j] * 2.0 *
                        (realCoord*coskR - imaginaryCoord*sinkR);
                    deltaU[j] += wvc*eigenVectors[convertedWV][i][j] * 2.0 *
                        (rRand[i]*coskR - iRand[i]*sinkR);
                    deltaU[j] -= coordToChange[j];
                }
            }
            double normalization = 1/Math.sqrt(cells.length);
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
        System.out.println("accept"+ count++);
    }

    public double energyChange() {
        return energyNew - energyOld;
    }

    public void rejectNotify() {
        System.out.println("reject");
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
    public void setWaveVectors(IVector[] wv){
        waveVectors = new IVector[wv.length];
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
    public void setConvertedWaveVector(int wv){
        if(wv == 1) {System.out.println("System is now entirely Gaussian!");};
        convertedWV = wv;
        wvc = waveVectorCoefficients[wv];
        coordinateDefinition.calcT(waveVectors[wv], realT, imagT);
    }
    public int getConvertedWaveVector(){
        return convertedWV;
    }
    public double[] getLastGaussian(){
        return gaussian;
    }

}