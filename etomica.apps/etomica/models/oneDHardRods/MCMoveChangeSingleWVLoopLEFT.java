package etomica.models.oneDHardRods;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.atom.Atom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinition.BasisCell;

public class MCMoveChangeSingleWVLoopLEFT extends MCMoveBoxStep {
    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    private final AtomIteratorLeafAtoms iterator;
    protected double[][] uOld;
    protected double[] deltaU;
    protected final IRandom random;
    protected double energyOld, energyNew /*, latticeEnergy*/;
    protected final MeterPotentialEnergy energyMeter;
    private double[][][] eigenVectors;
    private double[][] omega2;
    private IVectorMutable[] waveVectors;
    private double[] waveVectorCoefficients;
    int changedWV, harmonicWV;  //all wvs from the harmonic wv and up are not changed.
    
    
    public MCMoveChangeSingleWVLoopLEFT(IPotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster);
        
        this.random = random;
        iterator = new AtomIteratorLeafAtoms();
        energyMeter = new MeterPotentialEnergy(potentialMaster);
    }

    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        coordinateDefinition = newCoordinateDefinition;
        deltaU = new double[coordinateDefinition.getCoordinateDim()];
        uOld = null;
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }
    
    /**
     * The harmonic wavevector and all wavevectors with lower numbers are not
     * able to be changed by this MCMove.
     */
    public void setHarmonicWV(int hwv){
        harmonicWV = hwv;
    }

    /**
     * Set the wave vectors used by the move.
     * 
     * @param wv
     */
    public void setWaveVectors(IVectorMutable[] wv){
        waveVectors = new IVectorMutable[wv.length];
        waveVectors = wv;
    }
    public void setWaveVectorCoefficients(double[] coeff){
        waveVectorCoefficients = coeff;
    }
    /**
     * Informs the move of the eigenvectors for the selected wave vector.  The
     * actual eigenvectors used will be those specified via setModes
     */
    public void setEigenVectors(double[][][] newEigenVectors) {
        eigenVectors = newEigenVectors;
    }
    
    public void setOmegaSquared(double[][] o2){
        omega2 = o2;
    }
    
    public void setBox(IBox newBox) {
        super.setBox(newBox);
        iterator.setBox(newBox);
        energyMeter.setBox(newBox);
    }

    public AtomIterator affectedAtoms() {
        return iterator;
    }

//    public void setWaveVectorAndEigenVectorsChanged(int wv, int[] evectors){
//        //we will need some flag to indicate that these were set, and not
          //randomly assign them.
//    }
//    
    public boolean doTrial() {
//        System.out.println("MCMoveChangeMode doTrial");
        energyOld = energyMeter.getDataAsScalar();
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        
        //nan These lines make it a single atom-per-molecule class, and
        // assume that the first cell is the same as every other cell.
        BasisCell cell = cells[0];
//        double[] calcedU = coordinateDefinition.calcU(cell.molecules);
        uOld = new double[cells.length][coordinateDim];
        
        // Select the wave vector whose eigenvectors will be changed.
        //The zero wavevector is center of mass motion, and is rejected as a 
        //possibility. 
        changedWV = random.nextInt(waveVectorCoefficients.length  - harmonicWV);
        changedWV += harmonicWV;
//        System.out.println(changedWV);
        
        //calculate the new positions of the atoms.
        //loop over cells
        double[] delta = new double[coordinateDim*2];
        for ( int i = 0; i < coordinateDim*2; i++) {
            delta[i] = (2*random.nextDouble()-1) * stepSize;
        }
        for(int iCell = 0; iCell < cells.length; iCell++){
            //store old positions.
            double[] uNow = coordinateDefinition.calcU(cells[iCell].molecules);
            System.arraycopy(uNow, 0, uOld[iCell], 0, coordinateDim);
            cell = cells[iCell];
            for(int i = 0; i< coordinateDim; i++){
                  deltaU[i] = 0;
            }
            
            //loop over the wavevectors, and sum contribution of each to the
            //generalized coordinates.  Change the selected wavevectors eigen-
            //vectors at the same time!
            double kR = waveVectors[changedWV].dot(cell.cellPosition);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);
            for(int iMode = 0; iMode < coordinateDim; iMode++){
                if ( !(Double.isInfinite(omega2[changedWV][iMode])) ){
                    for(int j = 0; j < coordinateDim; j++){
                        deltaU[j] += eigenVectors[changedWV][iMode][j]*2.0*delta[iMode]*coskR;
                        if(waveVectorCoefficients[changedWV] == 1.0){
                            deltaU[j] -= 2.0*eigenVectors[changedWV][iMode][j] *delta[iMode+coordinateDim]*sinkR;
                        }
                    }
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
        
        energyNew = energyMeter.getDataAsScalar();
        return true;
    }
    
    public double getA() {
        return 1;
    }

    public double getB() {
        return -(energyNew - energyOld);
    }
    
    public void acceptNotify() {
//        System.out.println("accept MCMoveChangeSingleWVLoopLEFT");
    }

    public double energyChange() {
        return energyNew - energyOld;
    }

    public void rejectNotify() {
//        System.out.println("reject MCMoveChangeSingleWVLoopLEFT");
        // Set all the atoms back to the old values of u
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        for (int iCell = 0; iCell<cells.length; iCell++) {
            BasisCell cell = cells[iCell];
            coordinateDefinition.setToU(cell.molecules, uOld[iCell]);
        }
    }
    
    private void printLocations(){
        IAtomList list = box.getLeafList();
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        int ats = box.getLeafList().getAtomCount();
        
        if(box.getBoundary().getEdgeVector(0).getD() == 1){
            for(int i = 0; i < ats; i++){
                System.out.println(i + "  " + ((Atom)list.getAtom(i)).getPosition().getX(0));
            }
        }
        
        if(box.getBoundary().getEdgeVector(0).getD() == 3){
            for(int i = 0; i < ats; i++){
                System.out.println("Atom " + i);
                for(int j = 0; j < 3; j++){
                    System.out.println(j + " " + ((Atom)list.getAtom(i)).getPosition().getX(j));
                }
            }
        }
    }
    
}