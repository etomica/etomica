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
public class MCMoveChangeMode extends MCMoveBoxStep{

    private static final long serialVersionUID = 1L;
    protected CoordinateDefinition coordinateDefinition;
    private final AtomIteratorAllMolecules iterator;
    protected double[][] uOld;
    protected double[] u;
    protected final IRandom random;
    protected double energyOld, energyNew /*, latticeEnergy*/;
    protected final MeterPotentialEnergy energyMeter;
    private double[][][] eigenVectors;
    private IVector[] waveVectors;
    int changedWV, changedEV;
    
    public MCMoveChangeMode(IPotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster);
        
        this.random = random;
        iterator = new AtomIteratorAllMolecules();
        energyMeter = new MeterPotentialEnergy(potentialMaster);
    }

    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        coordinateDefinition = newCoordinateDefinition;
        u = new double[coordinateDefinition.getCoordinateDim()];
        uOld = null;         //
    }
    
    public CoordinateDefinition getCoordinateDefinition() {
        return coordinateDefinition;
    }

    /**
     * Set the wave vectors used by the move.
     * 
     * @param wv
     */
    public void setWaveVectors(IVector[] wv){
        waveVectors = new IVector[wv.length];
        waveVectors = wv;
    }
    /**
     * Informs the move of the eigenvectors for the selected wave vector.  The
     * actual eigenvectors used will be those specified via setModes
     */
    public void setEigenVectors(double[][][] newEigenVectors) {
        eigenVectors = newEigenVectors;
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
//    	System.out.println("called");
        energyOld = energyMeter.getDataAsScalar();
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        
        //nan These lines make it a single atom-per-molecule class, and
        // assumes that the first cell is the same as every other cell.
//        BasisCell cell = cells[0];
        double sqrtCells = Math.sqrt(cells.length);
//        double[] calcedU = coordinateDefinition.calcU(cell.molecules);
        uOld = new double[cells.length][coordinateDim];
        
        // Select the wave vector whose eigenvectors will be changed.
        changedWV = random.nextInt(waveVectors.length);
        // Select the eigenvector that will be changed.
        changedEV = random.nextInt(eigenVectors[changedWV].length);
        
        //calculate the new positions of the atoms.
        //loop over cells
        for(int iCell = 0; iCell < cells.length; iCell++){
            //store old positions.
            double[] uNow = coordinateDefinition.calcU(cells[iCell].molecules);
            System.arraycopy(uNow, 0, uOld[iCell], 0, coordinateDim);
            BasisCell cell = cells[iCell];
            System.out.println("Old");
            for(int i = 0; i< coordinateDim; i++){
            	System.out.println(u[i]);
                u[i] = 0;
            }
            
            //loop over the wavevectors, and sum contribution of each to the
            //generalized coordinates.  Change the selected wavevector's eigen-
            //vectors at the same time!
//            for(int iVector = 0; iVector < waveVectors.length; iVector++){
            int iVector = 1;
                double kR = waveVectors[iVector].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                for(int i = 0; i < coordinateDim; i++){
                    for(int j = 0; j < coordinateDim; j++){
                        double delta1 = (2*random.nextDouble()-1) * stepSize;
                        double delta2 = (2*random.nextDouble()-1) * stepSize;
                        
                        delta1 = 5.0;
                        delta2 = 1.0;
                        u[j] += eigenVectors[iVector][i][j]*2.0*(delta1*coskR - delta2*sinkR);
//                        System.out.println("iCell: "+ iCell+" i: "+i+ " j: "+j+" u: "+u[j]+ " ev: "+ eigenVectors[iVector][i][j]);
                    }
                }
//            }
            double normalization = 1/Math.sqrt(cells.length);
            System.out.println("prenormal");
            for(int i = 0; i < coordinateDim; i++){
                System.out.println(u[i]);
                u[i] *= normalization;
            }
            
            coordinateDefinition.setToU(cells[iCell].molecules, u);
            
            System.out.println("Postnormal");
            for(int i = 0; i < coordinateDim; i++){
                System.out.println(u[i]);
            }
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
        System.out.println("accept");
    }

    public double energyChange() {
        return energyNew - energyOld;
    }

    public void rejectNotify() {
//        System.out.println("reject");
        // Set all the atoms back to the old values of u
        BasisCell[] cells = coordinateDefinition.getBasisCells();
        for (int iCell = 0; iCell<cells.length; iCell++) {
            BasisCell cell = cells[iCell];
            coordinateDefinition.setToU(cell.molecules, uOld[iCell]);
        }
    }

}