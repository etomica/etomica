package etomica.models.hexane;

import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.integrator.mcmove.MCMovePhase;
import etomica.integrator.mcmove.MCMoveTracker;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Vector;

public class MCMoveHarmonic extends MCMovePhase {

    public MCMoveHarmonic(PotentialMaster potentialMaster) {
        super(potentialMaster, new MCMoveTracker());
        atomIterator = new AtomIteratorLeafAtoms();
    }
    
    public void setLatticeSites(Vector[] newLatticeSites) {
        latticeSites = newLatticeSites;
    }
    
    public void setEigenValues(double[] newEigenValues) {
        eigenValues = new double[newEigenValues.length];
        for (int i=0; i<eigenValues.length; i++) {
            eigenValues[i] = 0.5/newEigenValues[i];
        }
    }
    
    public void setEigenVectors(double[][] newEigenVectors) {
        eigenVectors = newEigenVectors;
    }
    
    public void setPhase(Phase newPhase) {
        super.setPhase(newPhase);
        atomIterator.setPhase(newPhase);
    }

    public AtomIterator affectedAtoms() {
        return atomIterator;
    }

    public boolean doTrial() {
        atomIterator.reset();
        while (atomIterator.hasNext()) {
            AtomLeaf atom = (AtomLeaf)atomIterator.nextAtom();
            Vector latticeSite = latticeSites[atom.getGlobalIndex()];
            atom.coord.position().E(latticeSite);
        }
        int D = potential.getSpace().D();
        Vector temp = potential.getSpace().makeVector();
        for (int i=0; i<eigenValues.length; i++) {
            double rand = Simulation.random.nextGaussian() * eigenValues[i];
            
            atomIterator.reset();
            int atomCount = 0;
            while (atomIterator.hasNext()) {
                temp.E(0);
                AtomLeaf atom = (AtomLeaf)atomIterator.nextAtom();
                int offset = atomCount * D;
                for (int j=0; j<D; j++) {
                    temp.setX(j, eigenVectors[i][offset+j]);
                }
                atom.coord.position().PEa1Tv1(rand,temp);
            }
        }
        return true;
    }

    public double getA() {
        // return 1 to guarantee success
        return 1;
    }

    public double getB() {
        // return 0 to guarantee success
        return 0;
    }

    public void acceptNotify() {
    }

    public double energyChange() {
        return 0;
    }

    public void rejectNotify() {
        throw new RuntimeException("This move should never be rejected");
    }

    private static final long serialVersionUID = 1L;
    private final AtomIteratorLeafAtoms atomIterator;
    private Vector[] latticeSites;
    private double[] eigenValues;
    private double[][] eigenVectors;
}
