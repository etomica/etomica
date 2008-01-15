package etomica.virial;

import etomica.atom.AtomSet;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.potential.PotentialMaster;
import etomica.simulation.ISimulation;
import etomica.space.IVectorRandom;
import etomica.util.IRandom;

/**
 * @author kofke
 *
 * Extension of MCMoveAtom that does trial in which several atom positions are
 * perturbed.  However, position of first atom is never altered.  
 */
public class MCMoveClusterMoleculeMulti extends MCMoveMolecule {

    private static final long serialVersionUID = 1L;
    private final MeterClusterWeight weightMeter;

    public MCMoveClusterMoleculeMulti(ISimulation sim, PotentialMaster potentialMaster) {
    	this(potentialMaster,sim.getRandom(), 1.0);
    }
    
    /**
     * Constructor for MCMoveAtomMulti.
     * @param parentIntegrator
     * @param nAtoms number of atoms to move in a trial.  Number of atoms in
     * box should be at least one greater than this value (greater
     * because first atom is never moved)
     */
    public MCMoveClusterMoleculeMulti(PotentialMaster potentialMaster,
            IRandom random, double stepSize) {
        super(potentialMaster,random,stepSize,Double.POSITIVE_INFINITY,false);
        weightMeter = new MeterClusterWeight(potential);
    }

    public void setBox(Box p) {
        super.setBox(p);
        weightMeter.setBox(p);
        translationVectors = new IVectorRandom[box.getMoleculeList().getAtomCount()-1];
        for (int i=0; i<box.getMoleculeList().getAtomCount()-1; i++) {
            translationVectors[i] = (IVectorRandom)potential.getSpace().makeVector();
        }
    }
    
    //note that total energy is calculated
    public boolean doTrial() {
        uOld = weightMeter.getDataAsScalar();
        AtomSet moleculeList = box.getMoleculeList();
        for(int i=1; i<moleculeList.getAtomCount(); i++) {
            translationVectors[i-1].setRandomCube(random);
            translationVectors[i-1].TE(stepSize);
            groupTranslationVector.E(translationVectors[i-1]);
            moveMoleculeAction.actionPerformed(moleculeList.getAtom(i));
        }
        ((BoxCluster)box).trialNotify();
        return true;
    }
	
    public void rejectNotify() {
        AtomSet moleculeList = box.getMoleculeList();
        for(int i=1; i<moleculeList.getAtomCount(); i++) {
            groupTranslationVector.E(translationVectors[i-1]);
            groupTranslationVector.TE(-1);
            moveMoleculeAction.actionPerformed(moleculeList.getAtom(i));
        }
        ((BoxCluster)box).rejectNotify();
        if (weightMeter.getDataAsScalar() == 0) {
            throw new RuntimeException("oops oops, reverted to illegal configuration");
        }
    }

    public void acceptNotify() {
        if (uNew == 0) {
            throw new RuntimeException("oops, accepted illegal configuration");
        }
        ((BoxCluster)box).acceptNotify();
        if (weightMeter.getDataAsScalar() == 0) {
            throw new RuntimeException("oops oops, accepted illegal configuration");
        }
    }
    
    public double getB() {
        return 0.0;
    }
    
    public double getA() {
        uNew = weightMeter.getDataAsScalar();
        return uNew/uOld;
    }
	
    private IVectorRandom[] translationVectors;
}
