package etomica.virial;

import etomica.action.AtomAction;
import etomica.atom.AtomSet;
import etomica.atom.IAtomGroup;
import etomica.atom.IAtomPositioned;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.space.IVector;
import etomica.util.IRandom;

/**
 * MCMove for use in a Mayer sampling simulation that rotates all molecules in
 * a Phase except the first molecule, which is never moved.  The angle of
 * rotation is the step size and can be tuned for some acceptance rate.
 */
public class MCMoveClusterRotateMoleculeMulti extends MCMoveRotateMolecule3D {

    /**
     * @param potentialMaster
     * @param space
     */
    public MCMoveClusterRotateMoleculeMulti(PotentialMaster potentialMaster,
            IRandom random, int numMolecules) {
        super(potentialMaster, random);
        weightMeter = new MeterClusterWeight(potential);
        setName("MCMoveClusterMolecule");
        nMolecules = numMolecules;
        selectedMolecules = new IAtomGroup[nMolecules];
        oldPositions = new IVector[nMolecules][];
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        weightMeter.setPhase(p);
        selectMolecules();
        for (int i=0; i<nMolecules; i++) {
            molecule = selectedMolecules[i];
            oldPositions[i] = new IVector[molecule.getChildList().getAtomCount()];
            for (int j=0; j<oldPositions[i].length; j++) {
                oldPositions[i][j] = p.getSpace().makeVector();
            }
        }
    }

    public boolean doTrial() {
        uOld = weightMeter.getDataAsScalar();
        boolean doRelax = false;
        if (trialCount-- == 0) {
            doRelax = true;
            trialCount = relaxInterval;
        }
        for (int i=0; i<selectedMolecules.length; i++) {
            molecule = selectedMolecules[i];
            leafAtomIterator.setRootAtom(molecule);
            leafAtomIterator.reset();
            r0.E(molecule.getType().getPositionDefinition().position(molecule));
        
            double dTheta = (2*random.nextDouble() - 1.0)*stepSize;
            rotationTensor.setAxial(random.nextInt(3),dTheta);
            
            int j=0;
            for (IAtomPositioned a = (IAtomPositioned)leafAtomIterator.nextAtom(); a != null;
                 a = (IAtomPositioned)leafAtomIterator.nextAtom()) {
                oldPositions[i][j++].E(a.getPosition());
            }
            leafAtomIterator.reset();
            doTransform();
            
            
            if (doRelax && relaxAction != null) {
                relaxAction.setAtom(molecule);
                relaxAction.actionPerformed();
            }
        }

        ((PhaseCluster)phase).trialNotify();
        uNew = weightMeter.getDataAsScalar();
        return true;
    }
    
    public void selectMolecules() {
        AtomSet atomList = ((IAtomGroup)phase.getSpeciesMaster().getAgentList().getAtom(0)).getChildList();
        for (int i=0; i<selectedMolecules.length; i++) {
            selectedMolecules[i] = (IAtomGroup)atomList.getAtom(i+1);
        }
    }

    public double getB() {
        return 0.0;
    }
    
    public double getA() {
        return (uOld==0.0) ? Double.POSITIVE_INFINITY : uNew/uOld;
    }
    
    public void acceptNotify() {
        super.acceptNotify();
        ((PhaseCluster)phase).acceptNotify();
    }
    
    public void rejectNotify() {
        for (int i=0; i<selectedMolecules.length; i++) {
            molecule = selectedMolecules[i];
            leafAtomIterator.setRootAtom(molecule);
            leafAtomIterator.reset();
            int j=0;
            for (IAtomPositioned a = (IAtomPositioned)leafAtomIterator.nextAtom(); a != null;
                 a = (IAtomPositioned)leafAtomIterator.nextAtom()) {
                a.getPosition().E(oldPositions[i][j++]);
            }
        }
        ((PhaseCluster)phase).rejectNotify();
    }
    
    public void setRelaxAction(AtomAction action) {
        relaxAction = action;
    }
    
    private static final long serialVersionUID = 1L;
    private final MeterClusterWeight weightMeter;
    private final IAtomGroup[] selectedMolecules;
    private final IVector[][] oldPositions;
    private final int nMolecules;
    private int trialCount, relaxInterval = 100;
    private AtomAction relaxAction;
}
