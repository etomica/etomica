package etomica.virial;

import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.util.IRandom;

/**
 * Standard Monte Carlo molecule-displacement trial move for cluster integrals.
 */
public class MCMoveClusterMolecule extends MCMoveMolecule {
    
    private static final long serialVersionUID = 1L;
    private final MeterClusterWeight weightMeter;

    public MCMoveClusterMolecule(Simulation sim, PotentialMaster potentialMaster) {
    	this(potentialMaster,sim.getRandom(), 1.0);
    }
    
    public MCMoveClusterMolecule(PotentialMaster potentialMaster, IRandom random, double stepSize) {
        super(potentialMaster,random,stepSize,Double.POSITIVE_INFINITY,false);
        weightMeter = new MeterClusterWeight(potential);
        setName("MCMoveClusterMolecule");
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        weightMeter.setPhase(p);
    }
    
    public boolean doTrial() {
        if(phase.moleculeCount()==1) return false;
        
        atom = atomSource.getAtom();
        while (atom.getIndex() == 0) {
            atom = atomSource.getAtom();
        }
        
        uOld = weightMeter.getDataAsScalar();
        groupTranslationVector.setRandomCube(random);
        groupTranslationVector.TE(stepSize);
        moveMoleculeAction.actionPerformed(atom);
        uNew = Double.NaN;
//        System.out.println(((AtomTreeNodeGroup)((AtomTreeNodeGroup)phases[0].speciesMaster.node.childList.getFirst().node).childList.getLast().node).childList.getFirst().coord.position());
        ((PhaseCluster)phase).trialNotify();
        return true;
    }
    
    public double getB() {return 0.0;}
    
    public double getA() {
        uNew = weightMeter.getDataAsScalar();
//        if (Simulation.random.nextInt(200000) == 5) {
//            System.out.println("uOld "+uOld+" uNew "+uNew);
//        }
        return (uOld==0.0) ? Double.POSITIVE_INFINITY : uNew/uOld;
    }
    
    public void acceptNotify() {
        super.acceptNotify();
        ((PhaseCluster)phase).acceptNotify();
//        System.out.println(atom+" accepted => "+atom.type.getPositionDefinition().position(atom));
    }
    
    public void rejectNotify() {
        super.rejectNotify();
        ((PhaseCluster)phase).rejectNotify();
//        System.out.println(atom+" rejected => "+atom.type.getPositionDefinition().position(atom));
    }
        
}