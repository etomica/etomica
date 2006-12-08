package etomica.virial;

import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;

/**
 * Standard Monte Carlo molecule-displacement trial move for cluster integrals.
 */
 
 /* History of changes
  * 7/9/02 Added energyChange() method
  */
public class MCMoveClusterMolecule extends MCMoveMolecule {
    
    private static final long serialVersionUID = 1L;
    private final MeterClusterWeight weightMeter;

    public MCMoveClusterMolecule(Simulation sim) {
    	this(sim.getPotentialMaster(),sim.getDefaults().atomSize);
    }
    
    public MCMoveClusterMolecule(PotentialMaster potentialMaster, double stepSize) {
        super(potentialMaster,stepSize,Double.POSITIVE_INFINITY,false);
        weightMeter = new MeterClusterWeight(potential);
        setName("MCMoveClusterMolecule");
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        weightMeter.setPhase(p);
    }
    
    public boolean doTrial() {
        if(phase.moleculeCount()==1) return false;
        
        atom = phase.randomMolecule();
        while (atom.getNode().getIndex() == 0) {
            atom = phase.randomMolecule();
        }
        
        uOld = weightMeter.getDataAsScalar();
        groupTranslationVector.setRandomCube();
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