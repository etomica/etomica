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
public class MCMoveClusterMolecule extends MCMoveMolecule implements MCMoveCluster {
    
    private final MeterClusterWeight weightMeter;

    public MCMoveClusterMolecule(Simulation sim) {
    	this(sim.potentialMaster,sim.getDefaults().atomSize);
    }
    
    public MCMoveClusterMolecule(PotentialMaster potentialMaster, double stepSize) {
        super(potentialMaster,stepSize,Double.POSITIVE_INFINITY,false);
        weightMeter = new MeterClusterWeight(potential);
        setName("MCMoveClusterMolecule");
    }
    
    public void setPhase(Phase[] p) {
        super.setPhase(p);
        weightMeter.setPhase(p[0]);
    }
    
    public boolean doTrial() {
//        System.out.println(((AtomTreeNodeGroup)((AtomTreeNodeGroup)phases[0].speciesMaster.node.childList.getFirst().node).childList.getLast().node).childList.getFirst().coord.position());
        Phase phase = phases[0];
        if(phase.moleculeCount()==1) return false;
        
        atom = phase.randomMolecule();
        while (atom.node.getOrdinal() == 1) {
            atom = phase.randomMolecule();
        }
        
        energyMeter.setTarget(atom);
        uOld = weightMeter.getDataAsScalar();
        groupTranslationVector.setRandomCube();
        groupTranslationVector.TE(stepSize);
        moveMoleculeAction.actionPerformed(atom);
        uNew = Double.NaN;
//        System.out.println(((AtomTreeNodeGroup)((AtomTreeNodeGroup)phases[0].speciesMaster.node.childList.getFirst().node).childList.getLast().node).childList.getFirst().coord.position());
        ((PhaseCluster)phases[0]).trialNotify(atom);
        return true;
    }
    
    public double trialRatio() {return 1.0;}
    
    public double lnProbabilityRatio() {
        return Math.log(probabilityRatio());
    }
    
    public double probabilityRatio() {
        uNew = weightMeter.getDataAsScalar();
//        if (Simulation.random.nextInt(200000) == 5) {
//            System.out.println("uOld "+uOld+" uNew "+uNew);
//        }
        return (uOld==0.0) ? Double.POSITIVE_INFINITY : uNew/uOld;
    }
    
    public void acceptNotify() {
        super.acceptNotify();
        ((PhaseCluster)phases[0]).acceptNotify();
    }
    
    public void rejectNotify() {
        super.rejectNotify();
        ((PhaseCluster)phases[0]).rejectNotify();
    }
        
}