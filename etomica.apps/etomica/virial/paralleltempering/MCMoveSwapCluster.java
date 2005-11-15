package etomica.virial.paralleltempering;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.integrator.IntegratorPT;
import etomica.integrator.IntegratorPhase;
import etomica.integrator.MCMove;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.virial.IntegratorClusterMC;
import etomica.virial.MCMoveCluster;
import etomica.virial.PhaseCluster;

/**
 * Swaps configurations and pairSet between phases for a virial clustering simulation. 
 */
public class MCMoveSwapCluster extends MCMove implements MCMoveCluster, IntegratorPT.MCMoveSwap {

    private IntegratorClusterMC integrator1, integrator2;	
    private AtomIteratorLeafAtoms iterator1 = new AtomIteratorLeafAtoms();
    private AtomIteratorLeafAtoms iterator2 = new AtomIteratorLeafAtoms();
    private AtomIteratorAllMolecules affectedAtomIterator = new AtomIteratorAllMolecules();
    private Vector r;
    private PhaseCluster phase1, phase2;
    private double weightOld1, weightOld2;
    private double weightNew1, weightNew2;
    private final Phase[] swappedPhases = new Phase[2];

    public MCMoveSwapCluster(PotentialMaster potentialMaster, 
                             IntegratorClusterMC integrator1, IntegratorClusterMC integrator2) {
        super(potentialMaster,2);
        r = potentialMaster.getSpace().makeVector();
        setTunable(false);
        this.integrator1 = integrator1;
        this.integrator2 = integrator2;		
    }

    public boolean doTrial() {
        if(phase1 == null || phase2 == null) {
            phase1 = (PhaseCluster)integrator1.getPhase();
            phase2 = (PhaseCluster)integrator2.getPhase();
            iterator1.setPhase(phase1);
            iterator2.setPhase(phase2);
        }

        weightOld1 = phase1.getSampleCluster().value(phase1.getCPairSet(), phase1.getAPairSet());
        weightOld2 = phase2.getSampleCluster().value(phase2.getCPairSet(), phase2.getAPairSet());
        
//        System.out.println("in trial "+integrator2.getWeight()+" "+weightOld2);
//        System.out.println("in trial "+integrator1.getWeight()+" "+weightOld1);
        iterator1.reset();
        iterator2.reset();

        while(iterator1.hasNext()) {
            Atom a1 = iterator1.nextAtom();
            Atom a2 = iterator2.nextAtom();

            //swap coordinates
            r.E(a1.coord.position());
            
            a1.coord.position().E(a2.coord.position());
            a2.coord.position().E(r);
        }

        //assumes energy will be determined using only pairSets in phases
        phase1.trialNotify(null);
        phase2.trialNotify(null);
		
        weightNew1 = weightNew2 = Double.NaN;
        return true;
    }
    
    public double lnTrialRatio() {return 0.0;}
    
    public double trialRatio() {return 1.0;}
    
    public double lnProbabilityRatio() {
        return Math.log(probabilityRatio());
    }
    
    public double probabilityRatio() {
        weightNew1 = phase1.getSampleCluster().value(phase1.getCPairSet(), phase1.getAPairSet());
        weightNew2 = phase2.getSampleCluster().value(phase2.getCPairSet(), phase2.getAPairSet());
//        System.out.println(weightOld1+" "+weightOld2+" "+weightNew1+" "+weightNew2);
        return  (weightNew1 * weightNew2) / (weightOld1 * weightOld2);
    }
	
    /**
     * Swaps positions of molecules in two phases.
     */
    public void acceptNotify() {
//        System.out.println("accepted");
		
        phase1.acceptNotify();
        phase2.acceptNotify();
    }
	
    public void rejectNotify() {
//        System.out.println("rejected");
        iterator1.reset();
        iterator2.reset();

        while(iterator1.hasNext()) {
            Atom a1 = iterator1.nextAtom();
            Atom a2 = iterator2.nextAtom();

            //swap coordinates
            r.E(a1.coord.position());
            
            a1.coord.position().E(a2.coord.position());
            a2.coord.position().E(r);
        }

        phase1.rejectNotify();
        phase2.rejectNotify();
    }
    
    public double energyChange(Phase phase) {
        if(phase == phase1) return weightNew1/weightOld1;
        if(phase == phase2) return weightNew2/weightOld2;
        return 0.0;
    }

    /**
     * Implementation of MCMoveSwap interface
     */
    public Phase[] swappedPhases() {
        swappedPhases[0] = phase1;
        swappedPhases[1] = phase2;
        return swappedPhases;
    }

    public AtomIterator affectedAtoms(Phase p) {
        if(p == phase1 || p == phase2) {
            affectedAtomIterator.setPhase(p);
            affectedAtomIterator.reset();
            return affectedAtomIterator;
        }
		    return AtomIterator.NULL;
    }
	
    public static SwapFactory FACTORY = new SwapFactory();
    
    private static class SwapFactory implements IntegratorPT.MCMoveSwapFactory, java.io.Serializable {
        public MCMove makeMCMoveSwap(PotentialMaster potentialMaster, 
                                     IntegratorPhase integrator1, IntegratorPhase integrator2) {
            return new MCMoveSwapCluster(potentialMaster, 
                                         (IntegratorClusterMC)integrator1, (IntegratorClusterMC)integrator2);
        }
    } 
	
}

