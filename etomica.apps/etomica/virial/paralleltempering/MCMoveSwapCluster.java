package etomica.virial.paralleltempering;

import etomica.*;
import etomica.virial.*;

/**
 * Swaps configurations and pairSet between phases.
 */
public class MCMoveSwapCluster extends MCMove implements IntegratorPT.MCMoveSwap {

	private IntegratorMC integrator1, integrator2;	
	private final IteratorDirective iteratorDirective = new IteratorDirective();
	private AtomIteratorListSimple iterator1 = new AtomIteratorListSimple();
	private AtomIteratorListSimple iterator2 = new AtomIteratorListSimple();
	private AtomIteratorMolecule affectedAtomIterator = new AtomIteratorMolecule();
	private Space.Vector r;
	private PhaseCluster phase1, phase2;
	private double buOld1, buOld2, buNew1, buNew2 = Double.NaN;
	private double beta1, beta2;
	private final Phase[] swappedPhases = new Phase[2];

	public MCMoveSwapCluster(IntegratorMC integratorMC, 
									Integrator integrator1, Integrator integrator2) {
		super(integratorMC);
		r = integratorMC.space.makeVector();
		setTunable(false);
		this.integrator1 = (IntegratorMC)integrator1;
		this.integrator2 = (IntegratorMC)integrator2;		
	}

	public boolean doTrial() {
		if(phase1 == null || phase2 == null) {
			phase1 = (PhaseCluster)integrator1.getPhase(0); 
			phase2 = (PhaseCluster)integrator2.getPhase(0);
			iterator1.setBasis(phase1.speciesMaster.atomList);
			iterator2.setBasis(phase2.speciesMaster.atomList);
		}

		beta1 = 1.0/integrator1.temperature();
		beta2 = 1.0/integrator2.temperature();
		
		buOld1 = beta1*potential.calculate(phase1, iteratorDirective, energy.reset()).sum();
		buOld2 = beta2*potential.calculate(phase2, iteratorDirective, energy.reset()).sum();
		
		//assumes energy will be determined using only pairSets in phases
		swapPairSets(phase1, phase2);
		
		buNew1 = buNew2 = Double.NaN;
		return true;
	}
    
	public double lnTrialRatio() {return 0.0;}
    
	public double lnProbabilityRatio() {
		buNew1 = beta1*potential.calculate(phase1, iteratorDirective, energy.reset()).sum();
		buNew2 = beta2*potential.calculate(phase2, iteratorDirective, energy.reset()).sum();
		return  -(buNew1 + buNew2) + (buOld1 + buOld2);
	}
	
	private void swapPairSets(PhaseCluster phase1, PhaseCluster phase2) {
		PairSet setA = phase1.getPairSet();
		phase1.setPairSet(phase2.getPairSet());
		phase2.setPairSet(setA);
	}
	
	/**
	 * Swaps positions of molecules in two phases.
	 */
	public void acceptNotify() {
		//put pairSets back and now implement swap by exchanges positions
		swapPairSets(phase1, phase2);
		
		iterator1.reset();
		iterator2.reset();

		while(iterator1.hasNext()) {
			Atom a1 = iterator1.next();
			Atom a2 = iterator2.next();//assumes N1 == N2

			//swap coordinates
			r.E(a1.coord.position());
				
			a1.coord.translateTo(a2.coord.position());
			a2.coord.translateTo(r);			
		}
	}
	
	/**
	 * Just re-exchanges pairSets between phases.
	 */
	public void rejectNotify() {
		swapPairSets(phase1, phase2);
	}
	
	public double energyChange(Phase phase) {
		if(phase == phase1) return (buNew1 - buOld1)/beta1;
		else if(phase == phase2) return (buNew2 - buOld2)/beta2;
		else return 0.0;
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
			affectedAtomIterator.setBasis(p);
			affectedAtomIterator.reset();
			return affectedAtomIterator;
		}
		else return AtomIterator.NULL;
	}
	
public static class Factory implements IntegratorPT.MCMoveSwapFactory {
	public MCMove makeMCMoveSwap(IntegratorMC integratorMC, 
									Integrator integrator1, Integrator integrator2) {
		return new MCMoveSwapCluster(integratorMC, 
										integrator1, integrator2);
	}
}//end of MCMoveSwapFactoryDefault 
	
}//end of MCMoveSwapPolydisperse

