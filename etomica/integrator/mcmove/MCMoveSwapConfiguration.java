/*
 * Created on Jul 17, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.integrator.mcmove;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.Integrator;
import etomica.integrator.MCMove;
import etomica.integrator.IntegratorPT.MCMoveSwap;
import etomica.integrator.IntegratorPT.MCMoveSwapFactory;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;


/**
	 * Basic MCMove for swapping coordinates of atoms in two phases.
	 * Requires same number of atoms in each phase.
	 */
public class MCMoveSwapConfiguration extends MCMove implements MCMoveSwap {

	private final Integrator integrator1, integrator2;	
	private final AtomIteratorLeafAtoms iterator1 = new AtomIteratorLeafAtoms();
	private final AtomIteratorLeafAtoms iterator2 = new AtomIteratorLeafAtoms();
	private final AtomIteratorLeafAtoms affectedAtomIterator = new AtomIteratorLeafAtoms();
	private final Vector r;
	private double u1, u2, temp1, temp2, deltaU1;
	private final Phase[] swappedPhases = new Phase[2];

	public MCMoveSwapConfiguration(PotentialMaster potentialMaster, 
	                                Integrator integrator1, Integrator integrator2) {
  		super(potentialMaster,2);
		r = potentialMaster.getSpace().makeVector();
		setTunable(false);
		this.integrator1 = integrator1;
		this.integrator2 = integrator2;
	}

	public boolean doTrial() {
		temp1 = integrator1.getTemperature();
		temp2 = integrator2.getTemperature();

        u1 = integrator1.getPotentialEnergy()[0];
        u2 = integrator2.getPotentialEnergy()[0];
        deltaU1 = Double.NaN;
        return true;
    }
    
    // NOOP
    public void setPhase(Phase[] p) {}
    
    public double lnTrialRatio() {return 0.0;}
    
    public double lnProbabilityRatio() {
        deltaU1 = u2 - u1;  //if accepted, energy of phase1 will be u2, and its old energy is u1
		return  -deltaU1*((1/temp1) - (1/temp2));
	}
	
	/**
	 * Swaps positions of molecules in two phases.
     * 
     * @throws RuntimeException wrapping a ConfigurationOverlapException if overlap is detected in either phase
	 */
	public void acceptNotify() {
		iterator1.setPhase(integrator1.getPhase()[0]);
		iterator2.setPhase(integrator2.getPhase()[0]);

		iterator1.reset();
		iterator2.reset();

		while(iterator1.hasNext()) {
			Atom a1 = iterator1.nextAtom();
			Atom a2 = iterator2.nextAtom();

			r.E(a1.coord.position());
				
			a1.coord.position().E(a2.coord.position());
			a2.coord.position().E(r);
		}
        ConfigurationOverlapException overlapException = null;
        try {
            integrator1.reset();
        } catch(ConfigurationOverlapException e) {
            overlapException = e;
        }
        try {
            integrator2.reset();
        } catch(ConfigurationOverlapException e) {
            overlapException = e;
        }
        if(overlapException != null) {
            throw new RuntimeException(overlapException);
        }
	}
	
	/**
     * Performs no action; nothing required when move rejected.
     */
	public void rejectNotify() {}
	
	/**
	 * Implementation of MCMoveSwap interface
	 */
	public Phase[] swappedPhases() {
	    swappedPhases[0] = integrator1.getPhase()[0];
	    swappedPhases[1] = integrator2.getPhase()[0];
	    return swappedPhases;
	}

	public double energyChange(Phase phase) {
	    if(phase == integrator1.getPhase()[0]) return +deltaU1;
	    if(phase == integrator2.getPhase()[0]) return -deltaU1;
	    return 0.0;
	}
	
	public AtomIterator affectedAtoms(Phase p) {
	    if(p == integrator1.getPhase()[0] || p == integrator2.getPhase()[0]) {
	        affectedAtomIterator.setPhase(p);
	        affectedAtomIterator.reset();
	        return affectedAtomIterator;
	    }
	    return AtomIterator.NULL;
	}
    
    public transient static SwapFactory FACTORY = new SwapFactory();
    
	private static class SwapFactory implements MCMoveSwapFactory, java.io.Serializable {
	    public MCMove makeMCMoveSwap(PotentialMaster potentialMaster, 
                                     Integrator integrator1, Integrator integrator2) {
	        return new MCMoveSwapConfiguration(potentialMaster, integrator1, integrator2);
	    }
	} 
}