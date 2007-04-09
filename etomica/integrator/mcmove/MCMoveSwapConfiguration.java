package etomica.integrator.mcmove;

import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorPhase;
import etomica.integrator.IntegratorPT.MCMoveSwap;
import etomica.integrator.IntegratorPT.MCMoveSwapFactory;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.space.IVector;


/**
 * Basic MCMove for swapping coordinates of atoms in two phases.
 * Requires same number of atoms in each phase.
 */
public class MCMoveSwapConfiguration extends MCMove implements MCMoveSwap {

    private static final long serialVersionUID = 1L;
	private final IntegratorPhase integrator1, integrator2;	
	private final AtomIteratorLeafAtoms iterator1 = new AtomIteratorLeafAtoms();
	private final AtomIteratorLeafAtoms iterator2 = new AtomIteratorLeafAtoms();
	private final AtomIteratorLeafAtoms affectedAtomIterator = new AtomIteratorLeafAtoms();
	private final IVector r;
	private double u1, u2, temp1, temp2, deltaU1;
	private final Phase[] swappedPhases = new Phase[2];

	public MCMoveSwapConfiguration(PotentialMaster potentialMaster, 
	                                IntegratorPhase integrator1, IntegratorPhase integrator2) {
  		super(potentialMaster);
		r = potentialMaster.getSpace().makeVector();
		this.integrator1 = integrator1;
		this.integrator2 = integrator2;
	}

	public boolean doTrial() {
		temp1 = integrator1.getTemperature();
		temp2 = integrator2.getTemperature();

        u1 = integrator1.getPotentialEnergy();
        u2 = integrator2.getPotentialEnergy();
        deltaU1 = Double.NaN;
        return true;
    }
    
    public double getA() {
    	// have to do this here since Integrator won't understand T dependence 
        deltaU1 = u2 - u1;  //if accepted, energy of phase1 will be u2, and its old energy is u1
        return Math.exp(-deltaU1*((1/temp1) - (1/temp2)));
    }
    
    public double getB() {
        return 0.0;
    }
	
	/**
	 * Swaps positions of molecules in two phases.
     * 
     * @throws RuntimeException wrapping a ConfigurationOverlapException if overlap is detected in either phase
	 */
	public void acceptNotify() {
		iterator1.setPhase(integrator1.getPhase());
		iterator2.setPhase(integrator2.getPhase());

		iterator1.reset();
		iterator2.reset();

		while(iterator1.hasNext()) {
			AtomLeaf a1 = (AtomLeaf)iterator1.nextAtom();
			AtomLeaf a2 = (AtomLeaf)iterator2.nextAtom();

			r.E(a1.getCoord().getPosition());
				
			a1.getCoord().getPosition().E(a2.getCoord().getPosition());
			a2.getCoord().getPosition().E(r);
		}
        ConfigurationOverlapException overlapException = null;
        try {
            //XXX grossly inefficient
            integrator1.reset();
        } catch(ConfigurationOverlapException e) {
            overlapException = e;
        }
        try {
            //XXX grossly inefficient
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
	    swappedPhases[0] = integrator1.getPhase();
	    swappedPhases[1] = integrator2.getPhase();
	    return swappedPhases;
	}

	public double energyChange(Phase phase) {
	    if(phase == integrator1.getPhase()) return +deltaU1;
	    if(phase == integrator2.getPhase()) return -deltaU1;
	    return 0.0;
	}
	
	public AtomIterator affectedAtoms(Phase p) {
	    if(p == integrator1.getPhase() || p == integrator2.getPhase()) {
	        affectedAtomIterator.setPhase(p);
	        affectedAtomIterator.reset();
	        return affectedAtomIterator;
	    }
	    return AtomIteratorNull.INSTANCE;
	}
    
    public final static SwapFactory FACTORY = new SwapFactory();
    
	protected static class SwapFactory implements MCMoveSwapFactory, java.io.Serializable {
	    public MCMove makeMCMoveSwap(PotentialMaster potentialMaster, 
                                     IntegratorPhase integrator1, IntegratorPhase integrator2) {
	        return new MCMoveSwapConfiguration(potentialMaster, integrator1, integrator2);
	    }
        private static final long serialVersionUID = 1L;
	} 
}