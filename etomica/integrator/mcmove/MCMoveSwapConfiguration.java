package etomica.integrator.mcmove;

import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorPT.MCMoveSwap;
import etomica.integrator.IntegratorPT.MCMoveSwapFactory;
import etomica.box.Box;
import etomica.space.IVector;


/**
 * Basic MCMove for swapping coordinates of atoms in two boxs.
 * Requires same number of atoms in each box.
 */
public class MCMoveSwapConfiguration extends MCMove implements MCMoveSwap {

    private static final long serialVersionUID = 1L;
	private final IntegratorBox integrator1, integrator2;	
	private final AtomIteratorLeafAtoms affectedAtomIterator = new AtomIteratorLeafAtoms();
	private final IVector r;
	private double u1, u2, temp1, temp2, deltaU1;
	private final Box[] swappedBoxs = new Box[2];

	public MCMoveSwapConfiguration(IntegratorBox integrator1, IntegratorBox integrator2) {
  		super(null);
		r = integrator1.getBox().getSpace().makeVector();
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
        deltaU1 = u2 - u1;  //if accepted, energy of box1 will be u2, and its old energy is u1
        return Math.exp(-deltaU1*((1/temp1) - (1/temp2)));
    }
    
    public double getB() {
        return 0.0;
    }
	
	/**
	 * Swaps positions of molecules in two boxs.
     * 
     * @throws RuntimeException wrapping a ConfigurationOverlapException if overlap is detected in either box
	 */
	public void acceptNotify() {
        AtomSet leafList1 = integrator1.getBox().getSpeciesMaster().getLeafList();
        AtomSet leafList2 = integrator2.getBox().getSpeciesMaster().getLeafList();
        int nLeaf = leafList1.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomPositioned a1 = (IAtomPositioned)leafList1.getAtom(iLeaf);
            IAtomPositioned a2 = (IAtomPositioned)leafList2.getAtom(iLeaf);

			r.E(a1.getPosition());
				
			a1.getPosition().E(a2.getPosition());
			a2.getPosition().E(r);
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
	public Box[] swappedBoxs() {
	    swappedBoxs[0] = integrator1.getBox();
	    swappedBoxs[1] = integrator2.getBox();
	    return swappedBoxs;
	}

	public double energyChange(Box box) {
	    if(box == integrator1.getBox()) return +deltaU1;
	    if(box == integrator2.getBox()) return -deltaU1;
	    return 0.0;
	}
	
	public AtomIterator affectedAtoms(Box p) {
	    if(p == integrator1.getBox() || p == integrator2.getBox()) {
	        affectedAtomIterator.setBox(p);
	        affectedAtomIterator.reset();
	        return affectedAtomIterator;
	    }
	    return AtomIteratorNull.INSTANCE;
	}
    
    public final static SwapFactory FACTORY = new SwapFactory();
    
	protected static class SwapFactory implements MCMoveSwapFactory, java.io.Serializable {
	    public MCMove makeMCMoveSwap(IntegratorBox integrator1, IntegratorBox integrator2) {
	        return new MCMoveSwapConfiguration(integrator1, integrator2);
	    }
        private static final long serialVersionUID = 1L;
	} 
}