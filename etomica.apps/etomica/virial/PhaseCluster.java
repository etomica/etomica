package etomica.virial;

import etomica.*;
import etomica.atom.AtomList;

/**
 * @author kofke
 *
 * Extension of Phase that forms and holds a PairSet instance for all of the
 * atoms in the phase.  Also instantiates phase with a NONE boundary type.
 */
public class PhaseCluster extends Phase {

	/**
	 * Constructor for PhaseCluster.
	 */
	public PhaseCluster() {
		this(Simulation.instance);
	}

	/**
	 * Constructor for PhaseCluster.
	 * @param parent
	 */
	public PhaseCluster(SimulationElement parent) {
		super(parent);
		setBoundary(space.makeBoundary(etomica.space3d.Boundary.NONE));	
	}
	
	public PairSet getPairSet() {
		if(pairSet == null && speciesMaster.atomList.size() > 0) {
//			pairSet = new PairSet(speciesMaster.atomList);
			pairSet = new PairSet(new AtomList(makeMoleculeIterator()));
		}
		return pairSet;
	}
	public void setPairSet(PairSet pairs) {
		this.pairSet = pairs;
	}

	private PairSet pairSet;
}
