/*
 * Created on Jan 16, 2004
 */
package etomica.chem;
import etomica.AtomFactory;
import etomica.Default;
import etomica.Potential;
import etomica.Space;
import etomica.potential.PotentialTruncation;
import etomica.potential.PotentialTruncationSimple;

/**
 * Top-level class for a molecular model.
 */
public abstract class Model {
	
	public Model() {
		this(Default.TRUNCATE_POTENTIALS ? 
						new PotentialTruncationSimple()
					  : PotentialTruncation.NULL);
	}
	public Model(PotentialTruncation truncation) {
		this.truncation = truncation;
	}
	
	public abstract AtomFactory makeAtomFactory(Space space);
	
	public abstract Potential makePotential(Space space);
	
	/**
	 * NeighborIteration indicates if atoms at this hierarchy in the model
	 * should be subject to neighbor-based iteration.  Default is such that
	 * neighbor iteration is performed at highest level of hierarchy.
	 * @return boolean true if this level is subject to neighbor iteration 
	 */
	public boolean doNeighborIteration() {return doNeighborIteration;}
	public void setDoNeighborIteration(boolean b) {doNeighborIteration = b;}
	
	private boolean doNeighborIteration = false;
	protected PotentialTruncation truncation;

	/**
	 * Returns the truncation.
	 * @return PotentialTruncation
	 */
	public PotentialTruncation getTruncation() {
		return truncation;
	}

	/**
	 * Sets the truncation.
	 * @param truncation The truncation to set
	 */
	public void setTruncation(PotentialTruncation truncation) {
		this.truncation = truncation;
	}

}
