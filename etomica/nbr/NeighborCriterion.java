/*
 * History
 * Created on Sep 22, 2004 by kofke
 */
package etomica.nbr;
import etomica.Atom;
import etomica.AtomsetFilter;
import etomica.Phase;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public abstract class NeighborCriterion implements AtomsetFilter {

	public NeighborCriterion() {}
	
	public abstract boolean needUpdate(Atom atom);
	
	public abstract void setPhase(Phase phase);
	
	public abstract boolean unsafe();
	
	public abstract void reset(Atom atom);
	
	//TODO consider ways to ensure this is removed from nbrmanager if no longer used
}
