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
public interface NeighborCriterion extends AtomsetFilter {

	public boolean needUpdate(Atom atom);
	
	public void setPhase(Phase phase);
	
	public boolean unsafe(Atom atom);
}
