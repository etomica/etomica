/*
 * History
 * Created on Dec 17, 2004 by kofke
 */
package etomica.lattice;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */
public interface SiteFactory {

    public void makeSites(AbstractLattice lattice, Object[] sites);
}
