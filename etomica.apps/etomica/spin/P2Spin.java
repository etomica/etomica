package etomica.spin;

import etomica.AtomPair;
import etomica.AtomSet;
import etomica.Space;
import etomica.potential.Potential2;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 22, 2005 by kofke
 */
public class P2Spin extends Potential2 {

    
    public P2Spin(Space space) {
        this(space, 1.0);
    }
    
    public P2Spin(Space space, double coupling) {
        super(space);
        setCoupling(coupling);
    }

    /* (non-Javadoc)
     * @see etomica.Potential#energy(etomica.AtomSet)
     */
    public double energy(AtomSet atoms) {
        AtomPair pair = (AtomPair)atoms;
        return -coupling * pair.atom0.coord.position().dot(pair.atom1.coord.position());
    }
    /* (non-Javadoc)
     * @see etomica.Potential#getRange()
     */
    public double getRange() {
        // TODO Auto-generated method stub
        return 0;
    }
    
    /**
     * @return Returns the coupling.
     */
    public double getCoupling() {
        return coupling;
    }
    /**
     * @param coupling The coupling to set.
     */
    public void setCoupling(double coupling) {
        this.coupling = coupling;
    }
    private double coupling;
}
