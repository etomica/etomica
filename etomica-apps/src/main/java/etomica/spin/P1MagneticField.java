/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin;

import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.potential.Potential1;
import etomica.space.Space;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */
public class P1MagneticField extends Potential1 {

    /**
     * @param space
     */
    public P1MagneticField(Space space) {
        super(space);
        direction = space.makeVector();
        direction.E(0.0);
        direction.setX(0,1.0);
    }

    /* (non-Javadoc)
     * @see etomica.Potential#energy(etomica.AtomSet)
     */
    public double energy(IAtomList atoms) {
        Vector r = atoms.getAtom(0).getPosition();
        return h * r.dot(direction);
    }
    
    
    /**
     * @return Returns the direction.
     */
    public Vector getDirection() {
        return direction;
    }
    /**
     * @param direction The direction to set.
     */
    public void setDirection(Vector direction) {
        this.direction.E(direction);
        this.direction.normalize();
    }
    /**
     * @return Returns the h.
     */
    public double getH() {
        return h;
    }
    /**
     * @param h The h to set.
     */
    public void setH(double h) {
        this.h = h;
    }

    private static final long serialVersionUID = 1L;
    private double h;
    private final Vector direction;
}
