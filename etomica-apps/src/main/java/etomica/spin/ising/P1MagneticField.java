/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.ising;

import etomica.atom.IAtom;
import etomica.atom.IAtomOriented;
import etomica.potential.IPotentialField;
import etomica.space.Space;
import etomica.space.Vector;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 */
public class P1MagneticField implements IPotentialField {

    private final Vector direction;
    private double h;

    /**
     * @param space
     */
    public P1MagneticField(Space space) {
        direction = space.makeVector();
        direction.E(0.0);
        direction.setX(0, 1.0);
    }

    /**
     * Returns the energy between IAtom atom and the field.
     */
    public double u(IAtom atom) {
        Vector r = ((IAtomOriented)atom).getOrientation().getDirection();
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
}
