/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.api.IElement;
import etomica.space.Space;
import etomica.space.Vector;


/**
 * Atom type for a sphere that has some feature depending upon an orientation coordinate.
 * For example an orientational dependent potential may be attached to an otherwise spherical atom
 */
public class AtomTypeOriented extends AtomType {

    protected final Vector momentOfIntertia;

    @Deprecated
    public AtomTypeOriented(IElement element, Space space) {
        super(element);
        momentOfIntertia = space.makeVector();
        momentOfIntertia.E(1);
    }

    public AtomTypeOriented(IElement element, Vector moment) {
        super(element);
        momentOfIntertia = moment;
    }

    /**
     * Returns the principal components of the moment of inertia of the
     * atom within the body-fixed frame. Do NOT modify the returned moment
     * of inertia returned.
     * @return the moment of inertia
     */
    public Vector getMomentOfInertia() {
        return momentOfIntertia;
    }

}
