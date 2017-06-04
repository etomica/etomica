/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.api.*;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.space.Space;

/**
 * Methods for a hard (impulsive), spherically-symmetric pair potential.
 * Subclasses must provide a concrete definition for the energy (method u(double)).
 */

public abstract class Potential2HardSpherical extends Potential2 implements PotentialHard, Potential2Spherical {
   
    public Potential2HardSpherical(Space space) {
	    super(space);
        dr = space.makeVector();
	}
	
	/**
    * The pair energy u(r^2) with no truncation applied.
    * @param r2 the square of the distance between the particles.
    */
    public abstract double u(double r2);

    /**
     * Energy of the pair as given by the u(double) method, with application
     * of any PotentialTruncation that may be defined for the potential.  This
     * does not take into account any false positioning that the Integrator may
     * be using.
     */
    public double energy(IAtomList pair) {
        IAtom atom0 = pair.getAtom(0);
        IAtom atom1 = pair.getAtom(1);

        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        boundary.nearestImage(dr);
        return u(dr.squared());
    }
    
    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    protected final Vector dr;
    protected IBoundary boundary;
}
