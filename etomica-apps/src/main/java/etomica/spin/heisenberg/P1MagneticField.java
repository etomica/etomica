/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.heisenberg;

import etomica.api.IAtomList;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.ISpace;
import etomica.space.Tensor;


/**
 *
 * this class is to introduce electric field to the 2D Heisenberg Model.
 * Right now the electric filed is set to be zero since we are only interested
 * in the secondDerivative of free energy w.r.t electric filed when electric field
 * is zero
 *
 * @author David Kofke & Weisong Lin
 *
 */
public class P1MagneticField extends Potential1 implements PotentialSoft {

    /**
     * @param space use to define vector
     * @param dipoleMagnitude the dipole strength of dipole
     */
    public P1MagneticField(ISpace space, double dipoleMagnitude) {
        super(space);
        direction = space.makeVector();
        direction.E(0.0);
        direction.setX(0,1.0);
        
        this.dipoleMagnitude = dipoleMagnitude;
        dr = space.makeVector();
        dr.E(0);
        gradient = new IVectorMutable[1];
        gradient[0] = space.makeVector();
    }

    /**
     *
     * @param atoms atomlist in the system
     * @return energy of dipole in electric field
     */
    public double energy(IAtomList atoms) {
        IVectorMutable r = atoms.getAtom(0).getPosition();
        return dipoleMagnitude*h*r.dot(direction);
    }
    
    
    /**
     * @return Returns the direction.
     */
    public IVectorMutable getDirection() {
        return direction;
    }
    /**
     * @param direction The direction to set.
     */
    public void setDirection(IVectorMutable direction) {
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
    public void setH(double h) {//I didn't set up the electric fields
        this.h = h;
    }

    private static final long serialVersionUID = 1L;
    private double h;
    private final IVectorMutable direction;
    
    private double dipoleMagnitude;
    private IVectorMutable dr;
    private IVectorMutable [] gradient;

    /**
     *
     * @param atoms
     * @return 0 since no virial is used here
     */
	public double virial(IAtomList atoms) {
		return 0;
	}

    /**
     *
     * @param atoms
     * @return gradient vector
     */
	public IVector[] gradient(IAtomList atoms) {
		gradient[0].Ea1Tv1(-dipoleMagnitude*h, direction);
		return gradient;
	}

	public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
		return gradient(atoms);
	}
}
