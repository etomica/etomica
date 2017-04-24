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
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */
public class P1MagneticField extends Potential1 implements PotentialSoft {

    /**
     * @param space
     */
    public P1MagneticField(ISpace space, double dipoleMagnitude) {//TODO I add dipoleMagnitude here 
        super(space);
        direction = space.makeVector();
        direction.E(0.0);
        direction.setX(0,1.0);//This one set the direction of atoms, should  I change it????
        
        //TODO
        this.dipoleMagnitude = dipoleMagnitude;
        dr = space.makeVector();
        dr.E(0);
        gradient = new IVectorMutable[1];
        gradient[0] = space.makeVector();
        //TODO
    }

    /* (non-Javadoc)
     * @see etomica.Potential#energy(etomica.AtomSet)
     */
    public double energy(IAtomList atoms) {
        IVectorMutable r = atoms.getAtom(0).getPosition();
        return dipoleMagnitude*h*r.dot(direction);//TODO Add mu here
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
    
    //TODO shoudl I add final
    private double dipoleMagnitude;
    private IVectorMutable dr;
    private IVectorMutable [] gradient;
    //TODO
	@Override
	public double virial(IAtomList atoms) {
		
		return 0;
	}

	@Override
	public IVector[] gradient(IAtomList atoms) {//TODO I thought I only need to return one vector not an array??
		gradient[0].Ea1Tv1(-dipoleMagnitude*h, direction);
		return gradient;
	}

	@Override
	public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
		return gradient(atoms);
	}
}
