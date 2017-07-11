/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

/**
 * Potential in which attaches a harmonic spring between each affected atom and
 * the nearest boundary in each direction.
 *
 * This class has not been used or checked for correctness.
 *
 * @author David Kofke
 */
 
public class P1Harmonic extends Potential1 implements PotentialSoft {
    
    private static final long serialVersionUID = 1L;
    private double w = 100.0;
    private final Vector[] force;
    private final Vector x0;
    
    public P1Harmonic(Space space) {
        super(space);
        force = new Vector[]{space.makeVector()};
        x0 = space.makeVector();
    }
    public void setSpringConstant(double springConstant) {
        w = springConstant;
    }
    
    public double getSpringConstant() {
        return w;
    }
    
    public void setX0(Vector x0) {
        this.x0.E(x0);
    }
    
    public Vector getX0() {
        return x0;
    }
    
    public Dimension getX0Dimension() {
        return Length.DIMENSION;
    }

    public Dimension getSpringConstantDimension() {
        return new CompoundDimension(new Dimension[]{Energy.DIMENSION,Length.DIMENSION},new double[]{1,-2});
    }

    public double energy(IAtomList a) {
        return 0.5*w*a.getAtom(0).getPosition().Mv1Squared(x0);
    }
    
    public double virial(IAtomList a) {
        return 0.0;
    }

    public Vector[] gradient(IAtomList a){
        Vector r = a.getAtom(0).getPosition();
        force[0].Ev1Mv2(r,x0);
        force[0].TE(w);
            
        return force;
    }
        
    public Vector[] gradient(IAtomList a, Tensor pressureTensor){
        return gradient(a);
    }
}
   
