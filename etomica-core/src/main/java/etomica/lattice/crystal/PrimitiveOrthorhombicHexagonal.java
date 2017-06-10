/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice.crystal;
import etomica.space.Vector;
import etomica.math.geometry.Polytope;
import etomica.math.geometry.Rectangle;
import etomica.space.Space;

/**
 * Primitive group for an orthorhombic system.  Primitive vectors orthogonal 
 * but with a length ratio to acommodate a hexagonal basis.
 * b = (2/sqrt(3))a ; alpha = 90deg.
 */
public class PrimitiveOrthorhombicHexagonal extends Primitive {
    
    public PrimitiveOrthorhombicHexagonal(Space space, double a) {
        super(space);
        //set up orthogonal vectors of unit size
        setSize(new double[]{a, Math.sqrt(3)*a});
        setAngles(new double[]{rightAngle});
    }
    
    public Primitive makeReciprocal() {
        return new PrimitiveOrthorhombicHexagonal(space, 2.0*Math.PI/size[0]);
    }
    
    /**
     * Sets A to the given value and B to 2.0/sqrt(3.0)*A
     */
    public void setSizeA(double newA) {
        if (size[0] == newA) {
            return;
        }
        double[] newSize = new double[D];
        newSize[0] = newA;
        newSize[1] = Math.sqrt(3)*newA;
        setSize(newSize);
    }
    public double getSizeA() {return size[0];}
        
    /**
     * Returns a new, identical instance of this primitive.
     */
    public Primitive copy() {
        return new PrimitiveOrthorhombicHexagonal(space, size[0]);
    }
    
    protected void update() {
        for (int i=0; i<D; i++) {
            latticeVectors[i].setX(i,size[i]);
        }
    }
    
    public void scaleSize(double scale) {
        double[] newSize = new double[D];
        for (int i=0; i<D; i++) {
            newSize[i] = scale*size[i];
        }
        setSize(newSize);
    }        
    
    public int[] latticeIndex(Vector q) {
        for(int i=0; i<D; i++) {
            double x = q.getX(i)/size[i];
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
    }

    public int[] latticeIndex(Vector q, int[] dimensions) {
        for(int i=0; i<D; i++) {
            double x = q.getX(i)/size[i];
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return idx;
    }
    
    public Polytope wignerSeitzCell() {
        throw new RuntimeException("method PrimitiveOrthorhombic.wignerSeitzCell not yet implemented");
    }
    
    /**
     * Returns a new Cuboid with edges of length given by the current
     * values of the primitive vectors.
     */
    public Polytope unitCell() {
        return new Rectangle(space, size[0], size[1]);
    }
    
    public String toString() {return "OrthorhombicHexagonal";}

    private static final long serialVersionUID = 1L;
}
    
