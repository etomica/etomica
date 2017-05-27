/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice.crystal;
import etomica.space.Vector;
import etomica.math.geometry.Cuboid;
import etomica.math.geometry.Polytope;
import etomica.space.Space;

/**
 * Primitive group for an orthorhombic system.  All primitive
 * vectors orthogonal but not necessarily of equal length.
 * a != b != c; alpha = beta = gamma = 90deg.
 */
public class PrimitiveOrthorhombic extends Primitive {
    
    public PrimitiveOrthorhombic(Space space) {
        this(space, 1.0, 1.0, 1.0);
    }
    public PrimitiveOrthorhombic(Space space, double a, double b, double c) {
        super(space);
        //set up orthogonal vectors of unit size
        setSize(new double[]{a, b, c});
        setAngles(new double[]{rightAngle, rightAngle, rightAngle});
    }
    
    public Primitive makeReciprocal() {
        return new PrimitiveOrthorhombic(space, 2.0*Math.PI/size[0], 2.0*Math.PI/size[1], 2.0*Math.PI/size[2]);
    }
    
    public void setSizeA(double newA) {
        if (size[0] == newA) {
            return;
        }
        double[] newSize = new double[D];
        for (int i=0; i<D; i++) {
            newSize[i] = size[i];
        }
        newSize[0] = newA;
        setSize(newSize);
    }
    public double getSizeA() {return size[0];}
    
    public void setSizeB(double newB) {
        if (size[1] == newB) {
            return;
        }
        double[] newSize = new double[D];
        for (int i=0; i<D; i++) {
            newSize[i] = size[i];
        }
        newSize[1] = newB;
        setSize(newSize);
    }
    public double getSizeB() {return size[1];}
        
    public void setSizeC(double newC) {
        if (size[2] == newC) {
            return;
        }
        double[] newSize = new double[D];
        for (int i=0; i<D; i++) {
            newSize[i] = size[i];
        }
        newSize[2] = newC;
        setSize(newSize);
    }
    public double getSizeC() {return size[2];}

    public double getAngleAlpha() {
        return angle[0];
    }

    public double getAngleBeta() {
        return angle[1];
    }

    public double getAngleGamma() {
        return angle[2];
    }

    /**
     * Returns a new, identical instance of this primitive.
     */
    public Primitive copy() {
        return new PrimitiveOrthorhombic(space, size[0], size[1], size[2]);
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
        return new Cuboid(space, size[0], size[1], size[2]);
    }
    
    public String toString() {return "Orthorhombic";}

    private static final long serialVersionUID = 1L;
}
    
