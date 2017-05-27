/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice.crystal;

import etomica.space.Vector;
import etomica.math.geometry.Cube;
import etomica.math.geometry.Polytope;
import etomica.math.geometry.Square;
import etomica.space.Space;

/**
 * Primitive group for a cubic system.  All primitive
 * vectors orthogonal and of equal length.
 */
public class PrimitiveCubic extends Primitive {
    
    private static final long serialVersionUID = 1L;
    private double aBC;
    
    public PrimitiveCubic(Space space) {
        this(space, 1.0);
    }
    public PrimitiveCubic(Space space, double latticeConstant) {
        super(space);
        //set up orthogonal vectors of unit size
        setSizeABC(latticeConstant); //also sets reciprocal via update
        double[] newAngles = new double[D];
        for (int i=0; i<D; i++) {
            newAngles[i] = rightAngle;
        }
        setAngles(newAngles);
    }
    
    public Primitive makeReciprocal() {
        return new PrimitiveCubic(space, 2.0*Math.PI/aBC);
    }
    
    /**
     * Returns a new PrimitiveCubic with the same size as this one.
     */
    public Primitive copy() {
        return new PrimitiveCubic(space, aBC);
    }
    
    /**
     * Sets the length of all primitive vectors to the given value.
     */
    public void setSizeABC(double newCubicSize) {
        if (newCubicSize == aBC) {
            // no change
            return;
        }
        double[] sizeArray = new double[D];
        for(int i=0; i<D; i++) {
            sizeArray[i] = newCubicSize;
        }
        setSize(sizeArray);
        aBC = newCubicSize;
    }
    
    protected void update() {
        for(int i=0; i<D; i++) latticeVectors[i].setX(i, size[0]);
    }
    
    /**
     * Returns the common length of all primitive vectors.
     */
    public double getSizeABC() {return aBC;}

    public double getAngleAlpha() {
        return angle[0];
    }

    public double getAngleBeta() {
        return angle[1];
    }

    public double getAngleGamma() {
        return angle[2];
    }

    public void scaleSize(double scale) {
        setSizeABC(scale*aBC);
    }

    public int[] latticeIndex(Vector q) {
        for(int i=0; i<D; i++) {
            double x = q.getX(i)/aBC;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
    }
    
    public int[] latticeIndex(Vector q, int[] dimensions) {
        for(int i=0; i<D; i++) {
            double x = q.getX(i)/aBC;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return idx;
    }
    
    /**
     * Returns a new Square (if primitive is 2D) or Cube (if 3D) with edges
     * given by the size of the primitive vectors.
     */
    public Polytope wignerSeitzCell() {
        return (D == 2) ? new Square(space,aBC) : new Cube(space,aBC);
    }
    
    /**
     * Returns a new Square (if primitive is 2D) or Cube (if 3D) with edges
     * given by the size of the primitive vectors.
     */
    public Polytope unitCell() {
        return (D == 2) ? new Square(space,aBC) : new Cube(space,aBC);
    }
    
    public String toString() {return "Cubic";}
    
}
    
