/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice.crystal;
import etomica.space.Vector;
import etomica.math.geometry.Polytope;
import etomica.space.Space;

/**
 * Primitive group for a trigonal system.  All primitive vectors have the same 
 * length, but there are no restrictions on primitive-vector angles.
 * a = b = c; alpha != gamma != beta
 */
public class PrimitiveTrigonal extends Primitive {

    private static final long serialVersionUID = 1L;

    public PrimitiveTrigonal(Space space) {
        this(space, 1.0, rightAngle, rightAngle, rightAngle);
    }
    public PrimitiveTrigonal(Space space, double cubicSize,
                             double alpha, double beta, double gamma) {
        super(space);
        setCubicSize(cubicSize);
        setAngles(new double[]{alpha, beta, gamma});
    }

    public Primitive makeReciprocal() {
        throw new RuntimeException("I don't actually know how to make a reciprocal primitive");
    }
    
    /**
     * Sets the length of all primitive vectors to the given value.
     */
    public void setCubicSize(double newCubicSize) {
        double[] sizeArray = new double[D];
        for(int i=0; i<D; i++) {
            sizeArray[i] = newCubicSize;
        }
        setSize(sizeArray);
    }
    
    protected void update() {
        double cosAlpha = Math.cos(angle[0]);
        double cosBeta = Math.cos(angle[1]);
        double cosGamma = Math.cos(angle[2]);
        double sinGamma = Math.sin(angle[2]);
        latticeVectors[0].setX(0,size[0]);
        latticeVectors[1].setX(0,size[0]*cosGamma);
        latticeVectors[1].setX(1,size[0]*sinGamma);
        latticeVectors[2].setX(0,size[0]*cosBeta);
        latticeVectors[2].setX(1,size[0]*(cosAlpha-cosBeta*cosGamma)/sinGamma);
        latticeVectors[2].setX(2,size[0]*Math.sqrt(1.0-cosAlpha*cosAlpha-cosBeta*cosBeta-cosGamma*cosGamma+2*cosAlpha*cosBeta*cosGamma)/sinGamma);
    }
    
    public void setAlpha(double t) {
        if (t == angle[0]) {
            return;
        }
        setAngles(new double[]{t, angle[1], angle[2]});
    }
    public double getAlpha() {return angle[0];}
    
    public void setBeta(double t) {
        if (t == angle[1]) {
            return;
        }
        setAngles(new double[]{angle[0], t, angle[2]});
    }
    public double getBeta() {return angle[1];}
    
    public void setGamma(double t) {
        if (t == angle[2]) {
            return;
        }
        setAngles(new double[]{angle[0], angle[1], t});
    }
    public double getGamma() {return angle[2];}
    
    /**
     * Returns a new, identical instance of this primitive.
     */
    public Primitive copy() {
        return new PrimitiveTrigonal(space, size[0], angle[0], angle[1], angle[2]);
    }
        
    public void scaleSize(double scale) {
        setCubicSize(size[0] * scale);
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
    
    public String toString() {return "Trigonal";}

}
