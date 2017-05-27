/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice.crystal;
import etomica.space.Vector;
import etomica.math.geometry.Polytope;
import etomica.space.Space;
import etomica.space3d.Space3D;

/**
 * Primitive group for a triclinic system.  No restrictions on
 * primitive-vector angles or lengths.
 * a != b != c; alpha != gamma != beta
 */
public class PrimitiveTriclinic extends Primitive {
    
    private static final long serialVersionUID = 1L;

    public PrimitiveTriclinic(Space space) {
        this(space, 1.0, 1.0, 1.0, rightAngle, rightAngle, rightAngle);
    }
    public PrimitiveTriclinic(Space space, double a, double b, double c,
                              double alpha, double beta, double gamma) {
        super(space);
        setSize(new double[]{a, b, c});
        setAngles(new double[]{alpha, beta, gamma});
    }

    //called by superclass constructor
    public Primitive makeReciprocal() {
        Vector aStar = space.makeVector();
        Vector bStar = space.makeVector();
        Vector cStar = space.makeVector();
        aStar.E(latticeVectors[1]);
        aStar.XE(latticeVectors[2]);
        double factor = 2.0*Math.PI/latticeVectors[0].dot(aStar); // a . (b X c)
        aStar.TE(factor);
        bStar.E(latticeVectors[2]);
        bStar.XE(latticeVectors[0]);
        factor = 2.0*Math.PI/latticeVectors[1].dot(bStar);
        bStar.TE(factor);
        cStar.E(latticeVectors[0]);
        cStar.XE(latticeVectors[1]);
        factor = 2.0*Math.PI/latticeVectors[2].dot(cStar);
        cStar.TE(factor);
        return new PrimitiveGeneral(space, new Vector[]{aStar, bStar, cStar});
    }
    
    public void setSizeA(double newA) {
        if (size[0] == newA) {
            return;
        }
        setSize(new double[]{newA, size[1], size[2]});
    }
    public double getSizeA() {return size[0];}
    
    public void setSizeB(double newB) {
        if (size[1] == newB) {
            return;
        }
        setSize(new double[]{size[0], newB, size[2]});
    }
    public double getSizeB() {return size[1];}
        
    public void setSizeC(double newC) {
        if (size[2] == newC) {
            return;
        }
        setSize(new double[]{size[0], size[1], newC});
    }
    public double getSizeC() {return size[2];}

    protected void update() {
        double cosAlpha = Math.cos(angle[0]);
        double cosBeta = Math.cos(angle[1]);
        double cosGamma = Math.cos(angle[2]);
        double sinGamma = Math.sin(angle[2]);
        latticeVectors[0].setX(0,size[0]);
        latticeVectors[1].setX(0,size[1]*cosGamma);
        latticeVectors[1].setX(1,size[1]*sinGamma);
        latticeVectors[2].setX(0,size[2]*cosBeta);
        latticeVectors[2].setX(1,size[2]*(cosAlpha-cosBeta*cosGamma)/sinGamma);
        latticeVectors[2].setX(2,size[2]*Math.sqrt(1.0-cosAlpha*cosAlpha-cosBeta*cosBeta-cosGamma*cosGamma+2*cosAlpha*cosBeta*cosGamma)/sinGamma);
    }
    
    public void setAngleAlpha(double t) {
        if (t == angle[0]) {
            return;
        }
        setAngles(new double[]{t, angle[1], angle[2]});
    }
    public double getAngleAlpha() {return angle[0];}
    
    public void setAngleBeta(double t) {
        if (t == angle[1]) {
            return;
        }
        setAngles(new double[]{angle[0], t, angle[2]});
    }
    public double getAngleBeta() {return angle[1];}
    
    public void setAngleGamma(double t) {
        if (t == angle[2]) {
            return;
        }
        setAngles(new double[]{angle[0], angle[1], t});
    }
    public double getAngleGamma() {return angle[2];}
    
    /**
     * Returns a new, identical instance of this primitive.
     */
    public Primitive copy() {
        return new PrimitiveTriclinic(space, size[0], size[1], size[2], angle[0], angle[1], angle[2]);
    }
        
    public void scaleSize(double scale) {
        setSize(new double[]{size[0]*scale, size[1]*scale, size[2]*scale});
    }        
    
    public int[] latticeIndex(Vector q) {
        throw new RuntimeException("nope");
    }

    public int[] latticeIndex(Vector q, int[] dimensions) {
        throw new RuntimeException("not this either");
    }
    
    public Polytope wignerSeitzCell() {
        throw new RuntimeException("method PrimitiveOrthorhombic.wignerSeitzCell not yet implemented");
    }
    
    public String toString() {return "Triclinic";}

    public static void main(String args[]) {
        PrimitiveTriclinic primitive = new PrimitiveTriclinic(Space3D.getInstance(), 1, 1.5, 2, Math.PI*0.4, Math.PI*0.45, Math.PI*0.6);
        Vector[] v = primitive.vectors();
        Primitive reciprocal = primitive.makeReciprocal();
        Vector[] vr = reciprocal.vectors();
        for (int i=0; i<v.length; i++) {
            for (int j=0; j<vr.length; j++) {
                System.out.println(i+" "+j+" "+v[i].dot(vr[j]));
            }
        }
    }
}
