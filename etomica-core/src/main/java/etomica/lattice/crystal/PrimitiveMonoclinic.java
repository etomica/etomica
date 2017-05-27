/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice.crystal;
import etomica.space.Vector;
import etomica.math.geometry.Polytope;
import etomica.space.Space;
import etomica.space3d.Space3D;

/**
 * Primitive group for a monoclinic system.  One primitive-vector
 * angle not normal, and vectors not necessarily of equal length.
 * a != b != c; alpha = gamma = 90deg, beta >= 90deg
 */
public class PrimitiveMonoclinic extends Primitive {
    
    private static final long serialVersionUID = 1L;

    public PrimitiveMonoclinic(Space space) {
        this(space, 1.0, 1.0, 1.0, rightAngle);
    }
    public PrimitiveMonoclinic(Space space, double a, double b, double c, double beta) {
        super(space);
        setSize(new double[]{a, b, c});//also sets reciprocal via update
        setAngleBeta(beta);
    }
    
    public Primitive makeReciprocal() {
        return new PrimitiveMonoclinicReciprocal(space, 2.0*Math.PI/(size[0]*Math.sin(angle[1])),
                2.0*Math.PI/size[1], 2.0*Math.PI/(size[2]*Math.sin(angle[1])), angle[1]);
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

    //direct lattice (ix = 0, iz = 2)
    // v[0] = (1,0,0); v[1] = (0,1,0); v[2] = (c,0,s)  (times a, b, c)
     
    //reciprocal lattice (ix = 2, iz = 0)
    // v[0] = (s,0,-c); v[1] = (0,1,0); v[2] = (0,0,1);  (times a, b, c)
    protected void update() {
        latticeVectors[0].setX(0,size[0]);
        latticeVectors[1].setX(1,size[1]);
        latticeVectors[2].setX(0,size[2]*Math.cos(angle[1]));
        latticeVectors[2].setX(2,size[2]*Math.sin(angle[1]));
    }
    
    public void setAngleBeta(double t) {
        setAngles(new double[]{rightAngle, t, rightAngle});
    }
    
    public double getAngleBeta() {
        return angle[1];
    }

    public double getAngleAlpha() {
        return angle[0];
    }

    public double getAngleGamma() {
        return angle[2];
    }

    /**
     * Returns a new, identical instance of this primitive.
     */
    public Primitive copy() {
        return new PrimitiveMonoclinic(space, size[0], size[1], size[2], angle[1]);
    }
        
    public void scaleSize(double scale) {
        setSize(new double[]{size[0]*scale, size[1]*scale, size[2]*scale});
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
    
    public String toString() {return "Monoclinic";}

    protected static class PrimitiveMonoclinicReciprocal extends PrimitiveMonoclinic {
        private static final long serialVersionUID = 1L;

        public PrimitiveMonoclinicReciprocal(Space space, double a, double b, double c, double beta) {
            super(space, a, b, c, beta);
        }

        protected void update() {
            latticeVectors[0].setX(0, size[0]*Math.sin(angle[1]));
            latticeVectors[0].setX(2,-size[0]*Math.cos(angle[1]));
            latticeVectors[1].setX(1, size[1]);
            latticeVectors[2].setX(2, size[2]);
        }
    }

    public static void main(String args[]) {
        PrimitiveMonoclinic primitive = new PrimitiveMonoclinic(Space3D.getInstance(), 1, 1, 1, Math.PI*100/180);
        Vector[] v = primitive.vectors();
        Primitive reciprocal = primitive.makeReciprocal();
        Vector[] vr = reciprocal.vectors();
        for (int i=0; i<v.length; i++) {
            for (int j=0; j<vr.length; j++) {
                System.out.println(i+" "+j+" "+v[i].dot(vr[j]));
            }
        }
        System.out.println(v[0]);
        System.out.println(vr[2]);
    }
}
