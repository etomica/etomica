/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.hexane;
import etomica.space.Vector;
import etomica.lattice.crystal.Primitive;
import etomica.math.geometry.Polytope;
import etomica.space.Space;
import etomica.space3d.Vector3D;

/**
 * Implements the triclinic crystal structure derived from Dr. Monson's data.
 * @author nancycribbin
 */

public class PrimitiveHexane extends Primitive {
    
    private static final long serialVersionUID = 1L;
    protected final Vector3D[] fixedLatticeVector = new Vector3D[]{
        new Vector3D(-0.0450226624118149, 0.533363447960696, 2.05308369446092),
        new Vector3D(0.355734889084923, 1.281874771828840, 0.146059929673340),
        new Vector3D(0.763064234689232, -0.640942796742683, -0.083271977067129)};
    protected final double[] fixedSizes;
    private Vector3D[] temp;
    
    public PrimitiveHexane(Space space) {
        super(space);

        fixedSizes = new double[3];
        for (int i=0; i<3; i++) {
            fixedSizes[i] = Math.sqrt(fixedLatticeVector[i].squared());
        }
        
        setSize(fixedSizes);
        
        //In units of sigma
//        super(space, 2.86842545, 1.338313769, 1.290909603, 0.779541414, 1.109209538,
//                0.768992016);
        
        //In units of Angstroms
//        super(space, 28.68425449956680, 13.38313769270750, 12.90909602969220, 
//                0.779541414, 1.109209538, 0.768992016);
    }
    
    
    
    public Primitive makeReciprocal() {
        // lattice vectors will be totally bogus.  we'll fix them later

        Vector3D temp = new Vector3D();
        temp.E(fixedLatticeVector[0]);
        temp.XE(fixedLatticeVector[1]);
        double consts = temp.dot(fixedLatticeVector[2]);
        consts = 2.0 * Math.PI / consts;
        
        PrimitiveHexane recip = new PrimitiveHexane(space);
        
        //construct the first fixed lattice vector of the reciprocal
        temp.E(fixedLatticeVector[1]);
        temp.XE(fixedLatticeVector[2]);
        temp.TE(consts);
        for(int i = 0; i < 3; i++){
            recip.fixedLatticeVector[0].setX(i, temp.getX(i));
        }
    
        //construct the second fixed lattice vector of the reciprocal
        temp.E(fixedLatticeVector[2]);
        temp.XE(fixedLatticeVector[0]);
        temp.TE(consts);
        for(int i=0; i < 3; i++){
            recip.fixedLatticeVector[1].setX(i, temp.getX(i));
        }
        
        
        //construct the third fixed lattice vector of the reciprocal
        temp.E(fixedLatticeVector[0]);
        temp.XE(fixedLatticeVector[1]);
        temp.TE(consts);
        for(int i = 0; i < 3; i++){
            recip.fixedLatticeVector[2].setX(i, temp.getX(i));
        }
        
        return recip;
    }

    /**
     * Sets A and scales B and C to maintain relative size
     */
    public void setA(double newA) {
        scaleSize(newA/size[0]);
    }
    
    protected void update() {
        double scale = size[0]/fixedSizes[0];
        for (int i=0; i<3; i++) {
            latticeVectors[i].Ea1Tv1(scale, fixedLatticeVector[i]);
        }
    }
    
    
    /**
     * Returns a new, identical instance of this primitive.
     */
    public Primitive copy() {
        return new PrimitiveHexane(space);
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
    
    public String toString() {return "Hexane";}

}
