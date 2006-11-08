package etomica.models.hexane;
import etomica.lattice.crystal.Primitive;
import etomica.math.geometry.Polytope;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

/**
 * Implements the triclinic crystal structure derived from Dr. Monson's data.
 * @author nancycribbin
 */

public class PrimitiveHexane extends Primitive {
    
    public PrimitiveHexane(Space space) {
        this(space, true);
    }
    
    protected PrimitiveHexane(Space space, boolean makeReciprocal) {
        super(space, makeReciprocal);
        latticeVectors[0] = new Vector3D(-0.0450226624118149, 0.533363447960696, 2.05308369446092);
        latticeVectors[1] = new Vector3D(0.355734889084923, 1.281874771828840, 0.146059929673340);
        latticeVectors[2] = new Vector3D(0.763064234689232, -0.640942796742683, -0.083271977067129);

        double[] newSize = new double[3];
        newSize[0] = Math.sqrt(latticeVectors[0].squared());
        newSize[1] = Math.sqrt(latticeVectors[1].squared());
        newSize[2] = Math.sqrt(latticeVectors[2].squared());
        
        setSize(newSize);
        
        //In units of sigma
//        super(space, 2.86842545, 1.338313769, 1.290909603, 0.779541414, 1.109209538,
//                0.768992016);
        
        //In units of Angstroms
//        super(space, 28.68425449956680, 13.38313769270750, 12.90909602969220, 
//                0.779541414, 1.109209538, 0.768992016);
    }
    
    
    
    protected Primitive makeReciprocal() {
        // lattice vectors will be totally bogus.  we'll fix them later
        return new PrimitiveHexane(space, false);
    }
    
    //called by update method of superclass
    protected void updateReciprocal() {
        //XXX this does not update the reciprocal's size
        PrimitiveHexane recip = (PrimitiveHexane)reciprocal;
        Vector3D aStar = (Vector3D)recip.latticeVectors[0];
        Vector3D bStar = (Vector3D)recip.latticeVectors[1];
        Vector3D cStar = (Vector3D)recip.latticeVectors[2];
        Vector3D aVec = (Vector3D)latticeVectors[0];
        Vector3D bVec = (Vector3D)latticeVectors[1];
        Vector3D cVec = (Vector3D)latticeVectors[2];
        aStar.E(bVec);
        aStar.XE(cVec);
        double factor = 2.0*Math.PI/aVec.dot(aStar); // a . (b X c)
        aStar.TE(factor);
        bStar.E(cVec);
        bStar.XE(aVec);
        bStar.TE(factor);
        cStar.E(aVec);
        cStar.XE(bVec);
        cStar.TE(factor);
    }
    
    /**
     * Sets A and scales B and C to maintain relative size
     */
    public void setA(double newA) {
        scaleSize(newA/size[0]);
    }
    
    public void setSize(double[] newSize) {
        if (size[0] == newSize[0] && size[1] == newSize[1] && size[2] == newSize[2]) {
            // no change
            return;
        }
        if (size[0] == 0) {
            //initialization
            super.setSize(newSize);
        }

        // size can be scaled, but not the relative values cannot change
        double scale = newSize[0] / size[0];
        if (newSize[1] - size[1]*scale > 0.0000001*(newSize[1]+size[1]*scale) ||
            newSize[2] - size[2]*scale > 0.0000001*(newSize[1]+size[2]*scale)) {
            throw new IllegalArgumentException("You can't change the relative size of the primitive vectors");
        }
        double[] newNewSize = new double[3];
        newNewSize[0] = newSize[0];
        newNewSize[1] = scale*size[1];
        newNewSize[2] = scale*size[2];
        super.setSize(newNewSize);
    }
    
    protected void update() {
        super.update();
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
            double x = q.x(i)/size[i];
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
    }

    public int[] latticeIndex(Vector q, int[] dimensions) {
        for(int i=0; i<D; i++) {
            double x = q.x(i)/size[i];
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return idx;
    }
    
    public Polytope wignerSeitzCell() {
        throw new RuntimeException("method PrimitiveOrthorhombic.wignerSeitzCell not yet implemented");
    }
    
    public Polytope unitCell() {
        throw new RuntimeException("method PrimitiveOrthorhombic.unitCell not yet implemented");
    }
    
    public String toString() {return "Hexane";}

    private static final long serialVersionUID = 1L;
}
