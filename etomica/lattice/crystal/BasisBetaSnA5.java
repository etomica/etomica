package etomica.lattice.crystal;
import etomica.lattice.Basis;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

/**
 * Basis that helps construct a beta-tin crystal on a BravaisLattice.  
 * 
 * The beta-tin crystal structure has a tetragonal-shaped unit cell, but, unlike the 
 * simple tetragonal structure, it has four atoms per unit cell, i.e., it has a four-
 * atom basis.  One atom is at the center, and four atoms are in each of the faces
 * formed by a*c.  There is not an atom in either of the a*a faces.  The fourth atom
 * in the basis is formed by the atoms at the corners of the unit cell. A good
 * visualization of the structure is available at http://cst-www.nrl.navy.mil/lattice/
 * struk/a5.html, in the 1998 communication by T.F. Fassler and C. Kronseder, and in
 * the 2003 review by A. Mujica, A. Rubio, A. Munoz, and R.J. Needs.  
 * 
 * It should be noted that the beta-tin structure is not specific to the element tin.
 * Silicon, germanium and tin have several solid phases, and each experience the 
 * diamond-cubic crystal structure and the beta-tin crystal structure. The diamond 
 * crystal structure for tin is often refered to as the alpha-tin phase.
 * 
 * This class uses the primitive vectors as given in PrimitiveTetragonal.java.
 * 
 * This class was adapted from D.A. Kofke's BasisCubic.java by K.R. Schadel and A. 
 * Schultz July 2005.
 */
 
public class BasisBetaSnA5 implements Basis {
    
    /**
     * @param primitive Primitive of the tetragonal lattice housing this basis.
     * Needed to ensure that separation of basis atoms is consistent with
     * spacing of atoms on lattice.
     * @param scaledCoordinates basis coordinates for the case in which the
     * primitive lattice constant (getSize) is unity.  Given instance is kept
     * for interal representation of basis, so changes to scaledCoordinates
     * will affect the basis.
     */
    public BasisBetaSnA5(PrimitiveTetragonal primitive) {
        this.primitive = primitive;
        // Is size needed? It seems as though coordinates is defined twice.
        size = unscaledCoordinates.length;
        coordinates = new Vector[size];
        for(int i=0; i<size; i++) {
            this.coordinates[i] = (Vector)unscaledCoordinates[i].clone();
        }
        oldAB = 1.0;
        oldC = 1.0;
    }
    
    public int size() {
        return size;
    }
    
    /**
     * The following vectors are the *unscaled* position vectors of the four atom sites 
     * chosen to constitute the basis.  The x and y coordinates of these vectors are to
     * be multiplied by the lattice parameter a, and the z coordinates are to be
     * multiplied by the lattice parameter c.  These vectors, with their coordinates 
     * not multiplied by the appropriated lattice constant, by themselves, would not 
     * define the positions of the basis atoms, no matter what the size of the unit 
     * cell.  Each of these vectors must be multiplied by the vector <a,a,c>, as done 
     * in the next field, in order to make any sense. 
     * 
     * The first vector corresponds to the bottom corner of the unit cell, situated at
     * the origin.  The second and third are atoms situated in an a*c face of the unit
     * cell; note the differing values of z.  The final vector is for the position of 
     * the atom at the center of the unit cell.  
     */
    
    private static final Vector3D[] unscaledCoordinates = new Vector3D[] {
			new Vector3D(0.00, 0.00, 0.00),
			new Vector3D(0.5, 0.00, 0.25),
			new Vector3D(0.00, 0.5, 0.75),
			new Vector3D(0.5, 0.5, 0.5)
    };
    
    public Vector[] positions() {
        double AB = primitive.getA();
        double C = primitive.getC();
        if(AB != oldAB || C != oldC) {
        	Vector ABC = primitive.space.makeVector();
        	ABC.setX(0, AB);
        	ABC.setX(1, AB);
        	ABC.setX(2, C);
            for(int i=0; i<size; i++) {
                coordinates[i].E(unscaledCoordinates[i]);
                coordinates[i].TE(ABC);
            }
            oldAB = AB;
            oldC = C;
        }
        return coordinates;
    }
    
    private int size;
    private Vector[] coordinates;
    private PrimitiveTetragonal primitive;
    private double oldAB;
    private double oldC;
}