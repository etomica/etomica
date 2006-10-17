package etomica.lattice.crystal;
import etomica.lattice.Basis;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

/**
 * Basis that helps construct a beta-tin crystal.  
 * 
 * The beta-tin phase of tin has the beta-tin crystal structure, designated 
 * by the Pearson symbol tI4 (ASM Handbook).  The "t" indicates that the structure 
 * belongs to the tetragonal crystal system.  The "I" indicates that the space 
 * lattice, or lattice type, is body-centered (There is an atom at the center 
 * of the unit cell).  The "4" indicates that there are four atoms per unit 
 * cell.  An eighth of an atom is at each of the eight corners, one atom is 
 * at the center, and half of an atom is at each of the four a*c faces.  
 * Note:  These atoms are NOT centered in the a*c faces.  There is not an 
 * atom in either of the a*a faces. The Strukturbericht symbol of the beta-Sn 
 * crystal structure is A5.  A good visualization of the structure is available 
 * at http://cst-www.nrl.navy.mil .
 * 
 * Each atom is at the center of a distorted octahedron.  Consider the atom 
 * at the center of a unit cell.  The two tips of its octahedron (the top and
 * bottom of the two joined pyramids, if you will) are formed by
 * the atoms at the centers of the unit cells directly above and below the atom 
 * of interest's unit cell.  The four corners of what would be the pyriamids'
 * bases are formed by the atoms in the a*c faces of the center atom's unit cell.
 * The center atom is slightly closer to these four atoms than it is to the two 
 * atoms forming the tips of the octahedron.  The first-nearest-neighbor 
 * distance for beta-Sn is 3.02 Angstroms (ASM Handbook), while the second-
 * nearest-neighbor distance is c (3.1815 Angstroms).  After doing a little 
 * math, one will find that two of the four atoms in the a*c faces are each 
 * 1/4*c below the center of their faces, and the other two are 1/4*c above the 
 * center of their faces.
 * 
 * It should be noted that the beta-tin structure is not specific to the element tin.
 * Silicon, germanium and tin have several solid phases, and each experience the 
 * diamond-cubic crystal structure and the beta-tin crystal structure (Mujica,
 * Rubio, Munoz, and Needs, 2003). The diamond crystal structure for tin is 
 * often refered to as the alpha-tin phase (gray tin).
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
     * multiplied by the lattice parameter c.  Each of these vectors are multiplied 
     * by the vector <a,a,c> in the next field.
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
        double AB = primitive.getAB();
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