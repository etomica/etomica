package etomica;
/**
 * An association of three atoms, with methods related to their geometry.
 * 
 * @author David Kofke
 * @author Jhumpa Adhikari
 */
public final class Atom3 implements AtomSet, java.io.Serializable {
    
    public static String getVersion() {return "Atom3:01.08.03";}
    
    public Atom atom1, atom2, atom3;
    public final Space.CoordinatePair c12, c13, c23;

   /**
    * Constructs an Atom3 for the given phase, but with no designated atoms.
    */
    public Atom3(Space space) {
        c12 = space.makeCoordinatePair();
        c13 = space.makeCoordinatePair();
        c23 = space.makeCoordinatePair();
    }
    
   /**
    * Redefines the Atom3 to correspond to the given atoms
    */
    public void reset(Atom a1, Atom a2, Atom a3) {
        atom1 = a1; 
        atom2 = a2;
        atom3 = a3;
        reset();
    }
    /**
     * Resets the coordinate pair for the current values of the atoms
     */
    public void reset() {
        c12.reset(atom1.coord, atom2.coord);
        c13.reset(atom1.coord, atom3.coord);
        c23.reset(atom2.coord, atom3.coord);
    }
    
    /**
     * Returns true if one of the three atoms is the given atom.
     */
	public boolean contains(Atom a) {
		return atom1 == a || atom2 == a || atom3 == a;
	}
	
	/**
	 * Returns 3, indicating this is an atom set of three atoms. Part of AtomSet
	 * interface.
	 */
	public final int nBody() {return 3;}
	
	

    
}  //end of  Atom3