package etomica;
/**
 * An association of two atoms.  Each AtomPair holds one CoordinatePair (obtained from a Space class),
 * which has all the methods needed to compute atom distance and momentum-difference vectors, dot products, etc.
 */

/* History
 * 07/17/03 (DAK) modified reset() methods to return this AtomPair
 * 08/21/03 (DAK) modified reset() method to invoke cPair.reset() rather than
 * cPair.reset(atom1.coord, atom2.coord); modified reset(Atom, Atom) accordingly
 * 08/29/03 (DAK) added reset2(Atom) method; made atom2 private
 */
 
public final class AtomPair implements AtomSet, java.io.Serializable {
    public static String getVersion() {return "AtomPair:01.06.25";}
    public Atom atom1;
    private Atom atom2;
    public final Space.CoordinatePair cPair;

   /**
    * Constructs an AtomPair for the given phase, but with no designated atoms.
    */
    public AtomPair(Space space) {
        cPair = space.makeCoordinatePair();
    }
    /**
     * Constructs an AtomPair using the given atoms.  Assumes that the atoms are in the same phase.
     * The method atom1() will return the first atom in the argument list here, and atom2() the second.
     */
    public AtomPair(Atom a1, Atom a2) {  //Assumes a1 and a2 are in same phase
        cPair = a1.node.parentSimulation().space().makeCoordinatePair();
        reset(a1, a2);
    }
    /**
     * Constructs an AtomPair using the given atoms and coordinate pair.  The coordinate pair
     * is assumed to correspond to the given atoms.  Passing it here can save on the overhead
     * of making it if it is already in place
     * The method atom1() will return the first atom in the argument list here, and atom2() the second.
     */
    public AtomPair(Atom a1, Atom a2, Space.CoordinatePair c) {atom1 = a1; atom2 = a2; cPair = c;}
    
    /**
     * Clones this atomPair without cloning the atoms or their coordinates.
     * The returned atomPair refers to the same pair of atoms as the original.
     * This can be used to make a working copy of an atomPair that is returned by an atomPair iterator.
     * Method is called "copy" instead of "clone" because whole object tree isn't cloned.
     */
    public AtomPair copy() {
        return new AtomPair(atom1, atom2, cPair.copy());  //cannot use super.clone() because cPair (being final) cannot then be changed to a clone of cPair
    }

	/**
	 * Returns true if the given atom is atom1 or atom2.
	 * @see etomica.AtomSet#contains(Atom)
	 */
	public boolean contains(Atom a) {
	//	return atom2 == a;
		return atom1 == a || atom2 == a;
	}
	/**
	 * Returns 2, indicating that this is an atom set containing two atoms.
	 * Part of the AtomSet interface.
	 * 
	 */
	public final int nBody() {return 2;}
	
   /**
    * Redefines the atom pair to correspond to the given atoms
    */
    public AtomPair reset(Atom a1, Atom a2) {
        atom1 = a1; 
        atom2 = a2;
        if(a2 != null) cPair.reset(a1.coord, a2.coord);
//        if(a2 != null) reset();
        return this;
    }
    /**
     * Defines the second atom as the given one, and resets the pair using it
     * and the presently defined first atom.
     * @param a2
     * @return AtomPair
     */
    public AtomPair reset2(Atom a2) {
    	atom2 = a2;
    	cPair.reset(atom1.coord, a2.coord);
    	return this;
    }
    /**
     * Resets the coordinate pair for the current values of the atoms
     */
    public AtomPair reset() {
//        cPair.reset(atom1.coord, atom2.coord); //08/21/03 (DAK) changed to reset()
		cPair.reset();
        return this;
    }
    
    public void setAtom1(Atom atom) {atom1 = atom;}
    
    public void setAtom2(Atom atom) {atom2 = atom;}
    
    public void setAtoms(Atom atom1, Atom atom2) {
    	this.atom1 = atom1;
    	this.atom2 = atom2;
    }
    /**
     * @return the square of the distance between the atoms, |r1 - r2|^2
     */
    public final double r2() {return cPair.r2();}
    
    /**
     * @return the square of the velocity-difference between the atoms, |p1/m1 - p2/m2|^2
     */
    public final double v2() {return cPair.v2();}
    
    /**
     * @return the dot product of the distance and the velocity difference, (r1 - r2).(p1/m1 - p2/m2)
     */
    public final double vDotr() {return cPair.vDotr();}
    
    /**
     * @return the vector distance between the atoms, r2 - r1
     */
    public final Space.Vector dr() {return cPair.dr();}
    
    /**
     * @return the i<sup>th</sup> component of the distance vector r2 - r1, where i=0 is the first component
     */
    public final double dr(int i) {return cPair.dr(i);}
    
    /**
     * @return the i<sup>th</sup> component of the velocity-difference vector p2/m2 - p2/m1, where i=0 is the first component
     */
    public final double dv(int i) {return cPair.dv(i);}
    
    /**
     * @return the first atom of the atom pair (as set when invoking the constructor or the reset method)
     */
    public final Atom atom1() {return atom1;}
    
    /**
     * @return the second atom of the atom pair (as set when invoking the constructor or the reset method)
     */
    public final Atom atom2() {return atom2;}
    
    /**
     * @return a String formed from the toString methods of atom1 and atom2.
     */
    public String toString() {return atom1.toString()+" "+atom2.toString();}
    /**
     * Sorts by separation distance all the atom pairs produced by an atomPair iterator.
     * Returns the first element of a linked list of atomPair(Linker)s, sorted by increasing distance
     * Perhaps better to do this using java.util.Collections (in java 1.2 API)
     */
    public static AtomPairLinker distanceSort(AtomPairIterator api) {
        if(!api.hasNext()) return null;
        AtomPairLinker firstLink = new AtomPairLinker(api.next().copy());
        while(api.hasNext()) {                      //loop through all pairs generated by api
            AtomPair nextPair = api.next().copy();  //make a copy of pair for use in ordered list
            //Insert pair into ordered list in proper location
            AtomPairLinker previous = null;  //holds value from previous iteration of this for-loop
            boolean inserted = false;
            for(AtomPairLinker link=firstLink; link!=null; link=link.next()) {
                if(nextPair.r2() < link.pair().r2()) {  //insert nextPair before current pair
                    if(previous == null) {firstLink = new AtomPairLinker(nextPair,firstLink);} //nextPair is new firstPair, to be followed by old firstPair
                    else {previous.setNext(new AtomPairLinker(nextPair,link));}  //place nextPair between last and pair
                    inserted = true;
                    break;  //break out of for-loop
                }
                previous = link;
            }  //end of for loop
            if(!inserted) //reached end of list without inserting;
                previous.setNext(new AtomPairLinker(nextPair));   //insert after last link
        }
        return firstLink;
    }//end of distanceSort
    
}  //end of  AtomPair