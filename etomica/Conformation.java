package etomica;

import etomica.atom.AtomList;

/**
 * General class for assignment of coordinates to a group of atoms.  Puts
 * a list of atoms in a standard conformation, which then can be manipulated
 * further by a Configuration class to place the molecules in a phase, or by
 * a super-conformation class that arranges these atoms with other ones in a 
 * molecule.<br>  
 * This class is used by an AtomFactory to arrange into a standard configuration
 * each atom group that it builds.  
 */
 
public abstract class Conformation implements java.io.Serializable {

    public Conformation(Space space) {
        this.space = space;
    }

    /**
     * Defined by subclass to assign coordinates to the atoms in the given list.
     */
    public abstract void initializePositions(AtomList atomList);

    protected final Space space;
    
    /**
     * Conformation that does nothing to the atom positions.
     */
    public final static Conformation NULL = new Conformation(null) {
        public void initializePositions(AtomList list) {}
    };
}//end of Configuration
