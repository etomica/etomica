package etomica.config;

import etomica.api.IAtomSet;
import etomica.space.Space;

/**
 * General class for assignment of coordinates to a group of atoms.  Puts
 * a list of atoms in a standard conformation, which then can be manipulated
 * further by a Configuration class to place the molecules in a box, or by
 * a super-conformation class that arranges these atoms with other ones in a 
 * molecule.
 * <p>
 * Specific arrange of atoms performed by the class is defined by the subclass
 * via its implementation of the <tt>initializePositions</tt> method.
 * <p>  
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
    public abstract void initializePositions(IAtomSet atomList);

    protected final Space space;
    
}
