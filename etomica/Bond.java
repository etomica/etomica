package etomica;

/**
 * Class used to hold information about bonds between atoms.  Allow atoms
 * to form permanent links to other atoms, and provides a means to hold
 * information describing the link.  Bond classes may be used to describe
 * physical covalent bonds, or weaker associations.
 * The basic Bond class describes a two-atom bond.  Other types of bond
 * classes (as yet undeveloped) would be needed to apply to multi-atom 
 * bonds (e.g. pi bonds).
 *
 * @author David Kofke
 */
 
public class Bond implements java.io.Serializable {

    public BondLinker link1, link2;
    public static final Bond instance = new Bond();
    
    //should have a bond reservoir
    
    private Bond() {}
    
    private Bond(Atom a1, Atom a2) {
        if(a1.preceeds(a2)) {
            link1 = BondLinker.getNew(this, a1, a1.firstUpBond, null);
            link2 = BondLinker.getNew(this, a2, a2.firstDownBond, null);
            a1.firstUpBond = link1;
            a2.firstDownBond = link2;
        }
        else {
            link1 = BondLinker.getNew(this, a1, a1.firstDownBond, null);
            link2 = BondLinker.getNew(this, a2, a2.firstUpBond, null);
            a1.firstDownBond = link1;
            a2.firstUpBond = link2;
        }
    }
    
    public final Atom atom1() {return link1.atom;}
    public final Atom atom2() {return link2.atom;}
    
    public static Bond makeBond(Atom a1, Atom a2) {
    //should have a bond reservoir
        return new Bond(a1, a2);
    }
    
    /**
     * Returns the atom partnered to the given atom via this bond.
     * Returns null if the given atom is not bonded via this bond.
     */
    public Atom partner(Atom atom) {
        if(atom == link1.atom) return link2.atom;
        else if(atom == link2.atom) return link1.atom;
        else return null;
    }
    
    public void breakBond() {
        link1.delete();
        link2.delete();
    }

//consider for multi-atom bonds
/*    protected Atom[] atoms;
    
    public Atom[] atoms() {return atoms;}
    public void setAtoms(Atom[] atoms) {
        this.atoms = atoms;
    }
*/    
}