package etomica;

/**
 * Class responsible for building new instances of the atoms (or atom groups)
 * that are collected in a given AtomGroup.
 */
public abstract class AtomFactory {
    
    private AtomReservoir reservoir = new AtomReservoir();
    
    public Atom makeAtom(AtomGroup parent, int index) {
        Atom atom = reservoir.removeAtom();
        if(atom == null) {
            atom = makeNewAtom(parent, index);
        }
        else {
            atom.setParentGroup(parent);
            atom.setIndex(index);
        }
        return atom;
    }
    
    public abstract Atom makeNewAtom(AtomGroup parent, int index);
    
    /**
     * Indicator of whether this factory makes atom groups, or if it 
     * makes atoms that are not groups of other atoms (instance of AtomGroup).
     */
    public abstract boolean producesAtomGroups();
    
    
    public abstract boolean vetoAddition(Atom a); //be sure to check that a is non-null
    
    public AtomReservoir reservoir() {return reservoir;}
    
}