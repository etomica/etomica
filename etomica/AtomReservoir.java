package etomica;

/**
 * Holder of atoms not being used in a phase.  Normally
 * associated with an AtomFactory.
 *
 * @author David Kofke
 */
public class AtomReservoir extends AtomGroup {
    
    private Atom[] atoms;
    private int count = 0;
    private int maximumCapacity;
    
    public AtomReservoir() {
        this(10,80);
    }
    public AtomReservoir(int initialCapacity, int maximumCapacity) {
        super(null, 0, AtomType.NULL);
        this.maximumCapacity = maximumCapacity;
        if(initialCapacity > maximumCapacity) initialCapacity = maximumCapacity;
        atoms = new Atom[initialCapacity];
    }
    
    public void addAtom(Atom atom) {
        if(count >= maximumCapacity || atom == null) return;
        atoms[count++] = atom;
        if(count == atoms.length) {
            int newCapacity = 2*atoms.length;
            if(newCapacity > maximumCapacity) newCapacity = maximumCapacity;
            Atom[] newAtoms = new Atom[newCapacity];
            for(int i=0; i<atoms.length; i++) {
                newAtoms[i] = atoms[i];
            }
            atoms = newAtoms;
        }
    }
    
    public Atom removeAtom() {
        Atom atom = null;
        while(atom == null && count > 0) {
            atom = atoms[count--];
            if(atom.parentGroup() != this) atom = null; //check that atom wasn't placed in another group without taking it from reservoir
        }
        return atom;
    }
    
    public boolean isEmpty() {return count == 0;}
    
    public Phase parentPhase() {return null;}
    public Species parentSpecies() {return null;}
    
    public void setMaximumCapacity(int i) {maximumCapacity = i;}
    public int getMaximumCapacity() {return maximumCapacity;}
    
}//end of AtomReservoir