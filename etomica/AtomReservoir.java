package etomica;

/**
 * Holder of atoms not being used in a phase.  Normally
 * associated with an AtomFactory.
 *
 * @author David Kofke
 */
public class AtomReservoir extends AtomGroup {
    
    private Atom[] atoms;
    private int count = 0;  //index of where next atom is to be placed in atoms array
    private int maximumCapacity;
    private Simulation parentSimulation;
    private Space space;
    
    public AtomReservoir(Space space) {
        this(space, 10,80);
        
    }
    public AtomReservoir(Space s, int initialCapacity, int maximumCapacity) {
        super(s, AtomType.NULL);
        space = s;
        if(initialCapacity < 1) initialCapacity = 1;
        if(maximumCapacity < 1) maximumCapacity = 1;
        this.maximumCapacity = maximumCapacity;
        if(initialCapacity > maximumCapacity) initialCapacity = maximumCapacity;
        atoms = new Atom[initialCapacity];
    }
    
    public void addAtom(Atom atom) {
        //safety check
        if(count >= maximumCapacity || atom == null) return;
        //restore atom to condition when built
        if(atom instanceof AtomGroup) ((AtomGroup)atom).creator().renew(atom);
        //add to reservoir
        atoms[count++] = atom;
        //check reservoir size and expand to accommodate next addition if now full
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
    
    /**
     * Override superclass method to terminate chain of notification to parent groups.
     * Otherwise performs no action.
     */
    public void addAtomNotify(Atom atom) {
        //perhaps throw an exception, since this shouldn't be called
    }
    
    /**
     * Override superclass method to terminate chain of notification to parent groups.
     * Otherwise performs no action.
     */
    public void removeAtomNotify(Atom atom) {
        //perhaps throw an exception, since this shouldn't be called
    }
    
    
    public Atom removeAtom() {
        Atom atom = null;
        while(atom == null && count > 0) {
            atom = atoms[--count];
            //check that atom wasn't placed in another group without taking it from reservoir
            if(atom.parentGroup() != this) atom = null; 
        }
        return atom;
    }
    
    public boolean isEmpty() {return count == 0;}
    
    public Phase parentPhase() {return null;}
    public Species parentSpecies() {return null;}
    public SpeciesAgent parentSpeciesAgent() {return null;}
    public Simulation parentSimulation() {return null;}
    
    public void setMaximumCapacity(int i) {
        maximumCapacity = i; 
        if(maximumCapacity < 1) maximumCapacity = 1;
    }
    public int getMaximumCapacity() {return maximumCapacity;}
    
}//end of AtomReservoir