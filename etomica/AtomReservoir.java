package etomica;

/**
 * Holder of atoms not being used in a phase.  Normally
 * associated with an AtomFactory.
 *
 * @author David Kofke
 */
public class AtomReservoir {
    
    private int maximumCapacity;
    private final AtomList atomList = new AtomList();
    
    /**
     * Constructs with default maximum capacity of 80.
     */
    public AtomReservoir() {
        this(80);
    }
    /**
     * Construct reservoir that will accept up to the given number of atoms.
     */
    public AtomReservoir(int maximumCapacity) {
        if(maximumCapacity < 0) maximumCapacity = 0;
        this.maximumCapacity = maximumCapacity;
    }
    
    /**
     * Sets parent of given atom to null and adds it to reservoir.
     */
    public void addAtom(Atom atom) {
        if(atom == null) return;
        atom.node.setParent((AtomTreeNodeGroup)null);
        if(atomList.size() >= maximumCapacity) return;
        //restore atom to condition when built
//        if(atom instanceof AtomGroup) ((AtomGroup)atom).creator().renew(atom);
        //add to reservoir
        atomList.add(atom.seq);
    }
    
    /**
     * Returns the most recently added atom.
     */
    public Atom removeAtom() {
        Atom atom = null;
        while(atom == null && atomList.size() > 0) {
            atom = atomList.removeLast();
            //check that atom wasn't placed in another group without taking it from reservoir
            if(atom.node.parentGroup() == null) return atom;
            else atom = null;
        }
        return null;//if reach here, list is empty
    }
    
    /**
     * Indicates whether reservoir contains any atoms.
     */
    public boolean isEmpty() {return atomList.size() == 0;}
    
    /**
     * Sets the maximum number of atoms that the reservoir will hold.
     * If current number is greater, it removes atoms (least recently added
     * removed first) until maximum is reached.
     */
    public void setMaximumCapacity(int i) {
        maximumCapacity = i; 
        if(maximumCapacity < 0) maximumCapacity = 0;
        while(atomList.size() > maximumCapacity) atomList.removeFirst();
    }
    /**
     * Returns the maximum number of atoms the reservoir will hold.
     */
    public int getMaximumCapacity() {return maximumCapacity;}
    
}//end of AtomReservoir