package etomica;

/**
 * Holder of atoms not being used in a phase.  Normally
 * associated with an AtomFactory.
 *
 * @author David Kofke
 */
public class AtomReservoir extends Atom {
    
    private int maximumCapacity;
    private final AtomList atomList = new AtomList();
    public final AtomFactory factory;
    public final AtomTreeNodeGroup node;//shadow superclass field of same name to avoid casts
    
    /**
     * Constructs with default maximum capacity of 80.
     */
    public AtomReservoir(AtomFactory factory) {
        this(factory, 80);
    }
    /**
     * Construct reservoir that will accept up to the given number of atoms.
     */
    public AtomReservoir(AtomFactory factory, int maximumCapacity) {
        super(factory.parentSimulation().space, AtomType.NULL, new NodeFactory(), null);
        this.factory = factory;
        if(maximumCapacity < 0) maximumCapacity = 0;
        this.maximumCapacity = maximumCapacity;
        node = (AtomTreeNodeGroup)super.node;
    }
    
    /**
     * Sets parent of given atom to null and adds it to reservoir.
     */
    public void addAtom(Atom atom) {
        if(atom == null) return;
        if(node.childList.size() >= maximumCapacity) atom.node.destroy();
        else atom.node.setParent(node);
        //restore atom to condition when built
//        if(atom instanceof AtomGroup) ((AtomGroup)atom).creator().renew(atom);
        //add to reservoir
    }
    
    /**
     * Returns the most recently added atom.
     */
    public Atom getAtom() {return node.childList.getLast();}
    
    /**
     * Indicates whether reservoir contains any atoms.
     */
    public boolean isEmpty() {return node.childList.size() == 0;}
    
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

    /**
     * Special AtomTreeNode class for AtomReservoir.
     */
    private static final class ReservoirAtomTreeNode extends AtomTreeNodeGroup {
        
        private final AtomReservoir reservoir;
        
        ReservoirAtomTreeNode(Atom atom) {
            super(atom, null);
            depth = 0;
            reservoir = (AtomReservoir)atom;
        }
        
        public void addAtomNotify(Atom atom) {
            if(childList.size() > reservoir.getMaximumCapacity()) {
                atom.node.destroy();
            }
        }
        
        public void removeAtomNotify(Atom atom) {
        }

        /**
        * Overrides super class method and terminates recursive call to identify
        * a constituent atom's species.
        */
        public final Species parentSpecies() {return reservoir.factory.parentSpecies();}
        /**
        * Overrides super class method and terminates recursive call to identify
        * a constituent atom's simulation.
        */
        public Simulation parentSimulation() {return reservoir.factory.parentSimulation();}
        
        /**
         * Returns null.
         */
        public final SpeciesAgent parentSpeciesAgent() {return null;}
        
        /**
         * Returns null, since a species agent is not contained within a molecule.
         */
        public final Atom parentMolecule() {return null;}
                
    }//end of ReservoirAtomTreeNode
    
    private static final class NodeFactory implements AtomTreeNode.Factory {
        public AtomTreeNode makeNode(Atom atom, AtomTreeNodeGroup dummy) {
            return new ReservoirAtomTreeNode(atom);
        }
    }//end of NodeFactory
    
}//end of AtomReservoir