package etomica;

/**
 * Holder of atoms not being used in a phase.  Normally
 * associated with an AtomFactory.
 *
 * @author David Kofke
 */
public class AtomReservoir extends Atom {
    
    private int capacity;
    private final AtomList atomList = new AtomList();
    public final AtomFactory factory;
//    public final AtomTreeNodeGroup node;//shadow superclass field of same name to avoid casts
    
    /**
     * Constructs with default maximum capacity of 80.
     */
    public AtomReservoir(AtomFactory factory) {
        this(factory, 80);
    }
    /**
     * Construct reservoir that will accept up to the given number of atoms.
     */
    public AtomReservoir(AtomFactory factory, int capacity) {
        super(factory.parentSimulation().space, AtomType.NULL, new NodeFactory(), null);
        this.factory = factory;
        if(capacity < 0) capacity = 0;
        this.capacity = capacity;
 //       node = (AtomTreeNodeGroup)super.node;
    }
    
    /**
     * Sets parent of given atom to the reservoir, or disposes
     * of atom if reservoir is full.
     */
    public void addAtom(Atom atom) {
        if(atom == null) return;
        if(((AtomTreeNodeGroup)node).childList.size() >= capacity) atom.node.dispose();
        else atom.node.setParent(((AtomTreeNodeGroup)node));
        //restore atom to condition when built
//        if(atom instanceof AtomGroup) ((AtomGroup)atom).creator().renew(atom);
        //add to reservoir
    }
    
    /**
     * Returns the most recently added atom.
     */
    public Atom getAtom() {return ((AtomTreeNodeGroup)node).childList.getLast();}
    
    /**
     * Indicates whether reservoir contains any atoms.
     */
    public boolean isEmpty() {return ((AtomTreeNodeGroup)node).childList.size() == 0;}
    
    /**
     * Sets the maximum number of atoms that the reservoir will hold.
     * If current number is greater, it removes atoms (least recently added
     * removed first) until maximum is reached.
     */
    public void setCapacity(int i) {
        capacity = i; 
        if(capacity < 0) capacity = 0;
        while(atomList.size() > capacity) atomList.removeFirst();
    }
    /**
     * Returns the maximum number of atoms the reservoir will hold.
     */
    public int getCapacity() {return capacity;}
    
    /**
     * Returns "Reservoir"
     */
     public String signature() {return "Reservoir";}

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
            if(childList.size() > reservoir.getCapacity()) {
                atom.node.dispose();
            }
        }
        
        public void removeAtomNotify(Atom atom) {
        }
        
        /**
         * Returns null.
         */
         public final Phase parentPhase() {return null;}

        /**
        * Overrides super class method and terminates recursive call to identify
        * a constituent atom's species.
        */
        public final Species parentSpecies() {return reservoir.factory.species();}
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