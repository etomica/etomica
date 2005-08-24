package etomica.atom;
import etomica.species.Species;
import etomica.units.Dimension;

/**
 * The SpeciesAgent is a representative of the species in each phase.
 * The agent handles addition, deletion, link-list ordering, counting, etc. of 
 * molecules in a phase.  Each phase has an agent from every species instance.
 * 
 * @author David Kofke
 */
 
public final class SpeciesAgent extends Atom {

    public SpeciesAgent(AtomType type, Species species) {
        super(null, type, NODE_FACTORY,  AtomLinker.FACTORY);
        type.setSpecies(species);
    }
        
    public final AtomFactory moleculeFactory() {return type.getSpecies().moleculeFactory();}
      
    public SpeciesAgent nextSpecies() {return (SpeciesAgent)seq.next.atom;}
    public int moleculeCount() {return ((AtomTreeNodeGroup)node).childAtomCount();}
    public Atom firstMolecule() {return ((AtomTreeNodeGroup)node).childList.getFirst();}
    public Atom lastMolecule() {return ((AtomTreeNodeGroup)node).childList.getLast();}
    public Atom randomMolecule() {return ((AtomTreeNodeGroup)node).childList.getRandom();}
            
    public Atom addNewAtom() {
        Atom aNew = moleculeFactory().makeAtom();
        aNew.node.setParent((AtomTreeNodeGroup)node);
        return aNew;
    }    
    
    /**
    * Sets the number of molecules for this species.  Makes the given number
    * of new molecules, linked-list orders and initializes them.
    * Any previously existing molecules for this species in this phase are abandoned
    * Any links to molecules of next or previous species are maintained.
    * Takes no action at all if the new number of molecules equals the existing number
    *
    * @param n  the new number of molecules for this species
    */
    public void setNMolecules(int n) {
        lastNMolecules = n;
        AtomTreeNodeGroup treeNode = (AtomTreeNodeGroup)node;
        while(treeNode.childList.size() > 0) treeNode.childList.getLast().node.dispose();
        if(n > treeNode.childAtomCount()) {
            for(int i=treeNode.childAtomCount(); i<n; i++) addNewAtom();
        }
        else if(n < treeNode.childAtomCount()) {
            if(n < 0) n = 0;
            for(int i=treeNode.childAtomCount(); i>n; i--) treeNode.childList.getLast().node.dispose();
        }
    }
    
    public void makeMolecules() {
        int nMolecules = (lastNMolecules>-1) ? lastNMolecules : type.getSpecies().getNMolecules(); 
        setNMolecules(nMolecules);
    }
    
    public int getNMolecules() {return moleculeCount();}
    public Dimension getNMoleculesDimension() {return Dimension.QUANTITY;}

    /**
     * Special AtomTreeNode class for SpeciesAgent.
     */
    public static final class AgentAtomTreeNode extends AtomTreeNodeGroup {
        
        private AgentAtomTreeNode(Atom atom) {
            super(atom);
        }

       /**
        * Overrides parent class method and terminates recursive call to identify this
        * as a constituent atom's species agent.
        */
        public final SpeciesAgent parentSpeciesAgent() {return (SpeciesAgent)this.atom;}
        
        /**
         * Throws a RuntimeException, because a species agent is not contained within a molecule.
         */
        public final Atom parentMolecule() {
            throw new RuntimeException("Error:  Unexpected call to parentMolecule in SpeciesAgent");
        }
    }
    
    private final static AtomTreeNodeFactory NODE_FACTORY = new AtomTreeNodeFactory() {
        public AtomTreeNode makeNode(Atom atom) {
            return new AgentAtomTreeNode(atom);
        }
    };
    public AtomLinker.Tab firstLeafAtomTab;
    private int lastNMolecules = -1;
    
} //end of SpeciesAgent
