package etomica;
import etomica.atom.AtomFactory;
import etomica.atom.AtomLinker;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomTreeNode;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTreeNodeFactory;
import etomica.units.Dimension;

/**
 * The SpeciesAgent is a representative of the species in each phase.
 * The agent handles addition, deletion, link-list ordering, counting, etc. of 
 * molecules in a phase.  Each phase has an agent from every species instance.
 * 
 * @author David Kofke
 */
 
public final class SpeciesAgent extends Atom {

    protected Integrator integrator;
    public AtomLinker.Tab firstLeafAtomTab;
    private final Phase phase;
    
    public SpeciesAgent(Space space, Species species, Phase phase, int nMolecules) {
        super(space, new AtomTypeGroup(null), new NodeFactory(species), 
                AtomSequencerFactory.SIMPLE);
        this.phase = phase;
        ((AtomTypeGroup)type).childrenAreGroups = species.moleculeFactory().isGroupFactory();
    }
        
    public final AtomFactory moleculeFactory() {return node.parentSpecies().moleculeFactory();}
      
    public SpeciesAgent nextSpecies() {return (SpeciesAgent)seq.next.atom;}
    public int moleculeCount() {return ((AtomTreeNodeGroup)node).childAtomCount();}
    public Atom firstMolecule() {return ((AtomTreeNodeGroup)node).firstChildAtom();}
    public Atom lastMolecule() {return ((AtomTreeNodeGroup)node).lastChildAtom();}
    public Atom randomMolecule() {return ((AtomTreeNodeGroup)node).randomAtom();}
    public Atom getMolecule(int i) {return ((AtomTreeNodeGroup)node).getAtom(i);}
            
    public Atom addNewAtom() {
        Atom aNew = moleculeFactory().makeAtom();
        phase.addMolecule(aNew, this);
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
    * @see #makeMolecule
    * @see #deleteMolecule
    * @see #addMolecule
    */
    public void setNMolecules(int n) {
        AtomTreeNodeGroup treeNode = (AtomTreeNodeGroup)node;
        while(treeNode.childList.size() > 0) treeNode.lastChildAtom().dispose();
        if(n > treeNode.childAtomCount()) {
            for(int i=treeNode.childAtomCount(); i<n; i++) addNewAtom();
        }
        else if(n < treeNode.childAtomCount()) {
            if(n < 0) n = 0;
            for(int i=treeNode.childAtomCount(); i>n; i--) treeNode.lastChildAtom().dispose();
        }
        
        //reconsider this
        //node.parentPhase().configuration.initializeCoordinates(this);
        node.parentPhase().reset();
    }
    
    public int getNMolecules() {return moleculeCount();}
    public Dimension getNMoleculesDimension() {return Dimension.QUANTITY;}

    /**
     * Special AtomTreeNode class for SpeciesAgent.
     */
    public static final class AgentAtomTreeNode extends AtomTreeNodeGroup {
        
        private final Species parentSpecies;//remove this 2/05
        private AgentAtomTreeNode(Species parentSpecies, Atom atom) {
            super(atom);
            this.parentSpecies = parentSpecies;
            depth = 1;
            setIndex(parentSpecies.getIndex());
        }

        /**
        * Overrides super class method and terminates recursive call to identify
        * a constituent atom's species.
        */
        public final Species parentSpecies() {return parentSpecies;}
        /**
        * Overrides parent class method and terminates recursive call to identify this
        * as a constituent atom's species agent.
        */
        public final SpeciesAgent parentSpeciesAgent() {return (SpeciesAgent)this.atom;}
        
        /**
        * Returns null, since a species agent is not contained within a molecule.
        */
        public final Atom parentMolecule() {return null;}
    }
    
    private static final class NodeFactory implements AtomTreeNodeFactory {
        Species species;
        NodeFactory(Species s) {
            species = s;
        }
        public AtomTreeNode makeNode(Atom atom) {
            return new AgentAtomTreeNode(species, atom);
        }
    }

} //end of SpeciesAgent
