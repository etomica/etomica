package etomica.atom;

import etomica.Atom;
import etomica.Configuration;
import etomica.Space;
import etomica.Species;
import etomica.Atom.AgentSource;
import etomica.atom.AtomTreeNode.Factory;

/**
 * Class responsible for building new instances of the atoms (or atom groups)
 * that are collected in a given AtomGroup.
 *
 * @author David Kofke
 */
 
 /* History
  * 10/18/02 (DAK) Modified to remove Simulation instance, using Space instance instead.
  *          Many other classes (other factories, lattice classes, modified to be consistent 
  *          with this change.
  * 08/12/03 (DAK) Re-introduced simulation field, allowing it to be null.  Used
  * by AtomType to determine simulation instance, which in turn is referenced by
  * AtomTreeNode.
  * 08/26/03 (DAK) Added constructors that take nodeFactory argument.  
  */
public abstract class AtomFactory {
    
    protected final AtomReservoir reservoir;
    public final Space space;
    private Species species;
    protected Configuration configuration;
    protected AtomSequencerFactory sequencerFactory;
    protected AtomTreeNodeGroup.Factory nodeFactory;
    private Atom.AgentSource[] agentSource = new Atom.AgentSource[0];
    protected final AtomTypeGroup groupType = new AtomTypeGroup(this);
    protected final AtomTypeSphere spheretype = new AtomTypeSphere(this);
    protected AtomType atomType;
    
    /**
     * Makes an atom factory with atoms having AtomSequencerSimple and
     * AtomTreeNodeGroup for sequencer and node, respectively.
     * @param space
     */
    public AtomFactory(Space space) {
        this(space, AtomSequencerFactory.SIMPLE);
    }
    
    public AtomFactory(Space space, AtomSequencerFactory sequencerFactory) {
    	this(space, sequencerFactory, AtomTreeNodeGroup.FACTORY);
    }
    
    public AtomFactory(Space space, AtomSequencerFactory sequencerFactory, AtomTreeNode.Factory nodeFactory) {
        this.space = space;
        this.sequencerFactory = sequencerFactory;
        this.nodeFactory = nodeFactory;
        this.atomType = groupType;
        reservoir = new AtomReservoir(this);
    }
    
	public AtomType getType() {
		return atomType;
	}
	
    /**
     * Makes an atom with the reservoir as its parent.
     */
    public Atom makeAtom() {
        return makeAtom((AtomTreeNodeGroup)reservoir.node);
    }
    
    /**
     * Makes an atom with the given node as its parent.  If reservoir
     * is not empty, takes atom from it; otherwise makes it from scratch
     * using the build method.
     */
    public Atom makeAtom(AtomTreeNodeGroup parent) {
        Atom atom = reservoir.getAtom();
        if(atom == null) atom = build(parent);//reservoir);//using reservoir causes problem for BravaisLattice
        atom.node.setParent(parent);
        
        //add agents from any registered sources
        if(agentSource.length > 0) atom.agents = new Object[agentSource.length];
        for(int i=0; i<agentSource.length; i++) {
            atom.agents[i] = agentSource[i].makeAgent(atom);
        }
        
        return atom;
    }
    
    /**
     * Constructs a new atomgroup having the given parent and sends it
     * to the build(Atom) method for addition of children.  This method 
     * is sometime overridden in subclasses to construct with atom that does
     * not use the group AtomTreeNode, or that is a subclass of Atom.
     */
    protected Atom build(AtomTreeNodeGroup parent) {
        Atom group = new Atom(space, atomType, nodeFactory, sequencerFactory, parent);
        return build(group);
    }
    
//    public static Atom build(Space space, AtomTreeNodeGroup parent, Model model) {
//    	if(model instanceof ModelAtomic) {
//    		return new Atom(space, atomType, nodeFactory, sequencerFactory, parent);
//    	}
//		return null;
//    }
    
    /**
     * Builds an atom from the one given, attaching child atoms as appropriate
     * for the definition of the concrete subclass.
     */
    public abstract Atom build(Atom atom);
    
    /**
     * Indicates if this factory produces atom groups or simple atoms.
     * If this method returns true, then the atoms made by this factory will have
     * a node of type (or derived from) AtomTreeNodeGroup.
     */
    public abstract boolean isGroupFactory();
        
    public AtomReservoir reservoir() {return reservoir;}
    
    public Space space() {return space;}
    
    public void setSpecies(Species species) {this.species = species;}
    public Species species() {return species;}
        
    public void setConfiguration(Configuration config) {configuration = config;}
    public Configuration getConfiguration() {return configuration;}
    
    /**
     * Adds given agent source to agent-source array and returns index
     * indicating where in atom agent-array the source's agent will
     * be placed.
     */
    public int requestAgentIndex(Atom.AgentSource aSource) {
        Atom.AgentSource[] newSource = new Atom.AgentSource[agentSource.length+1];
        for(int i=0; i<agentSource.length; i++) newSource[i] = agentSource[i];
        int index = agentSource.length;
        newSource[index] = aSource;
        agentSource = newSource;
        return index;
    }
}//end of AtomFactory