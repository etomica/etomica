package etomica;

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
    public final Simulation simulation;
    private Species species;
    protected Configuration configuration;
    protected AtomSequencer.Factory sequencerFactory;
    protected AtomTreeNodeGroup.Factory nodeFactory;
    protected BondInitializer bondInitializer = BondInitializer.NULL;
    private Atom.AgentSource[] agentSource = new Atom.AgentSource[0];
    protected final AtomType.Group groupType = new AtomType.Group(this);
    protected AtomType atomType;
    
    public AtomFactory(Simulation sim, AtomSequencer.Factory sequencerFactory) {
    	this(sim, sequencerFactory, AtomTreeNodeGroup.FACTORY);
    }

	public AtomFactory(Simulation sim, AtomSequencer.Factory sequencerFactory, AtomTreeNode.Factory nodeFactory) {    
    	this.simulation = sim;
    	this.space = sim.space;
    	this.sequencerFactory = sequencerFactory;
    	this.nodeFactory = nodeFactory;
    	this.atomType = groupType;
    	reservoir = new AtomReservoir(this);
    }
    
    public AtomFactory(Space space, AtomSequencer.Factory sequencerFactory) {
    	this(space, sequencerFactory, AtomTreeNodeGroup.FACTORY);
    }
    
    public AtomFactory(Space space, AtomSequencer.Factory sequencerFactory, AtomTreeNode.Factory nodeFactory) {
    	simulation = null;
        this.space = space;
        this.sequencerFactory = sequencerFactory;
        this.nodeFactory = nodeFactory;
        this.atomType = groupType;
        reservoir = new AtomReservoir(this);
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
    
//    protected abstract void renew(Atom a);
    
//    public abstract boolean vetoAddition(Atom a); //be sure to check that a is non-null
    
    public AtomReservoir reservoir() {return reservoir;}
    
    public Space space() {return space;}
    
    public Simulation simulation() {return simulation;}
    
    public void setSpecies(Species species) {this.species = species;}
    public Species species() {return species;}
        
    public void setConfiguration(Configuration config) {configuration = config;}
    public Configuration getConfiguration() {return configuration;}
    
    public void setBondInitializer(BondInitializer bonder) {bondInitializer = bonder;}
    public BondInitializer getBondInitializer() {return bondInitializer;}
    
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