package etomica;

 /**
  * Object corresponding to one physical atom or (in subclasses) group of atoms.
  * Each atom holds one instance of a Coordinate that is constructed by the governing
  * space class. 
  * all simulation kinetics and dynamics are performed by operating on the atom's coordinate.
  * In addition, an Atom has a type that holds information needed to draw the atom,
  * and possibly information to be used as parameters in a potential.  Each atom
  * holds an Integrator.Agent object that may be used to store information needed
  * by the integrator.
  * 
  * @author David Kofke
  * @author C. Daniel Barnes
  */
public class Atom implements java.io.Serializable {

  /*  public Atom(Space space, AtomType type, AtomTreeNodeGroup parent) {
        this(space, type, AtomTreeNodeGroup.FACTORY, parent);
    }
    public Atom(Space space, AtomType type, 
                    AtomTreeNode.Factory nodeFactory, AtomTreeNodeGroup parent) {
        this(space, type, nodeFactory, IteratorFactorySimple.INSTANCE.simpleSequencerFactory(), parent);
    }*/
    public Atom(Space space, AtomType type, 
                    AtomTreeNode.Factory nodeFactory,
                    AtomSequencer.Factory seqFactory, AtomTreeNodeGroup parent) {
        this.type = type;//do this first
        seq = seqFactory.makeSequencer(this);
        node = nodeFactory.makeNode(this, parent);//must follow setting of sequencer to permit addition to childlist of parent
                                                  //setting parent here in constructor bypasses the call to the parent's addAtomNotify method (which would be called if parent were set in setParent method)
        coord = space.makeCoordinate(this);//must follow setting of type field
        if(parent != null) {
        	parent.addAtomNotify(this);
        } else {
        	node.setIndex(NO_PARENT_INSTANCE_COUNT++);
        }
        
        if(agentSource.length > 0) allatomAgents = new Object[agentSource.length];
        for(int i=0; i<agentSource.length; i++) {
            allatomAgents[i] = agentSource[i].makeAgent(this);
        }
        seq.setParentNotify(parent);//added 08/11/03 (DAK)
    }
    
    /**
     * Makes a simple atom for the given space.  Node is for a leaf atom;
     * sequencer is simple; type is a sphere, unique to the new atom.
     * @param space
     */
    public Atom(Space space) {
    	this(space, new AtomType.Sphere(null), AtomTreeNodeLeaf.FACTORY, AtomSequencerSimple.FACTORY, null);                        
    }
    
    /**
     * Returns a string of digits that uniquely identifies this atom.  String is
     * formed by concatenating the index of this atom to the signature
     * given by the parent of this atom.
     */
    public String signature() {return node.parentGroup().signature() + " " + node.index();}

    /**
     * Returns a string formed by concatenating the signature of this atom
     * to a string that identifies it as a species master, species agent, 
     * molecule, group, or (leaf) atom.
     */
    public final String toString() {
    	if(this instanceof SpeciesMaster) return "Master(" + signature() + ")";
    	else if(this instanceof SpeciesAgent) return "Agent(" + signature() + ")";
    	if(node.parentGroup() instanceof SpeciesAgent) return "Molecule(" + signature() + ")";
    	else if(node.isLeaf()) return "Atom(" + signature() + ")";
    	else return "Group(" + signature() + ")";
    }    


    /**
     * Assigns the atom's integrator agent to the given instance.
     */
    public void setIntegratorAgent(Object ia) {this.ia = ia;}
            
	/**
	 *  Returns true if this is the given atom.  Part of AtomSet interface.
	 */
	public boolean contains(Atom a) {return this == a;}
	
	/**
	 * Returns 1, indicating this is an atom set of one atom. Part of
	 * AtomSet interface.
	 */
	public final int nBody() {return 1;}

//   linked lists of bonds
    public BondLinker firstUpBond;
    public BondLinker firstDownBond;
    
    public void sendToReservoir() {
        creator().reservoir().addAtom(this);
    }
    
    public AtomFactory creator() {return type.creator();}
    
    /**
     * Coordinates of this atom.
     * When the atom is constructed the coordinate class is provided by the 
     * governing Space for the simulation.
     */
    public final Space.Coordinate coord;
            
    public Object ia;//integrator agent
                
    public final AtomTreeNode node;
        
    public final AtomType type;
    
    public final AtomSequencer seq;
    
    /**
     * An array of agents, or objects made by an agent source and added to this
     * atom (and all others produced by this atom's creator factory) 
     * to perform some function or store data for the source.  
     * Placement of agents in this array is managed by the factory
     * producing this atom.  An agent source registers itself with the factory
     * via the requestAgentIndex method, and the integer returned with that call
     * indicates the place in this array where the source can find its agent
     * in this atom.
     */
    public Object[] agents;
    
    /**
     * An array of agents, or objects made by an agent source and added to this
     * atom (and all other atoms) to perform some function or store data for the source.
     * Placement of agents in this array is managed by the Atom class.  An agent
     * source registers itself with Atom via the Atom.requestAgentIndex method, and
     * the integer returned with that call indicates the place in this array where
     * the source can find its agent in this atom.
     * Differs from agents class in that all atoms get an instance of an agent if put
     * in the allatomAgents array, while only atoms produced by a particular factory
     * are given agents that are in the agents array.
     */
    public Object[] allatomAgents;
    
    /**
     * Array used to record all agent sources requesting to place an agent in every atom.
     */
    private static AgentSource[] agentSource = new AgentSource[0];
    
    /**
     * Counter for number of times an atom is instantiated without a parent.  Used
     * to assign a unique index to such atoms.
     */
    private static int NO_PARENT_INSTANCE_COUNT = 0;
    
    /**
     * Adds given agent source to allatomAgent-source array and returns index
     * indicating where in atom's allatomAgent-array the source's agent will
     * be placed.
     */
    public static int requestAgentIndex(AgentSource aSource) {
        AgentSource[] newSource = new AgentSource[agentSource.length+1];
        for(int i=0; i<agentSource.length; i++) newSource[i] = agentSource[i];
        int index = agentSource.length;
        newSource[index] = aSource;
        agentSource = newSource;
        return index;
    }    
    
    /**
     * Interface for an object that makes an agent to be placed in each atom
     * upon construction.  AgentSource objects register with the AtomFactory
     * the produces the atom.
     */
    public interface AgentSource {
        public Object makeAgent(Atom a);
    }
    
}//end of Atom