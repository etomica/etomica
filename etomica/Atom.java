package etomica;

import etomica.atom.AtomLinker;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomTypeSphere;
import etomica.space.Coordinate;
import etomica.space.CoordinateFactory;
import etomica.space.CoordinateFactorySphere;
import etomica.space.ICoordinate;
import etomica.utility.Arrays;

 /**
  * Object corresponding to one physical atom or group of atoms. Each atom holds
  * the following publicly accessible fields:
  * <ul>
  * <li>a Coordinate instance (fieldname: coord) that is constructed by the
  * governing space class; the coordinate stores information about the state of
  * the atom -- usually its position and momentum, but other definitions are possible
  * <li>an AtomType instance (fieldname: type) that holds information this atom
  * has in common with other atoms made by the same factory
  * <li>an instance of AtomTreeNode (fieldname: node) that is used to place it
  * in the species hierarchy
  * <li>an instance of AtomLinker (fieldname: seq, for sequencer) that places
  * it in the linked list of children held by its parenttom holds an object
  * <li>an object (fieldname: ia, for integrator agent) that may be used to
  * store information needed by the integrator
  * <li>an array of objects (field name: allAtomAgents) that can be used to
  * store in each atom any object needed by another class; such a class must
  * implement Atom.AgentSource and request its object be stored in every atom by
  * invoking Atom.requestAgentIndex before any atoms are constructed. The
  * integer returned by this method will indicate the location in the
  * allAtomAgents array where the agent-source's object will be held.
  * </ul>
  * @author David Kofke and C. Daniel Barnes
  * 
  * @see Coordinate
  * @see AtomType
  * @see AtomTreeNode
  */
public class Atom implements AtomSet, Comparable, java.io.Serializable {

    public Atom(ICoordinate coord, AtomType type, 
                    AtomTreeNodeFactory nodeFactory,
                    AtomSequencerFactory seqFactory) {
        this.type = type;
        this.coord = coord;
        seq = seqFactory.makeSequencer(this);
        node = nodeFactory.makeNode(this);
        node.setOrdinal(0,0); //-(++INSTANCE_COUNT));//default index; changed when added to parent after construction
        
        if(agentSource.length > 0) allatomAgents = new Object[agentSource.length];
        for(int i=0; i<agentSource.length; i++) {
            allatomAgents[i] = agentSource[i].makeAgent(this);
        }
    }
    
    /**
     * Makes a simple atom for the given space.  Coordinate is non-kinetic sphere;
     * node is for a leaf atom; sequencer is simple; type is a sphere, unique to the new atom; 
     * depth is 0.
     */
    public Atom(Space space) {
    	this(new CoordinateFactorySphere(space, false).makeCoordinate(), new AtomTypeSphere(null), AtomTreeNodeLeaf.FACTORY, AtomLinker.FACTORY);                        
        node.setOrdinal(0,++INSTANCE_COUNT);//default index; changed when added to parent after construction
    }
    
    /**
     * Returns 1, indicating that this AtomSet is an Atom.
     */
    public final int count() {return 1;}
    
    /**
     * Returns this if i==0, otherwise throws exception.
     * 
     * @throws IllegalArgumentException if i != 0
     */
    public final Atom getAtom(int i) {
        if(i == 0 ) return this;
        throw new IllegalArgumentException();
    }
    
    /**
     * Returns true if the given object is this atom instance, or if it is
     * a length-1 AtomSet holding this instance.
     */
    public boolean equals(Object object) {
        if(!(object instanceof AtomSet) || ((AtomSet)object).count() != 1) return false;
        return this == ((AtomSet)object).getAtom(0);
    }

    /**
     * Returns true if the given object is this atom instance, or if it is
     * a length-1 AtomSet holding this instance.
     */
    public final boolean equals(AtomSet atoms) {
        if (atoms == null || atoms.count() != 1) return false;
        return this == atoms.getAtom(0);
    }

    /**
     * Returns true if this atom is in the same species as the given atom.
     * 
     * @throws NullPointerException
     *             if the argument is null
     */
    public boolean inSameSpecies(Atom atom) {
        return type.getIndexManager().sameSpecies(node.index(), atom.node.index());
    }
    /**
     * Returns true if this atom is in the same molecule as the given atom.
     * 
     * @throws NullPointerException
     *             if the argument is null
     */
    public boolean inSameMolecule(Atom atom) {
        return type.getIndexManager().sameMolecule(node.index(), atom.node.index());
    }
    /**
     * Returns true if this atoms is in the same phase as the given atom.
     * 
     * @throws NullPointerException
     *             if the argument is null
     */
    public boolean inSamePhase(Atom atom) {
        return type.getIndexManager().samePhase(node.index(), atom.node.index());
    }
    
    /**
     * Returns a string of digits that uniquely identifies this atom.  String is
     * formed by concatenating the ordinal of this atom to the signature
     * given by the parent of this atom.  If atom has no parent, forms a string
     * from only the ordinal.
     */
    public String signature() {
        if(node.parentGroup() != null) {
            return node.parentGroup().signature() + " " + node.getOrdinal();
        }
        return Integer.toString(node.getOrdinal());
    }
    /**
     * Returns a string formed by concatenating the signature of this atom
     * to a string that identifies it as a species master, species agent, 
     * molecule, group, or (leaf) atom.
     */
    public final String toString() {
//        return Integer.toBinaryString(node.index());
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
     * Implementation of Comparable interface.  Returns -1, 0, 1 if given atom
     * is less, equal, or greater, respectively, than this atom.  
     * Order is determined by compareTo method of atoms' nodes.
     */
    public int compareTo(Object atom) {
        return node.compareTo(((Atom)atom).node);
    }
    
    /**
     * Coordinates of this atom.
     * When the atom is constructed the coordinate class is provided by the 
     * governing Space for the simulation.
     */
    public final ICoordinate coord;
    
    /**
     * Integrator agent, used by integrator to implement the simulation integration algorithm.
     */
    public Object ia;//integrator agent
    
    /**
     * Tree node, used to place the atom in the species tree.
     */
    public final AtomTreeNode node;
    
    /**
     * Atom type, holding properties held in common with other atoms made by this atom's
     * factory.
     */
    public final AtomType type;
    
    /**
     * Sequencer, used to access directly this atom's place in the childlist of its parent.
     */
    public final AtomLinker seq;
    
    /**
     * An array of agents, or objects made by an agent source and added to this
     * atom (and all other atoms) to perform some function or store data for the source.
     * Placement of agents in this array is managed by the Atom class.  An agent
     * source registers itself with Atom via the Atom.requestAgentIndex method, and
     * the integer returned with that call indicates the place in this array where
     * the source can find its agent in this atom.
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
    private static int INSTANCE_COUNT = 0;
    
    /**
     * Adds given agent source to allatomAgent-source array and returns index
     * indicating where in atom's allatomAgent-array the source's agent will
     * be placed.
     */
    public static int requestAgentIndex(AgentSource aSource) {
        agentSource = (AgentSource[])Arrays.addObject(agentSource, aSource);
        return agentSource.length - 1;
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