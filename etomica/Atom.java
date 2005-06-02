package etomica;

import etomica.atom.AtomLinker;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomTypeSphere;
import etomica.space.ICoordinate;
import etomica.utility.Arrays;

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
public class Atom implements AtomSet, Comparable, java.io.Serializable {

    public Atom(Space space, AtomType type, 
                    AtomTreeNodeFactory nodeFactory,
                    AtomSequencerFactory seqFactory) {
        this.type = type;//do this first
        seq = seqFactory.makeSequencer(this);
        node = nodeFactory.makeNode(this);
        coord = space.makeCoordinate(this);//must follow setting of type field
        node.setOrdinal(0,0); //-(++INSTANCE_COUNT));//default index; changed when added to parent after construction
        
        if(agentSource.length > 0) allatomAgents = new Object[agentSource.length];
        for(int i=0; i<agentSource.length; i++) {
            allatomAgents[i] = agentSource[i].makeAgent(this);
        }
    }
    
    /**
     * Makes a simple atom for the given space.  Node is for a leaf atom;
     * sequencer is simple; type is a sphere, unique to the new atom; depth is 0.
     * @param space
     */
    public Atom(Space space) {
    	this(space, new AtomTypeSphere(null), AtomTreeNodeLeaf.FACTORY, AtomSequencerFactory.SIMPLE);                        
        node.setOrdinal(0,++INSTANCE_COUNT);//default index; changed when added to parent after construction
    }
    
    public final int count() {return 1;}
    
    public final Atom getAtom(int i) {
        if(i == 0 ) return this;
        throw new IllegalArgumentException();
    }
    
    public boolean equals(Object object) {
        if(!(object instanceof AtomSet) || ((AtomSet)object).count() != 1) return false;
        return this == ((AtomSet)object).getAtom(0);
    }
    
    public final boolean equals(AtomSet atoms) {
        if (atoms.count() != 1) return false;
        return this == atoms.getAtom(0);
    }
    
    public boolean inSameSpecies(Atom atom) {
        return type.getIndexManager().sameSpecies(node.index(), atom.node.index());
    }
    public boolean inSameMolecule(Atom atom) {
        return type.getIndexManager().sameMolecule(node.index(), atom.node.index());
    }
    public boolean inSamePhase(Atom atom) {
        return type.getIndexManager().samePhase(node.index(), atom.node.index());
    }
    
    /**
     * Returns a string of digits that uniquely identifies this atom.  String is
     * formed by concatenating the index of this atom to the signature
     * given by the parent of this atom.
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
            
    public Object ia;//integrator agent
                
    public final AtomTreeNode node;
        
    public final AtomType type;
    
    public final AtomLinker seq;
    
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