package etomica.atom;

import etomica.Atom;
import etomica.Configuration;
import etomica.Space;
import etomica.Species;

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
    
    public final Space space;
    private Species species;
    protected Configuration configuration;
    protected AtomSequencerFactory sequencerFactory;
    protected AtomTreeNodeFactory nodeFactory;
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
    
    public AtomFactory(Space space, AtomSequencerFactory sequencerFactory, AtomTreeNodeFactory nodeFactory) {
        this.space = space;
        this.sequencerFactory = sequencerFactory;
        this.nodeFactory = nodeFactory;
        this.atomType = groupType;
    }
    
	public AtomType getType() {
		return atomType;
	}
	
    /**
     * Builds and returns the atom/atomgroup made by this factory.
     * Implementation of this method in the subclass defines the 
     * product of this factory.
     */
    public abstract Atom makeAtom();
    
    /**
     * Method used by subclasses to make the root atom of the group it is building.
     */
    protected Atom newParentAtom() {
        Atom atom = new Atom(space, atomType, nodeFactory, sequencerFactory);
        //add agents from any registered sources
        if(agentSource.length > 0) atom.agents = new Object[agentSource.length];
        for(int i=0; i<agentSource.length; i++) {
            atom.agents[i] = agentSource[i].makeAgent(atom);
        }
        return atom;
    }
    
    /**
     * Indicates if this factory produces atom groups or simple atoms.
     * If this method returns true, then the atoms made by this factory will have
     * a node of type (or derived from) AtomTreeNodeGroup.
     */
    public abstract boolean isGroupFactory();
        
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