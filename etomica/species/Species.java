package etomica.species;
import etomica.atom.AtomFactory;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeGroup;
import etomica.atom.SpeciesAgent;
import etomica.atom.SpeciesMaster;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.util.NameMaker;


 /**
  * A Species holds information about how to construct a molecule, and
  * provides for the management of molecules in a phase.
  * 
  * These are the important features of a Species: <br>
  * <ol>
  * <li>It holds an AtomFactory instance that constructs molecules when
  * needed.
  * <li>It makes a SpeciesAgent class that is placed in each phase (each
  * phase has one species agent from each species). This agent manages the
  * molecules of that species in that phase (counting, adding, removing,
  * etc.) The agent for a given phase may be obtained through the getAgent
  * method. <br>
  * <li>Each Species has a unique species index assigned when it is
  * constructed. The index assignment begins at 1 and is incremented after
  * each Species construction. This index is useful when collecting things
  * in reference to the species (for example, in the use of neighbor lists).
  * </ol>
  * The number of molecules of a species in a phase may be changed at run
  * time. Interactions among all molecules in a phase are defined by
  * associating an intermolecular potential to one or more Species via a
  * call to the addPotential method of the PotentialMaster for the simulation.
  * 
  * @author David Kofke
  * @author C. Daniel Barnes
  * @see SpeciesMaster
  * @see SpeciesAgent
  * @see PotentialMaster
  */
public class Species implements java.io.Serializable {

    /**
     * Constructs species with molecules built by the given atom factory.
     * Species agents made by this species will have the given type for 
     * their (common) AtomType.
     */
    public Species(AtomFactory factory, AtomType agentType) {
        this.factory = factory;
        this.agentType = agentType;
        setName(NameMaker.makeName(this.getClass()));
        agentType.setSpecies(this);
    }
    
    /**
     * Accessor method of the name of this species.
     * 
     * @return The given name of this species
     */
    public final String getName() {return name;}
    
    /**
     * Method to set the name of this species. The species' name
     * provides a convenient way to label output data that is associated with
     * it.  This method might be used, for example, to place a heading on a
     * column of data. Default name is the base class followed by the integer
     * index of this element.
     * 
     * @param name The name string to be associated with this species.
     */
    public void setName(String name) {this.name = name;}

    /**
     * Overrides the Object class toString method to have it return the output of getName.
     * 
     * @return The name given to the phase
     */
    public String toString() {return getName();}
    
    public AtomFactory moleculeFactory() {return factory;}
    
    /**
     * Constructs an Agent of this species and sets its parent phase.
     * The agent's type is in common with all other agents of this species.
     * 
     * @param p The given parent phase of the agent
     * @return The new agent.
     */
    public SpeciesAgent makeAgent(SpeciesMaster parent) {
        SpeciesAgent agent = new SpeciesAgent(agentType, this);
        agent.getNode().setParent((AtomTreeNodeGroup)parent.getNode());
        return agent;
    }

    /**
     * Returns the agent of this species in the given phase
     * Hashmap is used to connect phase(key)-agent(value) pairs
     * 
     * @param p The phase for which this species' agent is requested
     * @return The agent of this species in the phase
     */
    public final SpeciesAgent getAgent(Phase p) {return p.getAgent(this);}

    public final AtomType getMoleculeType() {
        return factory.getType();
    }
    
    public final AtomFactory getFactory() {
        return factory;
    }
    
    /**
     * Returns an AtomType that is appropriate for passing to the constructor.
     * Method is used by subclasses to generate the AtomType for the agent.
     * This is needed for it to make index managers for the AtomTypes needed
     * to make the factory (which is also passed to the Species constructor).
     * PositionDefinition for the agent type is set to null.
     */
    public static AtomTypeGroup makeAgentType(Simulation sim) {
        AtomTypeGroup agentType = new AtomTypeGroup(null);
        agentType.setParentType(sim.getSpeciesRoot().getSpeciesMasterType());
        return agentType;
    }
    
    /**
     * Returns a SpeciesSignature for this Species.  Subclasses must override
     * this method.
     */
    public SpeciesSignature getSpeciesSignature() {
        // AtomFactories can't be serialized (or rather, they'll serialize too 
        // many other things)
        return null;
    }
    
    private static final long serialVersionUID = 2L;
    final AtomType agentType;
    protected final AtomFactory factory;
    private String name;
}
