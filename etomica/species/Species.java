package etomica.species;
import etomica.atom.AtomFactory;
import etomica.atom.AtomManager;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSpeciesAgent;
import etomica.atom.ISpeciesAgent;
import etomica.atom.SpeciesAgent;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;


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
  * @see AtomManager
  * @see SpeciesAgent
  * @see PotentialMaster
  */
public class Species implements java.io.Serializable {

    /**
     * Constructs species with molecules built by the given atom factory.
     * Species agents made by this species will have the given type for 
     * their (common) AtomType.
     */
    public Species(AtomFactory factory) {
        this.factory = factory;
    }
    
    /**
     * Sets the AtomType for SpeciesAgents associated with this Species.
     * This method should only be called by the SpeciesManager
     */
    public void setAgentType(AtomTypeSpeciesAgent agentType) {
        this.agentType = agentType;
        factory.getType().setParentType(agentType);
    }

    public AtomFactory moleculeFactory() {return factory;}
    
    /**
     * Constructs an Agent of this species and sets its SpeciesMaster.
     * The agent's type is in common with all other agents of this species.
     * 
     * @param atomManager The SpeciesMaster that will hold this SpeciesAgent
     * @return The new agent.
     */
    public ISpeciesAgent makeAgent(AtomManager atomManager) {
        SpeciesAgent agent = new SpeciesAgent(agentType, atomManager);
        atomManager.addSpeciesAgent(agent);
        return agent;
    }

    /**
     * Returns the agent of this species in the given phase
     * Hashmap is used to connect phase(key)-agent(value) pairs
     * 
     * @param p The phase for which this species' agent is requested
     * @return The agent of this species in the phase
     */
    public final ISpeciesAgent getAgent(Phase p) {return p.getAgent(this);}

    public final AtomType getMoleculeType() {
        return factory.getType();
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
//    final AtomType agentType;
    protected final AtomFactory factory;
    protected AtomTypeSpeciesAgent agentType;
}
