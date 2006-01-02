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
import etomica.units.Dimension;
import etomica.units.Quantity;
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
 
public abstract class Species implements Comparable, java.io.Serializable {

    /**
     * Constructs species with molecules built by the given atom factory.
     * Species agents made by this species will have the given type for 
     * their (common) AtomType.
     */
    public Species(Simulation sim, AtomFactory factory, AtomType agentType) {
        this.factory = factory;
        this.agentType = agentType;
        nMolecules = sim.getDefaults().moleculeCount;
        setName(NameMaker.makeName(this.getClass()));
        index = sim.speciesRoot.addSpecies(this);
        factory.setSpecies(this);
    }
    
    /**
     * Constructs species with an atom factory that makes molecules from
     * the given model for the given space.
     */
//    public Species(Simulation sim, Model model) {
//    	this(sim, model.makeAtomFactory(sim.space));
//    }

    public int getIndex() {
        return index;
    }
    
    /**
     * Implementation of Comparable interface according to the relative order of this and the
     * given species' indexes.  If this index is greater, returns +1; if they are equal, returns 0;
     * if the given species index is greater, returns -1.
     * 
     * @throws ClassCastException if the given object is not an instance of Species
     */
    public int compareTo(Object species) {
        return ((Species)species).index > index ? -1 : (((Species)species).index == index ? 0 : +1);       
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
     * Nominal number of molecules of this species in each phase.
     * Actual number may differ if molecules have been added or removed to/from the phase
     */
    protected int nMolecules;;
    
    /**
     * Accessor method for nominal number of molecules in each phase.  Actual 
     * number of molecules of this species in a given phase is obtained via
     * the getNMolecules method of Species.Agent
     * 
     * @return Nominal number of molecules in each phase
     * @see SpeciesAgent#getNMolecules
     */
    public int getNMolecules() {return nMolecules;}
    public Dimension getNMoleculesDimension() {return Quantity.DIMENSION;}
    
    /**
     * Sets the number of molecules of this species for each phase.
     * Propagates the change to agents of this species in all phases, so
     * it creates the given number of molecules in every phase.
     * 
     * @param n The new number of molecules of this species in each phase
     * @see SpeciesAgent#setNMolecules
     */
    public void setNMolecules(int n) {
        nMolecules = n;
    }
    
    /**
     * Constructs an Agent of this species and sets its parent phase.
     * The agent's type is in common with all other agents of this species.
     * 
     * @param p The given parent phase of the agent
     * @return The new agent.
     */
    public SpeciesAgent makeAgent(SpeciesMaster parent) {
        SpeciesAgent agent = new SpeciesAgent(agentType, this);
        agent.node.setParent((AtomTreeNodeGroup)parent.node);
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
        return new AtomTypeGroup(sim.speciesRoot.getChildType(), null);
    }
    
    public abstract SpeciesSignature getSpeciesSignature();
    
    final AtomType agentType;
    protected final AtomFactory factory;
    private String name;
    private final int index;
    
}
