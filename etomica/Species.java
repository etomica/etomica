package etomica;
import etomica.units.Dimension;
import etomica.chem.*;
import etomica.chem.models.*;
//Java2 imports
//import java.util.HashMap;
//import java.util.Iterator;

//import etomica.utility.HashMap;
//import etomica.utility.Iterator;


     /**
      *
      * This description is out-of-date.
      *
      * A Species is a collection of identically formed Molecules.  Subclasses of
      * Species differ in the type of molecules that they collect.  Subclassing
      * Species is the primary way to define new types of molecules (Molecule is
      * subclassed only to define molecules formed from unusual (e.g., anisotropic)
      * atoms).  The type of a molecule is determined by the AtomType array passed 
      * to its constructor, and by the ConfigurationMolecule object that is used to 
      * define the internal arrangement of the atoms in the molecule.  Species differ
      * in how they define these classes when constructing molecules.
      * 
      * These are the important features of a Species:<br>
      * 
      * 1. It holds the definition of the AtomType and ConfigurationMolecule objects that
      * determine the type of molecule defined by the species<br>
      * 
      * 2. It makes a Species.Agent class that is placed in each phase (each phase has one
      * species agent from each species).  This agent manages the molecules of that species
      * in that phase (counting, adding, removing, etc.)  More information about this class 
      * is provided in its documentation comments.<br>
      * 
      * 3. It maintains a reservoir of molecules so that molecules can be added and removed 
      * to/from a phase without requiring (expensive) construction of a new molecule each time.
      * Correspondingly, it provides a makeMolecule method that can be used to obtain a molecule
      * of the species type.  This molecule is taken from the reservoir if available, or made fresh
      * if the reservoir is empty.<br>
      * 
      * 4. Each Species is part of a linked list of species objects.<br> 
      * 
      * 5. Each Species has a unique species index assigned when it is registered with the
      * Simulation.  This index is used to determine the potential for the interaction of 
      * atoms of the species.<br>
      * 
      * The number of Molecules of a Species in a Phase may be changed at run time.
      * One or more Species are collected together to form a Phase.
      * Interactions among all molecules in a Phase are defined by associating an
      * intermolecular potential to pairs of Species (or to a single Species to
      * define interactions among its molecules). Intramolecular potentials are
      * also associated with a Species, and are used to describe intra-molecular 
      * atomic interactions.
      * 
      * @author David Kofke
      * @author C. Daniel Barnes
      * @see SpeciesAgent
      */
     
public class Species extends SimulationElement {

    public static final String VERSION = "Species:03.01.25";
    
    protected final AtomFactory factory;
    protected Model model;
    
    /**
     * Constructs species and registers it as part of the given simulation, with
     * molecules built by the given atom factory.
     */
    public Species(Simulation sim, AtomFactory factory) {
        super(sim, Species.class);
        this.factory = factory;
    }
    
    public Species(Simulation sim, Model model) {
    	this(sim, model.makeAtomFactory(sim));
    }
                  
    public AtomFactory moleculeFactory() {return factory;}
    
    /**
     * Returns an iterator of this species (leaf) atoms in the given phase.
     */
/*    public AtomIterator makeAtomIterator(Phase p) {
        return getAgent(p).new LeafAtomIterator();
    }
    /**
     * Returns an iterator of this species agent's children (molecules) in the given phase.
     * /
    public AtomIterator makeMoleculeIterator(Phase p) {
        return getAgent(p).new ChildAtomIterator();
    }
*/        
    /**
     * Nominal number of molecules of this species in each phase.
     * Actual number may differ if molecules have been added or removed to/from the phase
     */
    protected int nMolecules = 20;
    
    /**
     * Accessor method for nominal number of molecules in each phase.  Actual 
     * number of molecules of this species in a given phase is obtained via
     * the getNMolecules method of Species.Agent
     * 
     * @return Nominal number of molecules in each phase
     * @see Species.Agent#getNMolecules
     */
    public int getNMolecules() {return nMolecules;}
    public Dimension getNMoleculesDimension() {return Dimension.QUANTITY;}
    
    /**
     * Sets the number of molecules of this species for each phase.
     * Propagates the change to agents of this species in all phases, so
     * it creates the given number of molecules in every phase.
     * 
     * @param n The new number of molecules of this species in each phase
     * @see Species.Agent#setNMolecules
     */
    public void setNMolecules(int n) {
        nMolecules = n;
        allAgents(new AtomAction() {public void actionPerformed(Atom a) {((SpeciesAgent)a).setNMolecules(nMolecules);}});
    }
    
    public Atom makeMolecule() {
        return factory.makeAtom();
    }
    
    /**
     * Performs the given action on all of this species agents in all phases.
     */
    public void allAgents(AtomAction action) {
        if(action == null) return;
        agents.doToAll(action);
        /*
        Iterator e = agents.values().iterator();
        while(e.hasNext()) {
            action.actionPerformed((Atom)e.next());
        }*/
    }
    
    /**
     * Performs the given action on all of this species molecules in all phases.
     */
    public void allMolecules(AtomAction action) {
        if(action == null) return;
        allMoleculeAction.setAction(action);//copy action to wrapper class
        agents.doToAll(allMoleculeAction);  //apply action wrapper over all agents
    }
    
    /**
     * Performs the given action on all of this species (leaf) atoms in all phases.
     */
    public void allAtoms(AtomAction action) {
        if(action == null) return;
        allAtomAction.setAction(action);//copy action to wrapper class
        agents.doToAll(allAtomAction);  //apply action wrapper over all agents
    }
    
    /**
     * Constructs an Agent of this species and sets its parent phase
     * 
     * @param p The given parent phase of the agent
     * @return The new agent.  Normally this is returned into the setAgent method of the given phase.
     * @see Species.Agent
     */
    public SpeciesAgent makeAgent(SpeciesMaster parent) {
        Phase phase = parent.node.parentPhase();
        SpeciesAgent a = new SpeciesAgent(this, nMolecules, parent.node);
        agents.put(phase,a);   //associate agent with phase; retrieve agent for a given phase using agents.get(p)
        return a;
    }

    /**
     * Returns the agent of this species in the given phase
     * Hashmap is used to connect phase(key)-agent(value) pairs
     * 
     * @param p The phase for which this species' agent is requested
     * @return The agent of this species in the phase
     */
    public final SpeciesAgent getAgent(Phase p) {return agents.get(p);}
    public final SpeciesAgent getAgent(SpeciesMaster s) {return agents.get(s);}
        

    /**
     * Hashtable to associate agents with phases
     */
 //   final HashMap agents = new HashMap();
    final AgentList agents = new AgentList();
    
    /**
     * Class that keeps a list of all agents in a way that
     * they can be referenced according to the phase they are in.
     * Uses the phase index to index them.
     * Mimics hash functionality.
     */
    private static final class AgentList {
        
        private SpeciesAgent[] agentArray = new SpeciesAgent[0];
        
        void put(Phase phase, SpeciesAgent agent) {
            int index = phase.index;
            //expand array size
            if(index >= agentArray.length) {
                SpeciesAgent[] newArray = new SpeciesAgent[index+1];
                for(int i=0; i<agentArray.length; i++) newArray[i] = agentArray[i];
                agentArray = newArray;
            }
            agentArray[index] = agent;
        }
        
        SpeciesAgent get(Phase phase) {return agentArray[phase.index];}
        SpeciesAgent get(SpeciesMaster s) {return agentArray[s.index];}
        
        void doToAll(AtomAction action) {
            for(int i=agentArray.length-1; i>=0; i--) action.actionPerformed(agentArray[i]);
        }
        
    }//end of AgentList
    
    /**
     * Wrapper that enables an action to be performed on all molecules under a species agent.
     */
    private final AllMoleculeWrapper allMoleculeAction = new AllMoleculeWrapper();
    private final class AllMoleculeWrapper extends AtomAction {
        private AtomAction action;
        void setAction(AtomAction a) {this.action = a;}
        public void actionPerformed(Atom a) {((SpeciesAgent)a).allMolecules(action);}
    }

    /**
     * Wrapper that enables an action to be performed on all leaf atoms under a species agent.
     */
    private final AllAtomWrapper allAtomAction = new AllAtomWrapper();
    private final class AllAtomWrapper extends AtomAction {
        private AtomAction action;
        void setAction(AtomAction a) {this.action = a;}
        public void actionPerformed(Atom a) {((SpeciesAgent)a).allAtoms(action);}
    }
}