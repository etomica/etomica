package etomica;
import etomica.units.Unit;
import java.util.Random;
import java.beans.Beans;

//Java2 imports
//import java.util.HashMap;
//import java.util.Iterator;

import etomica.utility.HashMap;
import etomica.utility.Iterator;
import etomica.utility.LinkedList;


     /**
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
      * @see Molecule
      * @see Simulation#getPotential
      * @see Species.Agent
      * @see Species.Reservoir
      */
     
public abstract class Species implements Simulation.Element, java.io.Serializable {

    public static final String VERSION = "Species:01.07.09";
    private Simulation parentSimulation;
    private boolean added = false;
    protected AtomFactory factory;
    private int index;
    
    public Species(Simulation sim, AtomFactory factory) {
        parentSimulation = sim;
        this.factory = factory;
        parentSimulation.register(this);
        index = ((LinkedList)sim.speciesList()).size();
        name = "Species "+Integer.toString(index);
    }
          
    public final Simulation parentSimulation() {return parentSimulation;}
    public final Class baseClass() {return Species.class;}
    public final boolean wasAdded() {return added;}
    public final void setAdded(boolean b) {added = b;}
    
    public AtomFactory moleculeFactory() {return factory;}
    
    public AtomIterator makeAtomIterator(Phase p) {
        return getAgent(p).new LeafAtomIterator();
    }
    public AtomIterator makeMoleculeIterator(Phase p) {
        return getAgent(p).new ChildAtomIterator();
    }
        
    /**
     * Nominal number of molecules of this species in each phase.
     * Actual number may differ if molecules have been added or removed to/from the phase
     */
    protected int nMolecules;
    
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
        Iterator e = agents.values().iterator();
        while(e.hasNext()) {
            SpeciesAgent a = (SpeciesAgent)e.next();
            a.setNMolecules(n);
        }
    }
    
    /**
     * Accessor method for nominal number of molecules in each phase.  Actual 
     * number of molecules of this species in a given phase is obtained via
     * the getNMolecules method of Species.Agent
     * 
     * @return Nominal number of molecules in each phase
     * @see Species.Agent#getNMolecules
     */
    public int getNMolecules() {return nMolecules;}
    
    /**
     * Accessor method of the name of this species
     * 
     * @return The given name of this species
     */
    public String getName() {return name;}
    
    public int index() {return index;}

    /**
     * Method to set the name of this species
     * The species' name provides a convenient way to label output data that is 
     * associated with this species.  This method might be used, for example, to place
     * a heading on a column of data.
     * Default name is "Species" followed by the integer species index of this species.
     * 
     * @param name The name string to be associated with this species
     */
    public final void setName(String name) {this.name = name;}

    /**
     * Overrides the Object class toString method to have it return the output of getName
     * 
     * @return The name given to the species
     */
    public String toString() {return getName();}  //override Object method
          
    /**
     * Constructs an Agent of this species and sets its parent phase
     * 
     * @param p The given parent phase of the agent
     * @return The new agent.  Normally this is returned into the setAgent method of the given phase.
     * @see Species.Agent
     */
    public SpeciesAgent makeAgent(SpeciesMaster parent) {
        Phase phase = parent.parentPhase();
        SpeciesAgent a = new SpeciesAgent(this, parent, nMolecules);
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
    public final SpeciesAgent getAgent(Phase p) {return (SpeciesAgent)agents.get(p);}
        
    /*
    * A name to be associated with the species.  Use is optional.
    */
    String name;
          
              
    /**
     * Hashtable to associate agents with phases
     */
    final HashMap agents = new HashMap();

}