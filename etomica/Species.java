package etomica;
import etomica.units.Unit;
import java.util.Random;
import java.beans.Beans;

//Java2 imports
//import java.util.HashMap;
//import java.util.Iterator;

import etomica.utility.HashMap;
import etomica.utility.Iterator;


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
    private boolean stationary;
    private Simulation parentSimulation;
    private boolean added = false;
    private AtomFactory factory;
    
    public Species(Simulation sim, AtomFactory factory) {
        parentSimulation = sim;
        this.factory = factory;
        parentSimulation.register(this);
    }
          
    public final Simulation parentSimulation() {return parentSimulation;}
    public final Class baseClass() {return Species.class;}
    public final boolean wasAdded() {return added;}
    public final void setAdded(boolean b) {added = b;}
          
        
           /**
            * Nominal number of molecules of this species in each phase.
            * Actual number may differ if molecules have been added or removed to/from the phase
            */
           //Temporary method to handle setting of nMolecules for all phases
           //Allows only one value for all phases
    protected int nMolecules;
    /**
     * Sets the number of molecules of this species for each phase
     * Propagates the change to agents of this species in all phases, so
     * it creates the given number of molecules in every phase.  Any existing molecules
     * of the species are discarded.
     * 
     * @param n The new number of molecules of this species in each phase
     * @see Species.Agent#setNMolecules
     */
    public void setN(int n) {
        nMolecules = n;
        Iterator e = agents.values().iterator();
        while(e.hasNext()) {
            Agent a = (Agent)e.next();
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
     * Creates new molecule of this species without associating it with a phase.
     */
    private Molecule makeMolecule() {return makeMolecule(null);}

    /**
     * Concrete instances of this abstract class invoke the Molecule constructor
     * with the atom type and molecule configuration appropriate to make a molecule
     * of the type defined by this species.  This method is invoked by the 
     * getMolecule method, which provides the means for a user to obtain a molecule
     * of this species.
     * 
     * @param p Phase in which the returned molecule is placed
     * @return A molecule of this species
     * @see #getMolecule
     */
    protected abstract Molecule makeMolecule(Phase p);

    /**
     * Specifies the number of atoms in a molecule of this species.
     * Since the number of atoms in a molecule cannot be changed once the molecule 
     * is constructed, to have the atom/molecule change take effect it is necessary to 
     * create new molecules of this species in each phase.  Thus the method invokes
     * the setNMolecules method of each agent of this species, replacing all the existing
     * molecules with ones having the newly prescribed value of atomsPerMolecule.
     * 
     * @param na The new number of atoms per molecule
     */
    public void setAtomsPerMolecule(int na) {
        if(na == atomsPerMolecule) return;  //do nothing is value isn't changings
        atomsPerMolecule = na;
//        if(parentSimulation == null) {return;}
        Iterator e = agents.values().iterator();
        while(e.hasNext()) {
            Agent a = (Agent)e.next();
            a.setNMolecules(a.nMolecules, true);//2nd argument indicates to make new molecules even though number of them is not changing
        }
    }
            
    /**
     * Accessor method for number of atoms per molecule
     * 
     * @return Present number of atoms in a molecule of this species.
     */
    public int getAtomsPerMolecule() {return atomsPerMolecule;}
            
            
    /**
     * Accessor method of the name of this species
     * 
     * @return The given name of this species
     */
    public String getName() {return name;}

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
     * Sets the value of the index of this species
     * Value is set by the register method of Simulation when the species is constructed
     * Index is used to associate species with the potentials governing interactions
     * of molecules of this species with each other and with the molecules of other species.
     * 
     * @param index The new index of this species
     * @see Simulation.register(Species)
     */
    public final void setSpeciesIndex(int index) {
        speciesIndex = index;
        if(name == null) name = "Species "+Integer.toString(index);
    }
                    
    /**
     * Constructs an Agent of this species and sets its parent phase
     * 
     * @param p The given parent phase of the agent
     * @return The new agent.  Normally this is returned into the setAgent method of the given phase.
     * @see Species.Agent
     */
    public SpeciesAgent makeAgent(SpeciesMaster parent, int index) {
        Agent a = new Agent(this, parent, index, nAtoms, 
                            parent.parentPhase().coordinateInitializer());
        agents.put(p,a);   //associate agent with phase; retrieve agent for a given phase using agents.get(p)
        return a;
    }

    /**
     * Returns the agent of this species in the given phase
     * Hashmap is used to connect phase(key)-agent(value) pairs
     * 
     * @param p The phase for which this species' agent is requested
     * @return The agent of this species in the phase
     */
    public final Agent getAgent(Phase p) {return (Agent)agents.get(p);}
    
    /**
     * Primary method for other classes to obtain a molecule of this species
     * Molecule is obtained from the species reservoir, if it is not empty, otherwise
     * a new molecule is constructed and returned.
     * The container of the returned molecule is this species' reservoir.  It is
     * removed from the reservoir when it is added to another molecule container 
     * (typically a phase).
     * 
     * @return A molecule of this species
     * @see Species.Reservoir
     */
    public final Molecule getMolecule() {
        if(reservoir.isEmpty()) {
            Molecule m = makeMolecule();
            reservoir.addNewMolecule(m);
        }
        return reservoir.nextMolecule();
    }

    /**
     * @return The molecule reservoir for this species
     * @see Species.Reservoir
     */
    public final Reservoir reservoir() {return reservoir;}
    
    /*
    * A name to be associated with the species.  Use is optional.
    */
    String name;
          
              
    /**
     * Hashtable to associate agents with phases
     */
    final HashMap agents = new HashMap();

    /**
     * Class for storing molecules that are not presently in any phase
     */
    private final Reservoir reservoir = new Reservoir(this);
          
    private final Random rand = new Random();  //for choosing molecule at random

    
    /**
     * Class to store molecules that are not presently in any particular phase
     * Useful for semigrand and grand canonical simulations
     * Keeps molecules in a linked list
     */
    protected static final class Reservoir implements Molecule.Container, java.io.Serializable {
        
        private Molecule next;
        private Species species;
        private boolean added = false;
        
        public Reservoir(Species s) {
            species = s;
        }
        
/*        //Simulation.Element interface extended by Molecule.Container
        public Simulation parentSimulation() {return species.parentSimulation();}
        public Class baseClass() {return Species.Reservoir.class;}
        public final boolean wasAdded() {return added;}
        public final void setAdded(boolean b) {added = b;}
*/
        public Molecule nextMolecule() {return next;}  //peek

      /**
       * Adds a freshly made molecule
       * Same as addMolecule, but does not attempt to remove from original container first
       */
        private void addNewMolecule(Molecule m) {
            m.setNextMolecule(next);
            m.setContainer(this);
            next = m;
        }
        
        //does not handle m==null
        //does not clear parentPhase
        //Molecule.Container method
        public void addMolecule(Molecule m) {
            m.container.removeMolecule(m);
            m.setNextMolecule(next);
            m.setContainer(this);
            next = m;
        }
        
        //Molecule.Container method
        public void removeMolecule(Molecule m) {
            if(m != next) {  //should throw an exception
                System.out.println("Error:  removal of inappropriate molecule from reservoir");
            }
            if(next != null) next = next.nextMolecule();
        }
        
        public boolean isEmpty() {return next==null;}
    }
}