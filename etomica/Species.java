package etomica;
import etomica.units.Unit;
import java.util.Random;
import java.util.HashMap;
import java.util.Iterator;
import java.beans.Beans;

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

    public static final String VERSION = "Species:01.02.14.0";
    private boolean stationary;
    private Simulation parentSimulation;
    private boolean added = false;
    
//    public static int count = 0;
    
    public Species(Simulation sim) {
        parentSimulation = sim;
//        setSpeciesIndex(count++);
        parentSimulation.register(this);
    }
          
    public final Simulation parentSimulation() {return parentSimulation;}
    public final Class baseClass() {return Species.class;}
    public final boolean wasAdded() {return added;}
    public final void setAdded(boolean b) {added = b;}
          
    /**
     * Changes the default molecule configuration associated with this species.  
     * Propagates the change to agents of this species in all phases.
     * 
     * @param mc The new Molecule.Configuration object
     * @see Molecule.Configuration
     */
    public void add(Molecule.Configuration mc) {
        mc.setParentSpecies(this);
        this.moleculeConfiguration = mc;
        Iterator e = agents.keySet().iterator();
        while(e.hasNext()) {
            Phase p = (Phase)e.next();
            mc.initializeCoordinates(p);
        }
    }
           
    /**
     * Sets the "stationary" field of all atoms of this species to the given value
     */
    public void setStationary(boolean b) {
        stationary = b;
        Iterator e = agents.values().iterator();
        while(e.hasNext()) {
            Agent a = (Agent)e.next();
            Atom.Iterator aIter = a.makeAtomIterator();
            while(aIter.hasNext()) {aIter.next().setStationary(b);}
        }
    }
    
    public boolean getStationary() {return stationary;}
        
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
    public void setNMolecules(int n) {
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
     * Sets the next species in the linked list of species
     * and sets this species to be the previous species of the given species
     * 
     * @param s The species to follow this one in the linked list of species
     */
    public final void setNextSpecies(Species s) {
        this.nextSpecies = s;
        if(s != null) s.previousSpecies = this;
    }

    /**
     * Accessor method of the species following this one in the linked list of species.
     * 
     * @return The next species in the linked list
     */
    public final Species nextSpecies() {return nextSpecies;}

    /**
     * Accessor method of the species preceding this one in the linked list of species
     * 
     * @return The previous species in the linked list of species
     */
    public final Species previousSpecies() {return previousSpecies;}
            
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
     * Accessor method for the index of this species
     * 
     * @param index The index of this species
     */
    public final int getSpeciesIndex() {return speciesIndex;}

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
    public Agent makeAgent(Phase p) {
        Agent a = new Agent(p);
        a.setNMolecules(nMolecules);
        Atom.Iterator aIter = a.makeAtomIterator();
        while(aIter.hasNext()) {aIter.next().setStationary(stationary);}
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
     * Returns an atom iterator for this species in the given phase
     */
     public final Atom.Iterator makeAtomIterator(Phase p) {
        if(p != null) return getAgent(p).makeAtomIterator();
        else return null;
     }
    
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
          
    public ColorScheme colorScheme;
         
    /**
    * A unique integer that identifies the species.  Used to associate
    * inter- and intra-molecular potentials to the Species.  Assigned
    * by Simulation when species is registered.  Indices are assigned starting
    * from zero and increase in the order that the species are constructed.
    * 
    * @see Phase#potential2
    * @see Phase#potential1
    */
    int speciesIndex;
          
    /**
    * Number of atoms in each molecule of this species
    */
    protected int atomsPerMolecule;
              
    /**
     * Hashtable to associate agents with phases
     */
    HashMap agents = new HashMap();

    /**
    * Object responsible for setting default configuration of atoms in molecule
    */
    public Molecule.Configuration moleculeConfiguration;
    
    Species previousSpecies, nextSpecies;
    
    /**
     * Class for storing molecules that are not presently in any phase
     */
    private final Reservoir reservoir = new Reservoir(this);
          
    private final Random rand = new Random();  //for choosing molecule at random

        /**
         * The Species.Agent is a representative of the species in each phase.
         * The agent handles addition, deletion, link-list ordering, counting, etc. of 
         * molecules in a phase.  Each phase has an agent from every species instance.
         * 
         * @author David Kofke
         * @see Species
         * @see Phase
         */
    public final class Agent implements java.io.Serializable {
        public Phase parentPhase;  //want final, but won't compile
        private Atom.Iterator atomIterator;  //iterator for all atoms of this species in this phase
        
        /**
         * @param p The phase in which this agent will be placed
         */
        public Agent(Phase p) {
            parentPhase = p;
            atomIterator = new AtomIterator();
        }
        
        public final Species parentSpecies() {return Species.this;}
        public final Phase parentPhase() {return parentPhase;}
        
        /**
         * Returns a new, unique instance of an atom iterator for the species in the phase
         */
        public final Atom.Iterator makeAtomIterator() {return new AtomIterator();}
        
        /**
        * @return the number of molecules for this species in this phase
        */
        public final int moleculeCount() {return nMolecules;}
              
        /**
        * Sets the number of molecules for this species.  Makes the given number
        * of new molecules, linked-list orders and initializes them.
        * Any previously existing molecules for this species in this phase are abandoned
        * Any links to molecules of next or previous species are maintained.
        * Takes no action at all if the new number of molecules equals the existing number
        *
        * @param n  the new number of molecules for this species
        * @see #makeMolecule
        * @see #deleteMolecule
        * @see #addMolecule
        */
        public void setNMolecules(int n) {
            setNMolecules(n, false);
        }
        
        /**
         * Same as setNMolecules, but takes a boolean argument that can indicate that all
         * new molecules should be made even if the number of molecules is not changing.
         * @param forceRebuild if true will cause making of new molecules even if number is not changing.
         */
        public void setNMolecules(int n, boolean forceRebuild) {
            if(n == nMolecules && !forceRebuild) return; //do nothing if number is not changing
            Integrator integrator = parentPhase.integrator();
            boolean wasPaused = false;
            if(integrator != null) {
                wasPaused = integrator.isPaused();//record pause state of integrator
                integrator.pause();
            }
            Molecule previous = null;
            Molecule next = null;
            if(firstMolecule != null) previous = firstMolecule.previousMolecule();
            if(lastMolecule != null) next = lastMolecule.nextMolecule();
            nMolecules = n;
            if(nMolecules == 0) {
                firstMolecule = null;
                lastMolecule = null;
                if(previous != null) previous.setNextMolecule(next);
            }
            else {
                firstMolecule = makeMolecule(parentPhase);
                lastMolecule = firstMolecule;
                for(int i=1; i<nMolecules; i++) {
                    lastMolecule.setNextMolecule(makeMolecule(parentPhase));
                    lastMolecule = lastMolecule.nextMolecule();
                }
                lastMolecule.setNextMolecule(next);
                if(previous != null) previous.setNextMolecule(firstMolecule);
            }
            parentPhase.updateCounts();
            parentPhase.configuration.initializeCoordinates(parentPhase);
            parentPhase.iteratorFactory().reset();
            if(integrator != null) {
                if(integrator.isInitialized()) integrator.initialize();//reinitialize only if initialized already
                if(!wasPaused) integrator.unPause();//resume if was not paused originally
            }
        }
              
        /**
        * Chooses a molecule randomly from Species
        *
        * @return the randomly seleted molecule
        */
        public final Molecule randomMolecule() {
            int i = (int)(rand.nextDouble()*nMolecules);
            Molecule m = firstMolecule;
            for(int j=i; --j>=0; ) {m = m.nextMolecule();}
            return m;
        }
                
        /**
        * Chooses a molecule randomly from Species, and deletes it (removes it from the linked list)
        *
        * @return the deleted molecule
        */
 /*       public final Molecule deleteMolecule() {
            Molecule m = this.randomMolecule();
            deleteMolecule(m);
            return m;
        }
*/                
        /**
        * Removes molecule from species, and updates atom and molecule linked lists.
        * Updates all values of first/last Molecule/Atom for species and
        * phase, if appropriate.  Also updates number-of-atom/molecule variables
        * for species and phase.
        * No measures are taken to remove this species if it holds zero molecules
        * after this molecule is deleted.
        *
        * @param m the molecule being deleted
        * @see #addMolecule
        */
              
        public final void deleteMolecule(Molecule m) {
            if(m.parentSpecies != Species.this) {
                System.out.println("Error:  attempt to delete molecule from incorrect species");
                return;
            }
            Molecule next = m.nextMolecule();
            Molecule previous = m.previousMolecule();
            if(m == firstMolecule) {
                if(nMolecules == 1) {firstMolecule = null;}  //deleting the first and only molecule of the species
                else {firstMolecule = next;}                 //deleting first molecule, but others are present
            }
            if(m == lastMolecule) {
                if(nMolecules == 1) {lastMolecule = null;}
                else {lastMolecule = previous;}
            }
            if(previous != null) {previous.setNextMolecule(next);} //reconnect linked list if not at beginning
            else if(next != null) {next.clearPreviousMolecule();}  //beginning of list; no previous molecule for next
            nMolecules--;
       //     m.parentPhase = null;        //line deleted because of spareMolecule
            m.setNextMolecule(null);
            m.clearPreviousMolecule();
        //    parentPhase.moleculeCount--;
        //    parentPhase.atomCount -= m.nAtoms;
        //    if(spareMolecule == null) spareMolecule = m;
        }

        /**
        * Adds a molecule to this species and updates linked lists.  Does not handle
        * creation of molecule.  New molecule
        * becomes last molecule of species.  Updates first/last Molecule
        * for species, if appropriate.
        * Not yet correctly implemented for use in a phase containing multiple
        * species (i.e., mixtures).  Does not adjust total molecule and atom count in parent phase.
        *
        * @param m the molecule being added
        * @see deleteMolecule
        */
        public final void addMolecule(Molecule m) {
            if(nMolecules > 0) {
                m.setNextMolecule(lastMolecule.nextMolecule());
                lastMolecule.setNextMolecule(m);
                lastMolecule = m;
            }
            else {  //m is going to be the only molecule in species
                firstMolecule = m;
                lastMolecule = m;
                m.setNextMolecule(null); 
                for(Species.Agent s=this.nextSpecies(); s!=null; s=s.nextSpecies()) { //loop forward in species, looking for next molecule
                    if(s.firstMolecule() != null) {
                        m.setNextMolecule(s.firstMolecule());
                        break;
                    }
                }
                m.clearPreviousMolecule();
                for(Species.Agent s=this.previousSpecies(); s!=null; s=s.previousSpecies()) { //loop backward in species, looking for previous molecule
                    if(s.lastMolecule() != null) {
                        s.lastMolecule.setNextMolecule(m);
                        break;
                    }
                }
            }
            nMolecules++;
//            m.parentPhase = parentPhase;
//            parentPhase.moleculeCount++;
//            parentPhase.atomCount += m.nAtoms;
//            colorScheme.initializeMoleculeColor(m);  //need to deal with this for changes in colorscheme
        }
              
        /**
        * @return the next species in the linked list of species.  Returns null if this is the last species.
        */
        public final Agent nextSpecies() {return nextSpecies;}
             
        /**
        * Sets the species following this one in the linked list of species.
        * Also links last Molecule/Atom of this species to the corresponding
        * first Molecule/Atom of the next species
        *
        * @param s the species to be designated as this species nextSpecies
        */
        public final void setNextSpecies(Agent s) {
            this.nextSpecies = s;
            Molecule last = lastMolecule();
            if(s==null) {
                if(last!=null) last.setNextMolecule(null); 
                return;
            }
            s.previousSpecies = this;
            if(last != null) {last.setNextMolecule(s.firstMolecule);}
        }
        /**
        * @return the species preceding this one in the linked list of species.  Returns null if this is the first species.
        */
        public final Agent previousSpecies() {return previousSpecies;}   
              
        /**
        * Method that sets the initial coordinates of the molecules in the
        * species.  If fillVolume is <code>true</code>, will set coordinate
        * origin and value to fill entire space in Phase; if fillVolume is <code>
        * false</code>, will not do this, and thereby permit adjustment of 
        * origin and size of occupied volume at design time.
        * Called by paint at design time, and by Phase.add at run time, and by this.setBounds (override of superclass)
        */
        //  public abstract void initializeSpecies(Phase phase);
              
        public void setBounds(int x, int y, int width, int height) {
        ///    Rectangle r = getBounds();
        ///    if(r.x!=x || r.y!=y || r.width!=width || r.height!=height) {  
        ///        super.setBounds(x, y, width, height);
        //        if(parentSimulation != null) initializeSpecies(parentSimulation);
        ///    }
        }

        public final Molecule firstMolecule() {return firstMolecule;}
        public final Molecule lastMolecule() {return lastMolecule;}
              
        /**
        * Used to terminate loops over molecules in species
        */
        public final Molecule terminationMolecule() {
            return (lastMolecule == null) ? null : lastMolecule.nextMolecule();
        }  
        public final Atom firstAtom() { //return firstAtom;
            return (firstMolecule == null) ? null : firstMolecule.firstAtom();
        }
        public final Atom lastAtom() { //return lastAtom;
            return (lastMolecule == null) ? null : lastMolecule.lastAtom();
        }
              
        /**
        * Used to terminate loops over atoms in species
        */
        public final Atom terminationAtom() {
            Atom last = lastAtom();
            return (last == null) ? null : last.nextAtom();
        }
              
        /**
        * Number of molecules for this species
        */
        protected int nMolecules;
              
        /**
        * The Species following this one in the linked list of Species
        */
        Agent nextSpecies;
             
        /**
        * The Species preceding this one in the linked list of Species
        */
        Agent previousSpecies;
              
        /**
        * First molecule in this Species
        *
        * @see Molecule
        */
        protected Molecule firstMolecule;
              
        /**
        * Last molecule in this Species
        *
        * @see Molecule
        */
        protected Molecule lastMolecule;
        
        /**
         * Iterator for all atoms of this species in this phase
         */
        public final class AtomIterator implements Atom.Iterator {
            private Atom atom, nextAtom;
            private boolean hasNext;
            public AtomIterator() {reset();}
            public boolean hasNext() {return hasNext;}
            public void reset() {
                atom = firstAtom();
                hasNext = (atom != null);
            }
            public void reset(Atom a) {reset();}
            public Atom next() {
                nextAtom = atom;
                if(atom == lastAtom()) {hasNext = false;}
                else {atom = atom.nextAtom();}
                return nextAtom;
            }
            public void allAtoms(AtomAction act) {
                Atom term = terminationAtom();
                for(Atom a=firstAtom(); a!=term; a=a.nextAtom()) {act.actionPerformed(a);}
            }
        } //end of AtomIterator

    } //end of Agent
    
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