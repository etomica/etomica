package etomica;

import java.util.Observable;
import java.util.Observer;

//Java2 imports
//import java.util.LinkedList;
//import java.util.Iterator;

import etomica.utility.LinkedList;
import etomica.utility.Iterator;


/**
 * This description is out of date.
 *
 * A Phase collects all atoms that interact with one another; atoms in different phases do
 * not interact.  These are the important features of a Phase:<p> 
 * 
 * 1. It maintains an IteratorFactory that can be used to obtain iterators that loop over all
 * or some atoms and atom pairs in the phase.<p>
 * 
 * 2. It provides methods for addition and removal of molecules in the phase.<p>
 * 
 * 3. It holds a Boundary object, obtained from the governing Space, that defines the behavior
 * of atoms as they move into or across the boundary of the phase.<p>
 * 
 * 4. It holds a set of Meter objects that can be used to measure properties of that atoms
 * held by the phase.  An EnergyMeter is constructed by default to permit evaluation of
 * potential and kinetic energies of the atoms in the phase.<p>
 * 
 * 5. It has a Configuration object that determines the default initial configuration of the
 * atoms in the phase.<p>
 *
 * 6. It holds a potential agent that manages all potentials governing the interactions
 * among the atoms in the phase.<p>
 *
 * 
 * @author David Kofke
 * @see IteratorFactory
 * @see Space.Boundary
 * @see MeterAbstract
 */
public final class Phase extends SimulationElement {
        
    private Space.Boundary boundary;
    private transient final LinkedList meterList = new LinkedList();
    private PhaseAction.Inflate inflater;
    public final SpeciesMaster speciesMaster;
    private boolean lrcEnabled = true;
    
    public Phase() {
        this(Simulation.instance);
    }
    
    public Phase(Simulation sim) {
        super(sim, Phase.class);
        
        speciesMaster = new SpeciesMaster(this);
        
        moleculeIterator = makeMoleculeIterator();//must come after speciesMaster assignment
        
        inflater = new PhaseAction.Inflate(this);
        //don't use setIteratorFactory here to remove any possibility of complications with Observers
        if(sim.space() instanceof IteratorFactory.Maker) {
            iteratorFactory = ((IteratorFactory.Maker)sim.space()).makeIteratorFactory(this);
 //           if(iteratorFactory instanceof PotentialField.Maker) {
 //               addField(((PotentialField.Maker)iteratorFactory).makePotentialField(this));
 //           }
        }
        else {
            iteratorFactory = new IteratorFactory(this);
        }

//        setBoundary(Space.Boundary.DEFAULT);
        setBoundary(parentSimulation().space().makeBoundary());

        if(parentSimulation().space().D() < 3) 
            setConfiguration(new ConfigurationSequential(parentSimulation().space()));  //default configuration
        else
            setConfiguration(new ConfigurationFcc(parentSimulation().space()));
    }//end of constructor
   
    /**
     * Mutator method for flag that enables or disables application of long-range
     * correction to truncated potentials.  Enabled by default.
     */
    public void setLrcEnabled(boolean b) {lrcEnabled = b;}
    /**
     * Accessor method for flag that enables or disables application of long-range
     * correction to truncated potentials.  Enabled by default.
     */
    public boolean isLrcEnabled() {return lrcEnabled;}
    
    public void addSpecies(Species s) {
        speciesMaster.addSpecies(s);
//        moleculeIterator.setBasis(speciesMaster); this is not needed because molecule iterator is a listener to species master
    }
    
    public Space.Vector randomPosition() {return boundary.randomPosition();}
        
    /**
     * Returns an array of Space.Vector with the center-of-mass position of each molecule
     * Not meant for computation-intensive use.
     */
/*    public Space.Coordinate[] getMoleculePosition() {
        Space.Coordinate[] positions = new Space.Coordinate[moleculeCount];
        Molecule molecule = firstMolecule();
        for(int i=0; i<moleculeCount; i++) {
            positions[i] = molecule.coordinate();
            molecule = molecule.nextMolecule();
        }
        return positions;
    }
*/    
    /**
     * Returns an array of the molecules currently in the phase, for those times
     * when you just need convenient indexed access to each molecule.
     * Use is discouraged.
     */
/*     public Molecule[] moleculeArray() {
        Molecule[] array = new Molecule[moleculeCount];
        int i=0;
        Molecule terminationMolecule = lastMolecule().nextMolecule();
        for(Molecule m=firstMolecule(); m!=terminationMolecule; m=m.nextMolecule()) {
            array[i++] = m;
        }
        return array;
     }
 */   
    /**
     * Returns the ith molecule in the linked list of molecules.
     * 0 returns the first molecule, and moleculeCount-1 returns the last.
     * An argument outside this range throws an IndexOutOfBoundsException
     */
    //could make more efficient by starting from first or last molecule, as appropriate
    public Atom molecule(int i) {
        if(i >= moleculeCount() || i < 0) 
            throw new IndexOutOfBoundsException("Index: "+i+
                                                ", Number of molecules: "+moleculeCount());
        int sum = 0;
        SpeciesAgent s;
        for(s=speciesMaster.firstSpecies(); s!=null; s=s.nextSpecies()) {
            sum += s.node.childAtomCount();
            if(sum > i) break;
        }
        return s.node.getAtom(i-(sum-s.node.childAtomCount()));
    }
    
    /**
     * Returns a randomly selected molecule from the phase.
     */
    public Atom randomMolecule() {
        int i = (int)(moleculeCount() * Simulation.random.nextDouble());
        return molecule(i);
    }
      
      public SpeciesMaster speciesMaster() {return speciesMaster;}
     
    /**
     * Sets the boundary object of the phase.
     */
     public void setBoundary(Space.Boundary b) {
        boolean hasIntegrator = integrator != null;
        if(hasIntegrator) integrator.pause();
        boundaryMonitor.notifyObservers(b);
        boundary = b;
        boundary.setPhase(this);
        if(hasIntegrator) {
            integrator.reset();
            integrator.unPause();
        }
     }

    /**
     * Indicates the Boundary object of the phase
     * 
     * @return An integer that codes for the boundary via a static variable in Space.Boundary
     */
//    public final int getBoundary() {return iBoundary;}
    public final Space.Boundary getBoundary() {return boundary;}
    public final Space.Boundary boundary() {return boundary;}
    
    /**
     * Accessor of the phase's iterator factory object
     * 
     * @return The phase's iterator factory
     * @see IteratorFactory
     */
    public final IteratorFactory iteratorFactory() {return iteratorFactory;}

    /**
     * Sets the iterator factory of the phase
     * 
     * @param it
     * @see IteratorFactory
     */
    public final void setIteratorFactory(IteratorFactory it) {
        //notify observers before updating iteratorFactory in case observers want to close business with old iteratorFactory
        iteratorFactoryMonitor.notifyObservers(it);
        iteratorFactory = it;
    }
    
    /**
     * Accessor of the integrator governing the movement of the atoms in the phase.
     * 
     * @return The phase's integrator.
     * @see Integrator
     */
    public final Integrator integrator() {return integrator;}

    /**
     * Sets the integrator of the phase.
     * 
     * @param i The new integrator.
     */
    public void setIntegrator(Integrator i) {
           //need to define an exception for this
        if(!i.addPhase(this)) return; //addPhase will return false if the phase was not successfully added to the integrator
        //notify observers before updating integrator in case observers want to close business with old integrator
        integratorMonitor.notifyObservers(i);
//        integratorMonitor.fireEvent(new PhaseIntegratorEvent(i,PhaseIntegratorEvent.NEW_INTEGRATOR));
        if(integrator != null) integrator.removePhase(this);
        integrator = i;
    }
                        
    public final Space.Vector dimensions() {return boundary.dimensions();}
    public final double volume() {return boundary.volume();}  //infinite volume unless using PBC
    
    public void setDensity(double rho) {
        double vNew = moleculeCount()/rho;
        double scale = Math.pow(vNew/boundary.volume(), 1.0/parentSimulation().space().D());
        inflater.actionPerformed(scale);
    }
    public double getDensity() {return moleculeCount()/boundary.volume();}

    /**
     * @return the first atom in the linked list of atoms in this Phase
     */
    public final Atom firstAtom() {
        return speciesMaster.node.firstLeafAtom();
    }
    
    /**
     * @return the last atom in the linked list of atoms in this Phase
     */
    public final Atom lastAtom() {
        return speciesMaster.node.lastLeafAtom();
    }
    
    /**
     * @return the first Agent in the linked list of Species.Agents in this phase
     */
    public final SpeciesAgent firstSpecies() {return (SpeciesAgent)speciesMaster.node.firstChildAtom();}
    
    /**
     * @return the last Agent in the linked list of Species.Agents in this phase
     */
    public final SpeciesAgent lastSpecies() {return (SpeciesAgent)speciesMaster.node.lastChildAtom();}
    
    public int moleculeCount() {return speciesMaster.moleculeCount();}
    
    public int atomCount() {return speciesMaster.node.leafAtomCount();}
    
    /**
    * Returns the temperature (in simulation units) of this phase as computed via the equipartition
    * theorem from the kinetic energy summed over all (atomic) degrees of freedom
    */  
    public double kineticTemperature() {
        return MeterTemperature.currentValue(this);
    }
    
    public void setConfiguration(Configuration c) {
        configuration = c;
        configuration.initializeCoordinates(speciesMaster.node.childAtomArray());
        iteratorFactory.reset();
    }
    
    public Configuration getConfiguration() {return configuration;}

	/**
	 * Adds a meter to the list of meters working in this phase.
	 */
	public void addMeter(MeterAbstract m) {meterList.add(m);}

	/**
	 * Removes a meter from the list of meters working in this phase.
	 */
	public void removeMeter(MeterAbstract m) {meterList.remove(m);}

	/**
	 * Returns the list of meters working in this phase.
	 */
	public LinkedList getMeterList() {return meterList;}
    
    /**
     * Deploys the agent of a species in this phase
     */
//    void addSpecies(SpeciesAgent species) {
 //       speciesMaster.addAtom(species);
        //set internal configuration of molecule
   //     if(species.parentSpecies().moleculeConfiguration != null) species.parentSpecies().moleculeConfiguration.initializeCoordinates(this);
        //add species to configuration for this phase and notify iteratorFactory
   //     configuration.initializeCoordinates(speciesMaster.childAtomArray());
   //     iteratorFactory.reset();  
 //   }
    
    /**
     * Adds the given molecule to this phase, placing it in the molecule/atom linked lists
     * and removing it from the container it previously resided.
     * Looks up the molecule's species agent in this phase and passes it with the 
     * molecule to the two-argument addMolecule method.  If the agent is known already, 
     * that form of the method should be used instead.
     */
    public void addMolecule(Atom a, Species s) {
        addMolecule(a, s.getAgent(this));
    }
    
    /**
     * Adds the given molecule to this phase.
     * @param molecule the molecule to be added
     * @param s the species agent in this phase for the molecule's species.  If known, this agent can be provided to save the effort of looking it up.
     */
    public void addMolecule(Atom a, SpeciesAgent s) {
        if(a == null || s == null) return;
        s.node.addAtom(a);
    }
    
    /**
     * Removes molecule from phase.  
     */
    public void removeMolecule(Atom a) {
        if(a == null) return;
        removeMolecule(a, a.node.parentSpeciesAgent());
    }
    public void removeMolecule(Atom a, SpeciesAgent s) {
        if(a == null || s == null) return;
        s.node.removeAtom(a);
    }

//need a better way
    /**
    * Synchronized version of deleteMolecule.  
    * Useful if molecules are being deleted by GUI events, rather than by integrator 
    */
//    public final synchronized void deleteMoleculeSafely(Molecule m) {  //will this make deleteMolecule synchronized?
 //       deleteMolecule(m);
//    }
    
    /**
    * Synchronized version of addMolecule
    * Useful if molecules are being added by GUI events, rather than by integrator 
    */
//    public final synchronized void addMoleculeSafely(Molecule m) {
//        addMolecule(m);
//    }
                                 
    public Configuration configuration;
          
    private Integrator integrator;
    
    private IteratorFactory iteratorFactory;
          
//    public transient MeterEnergy energy;
        
    public Phase.Monitor integratorMonitor = new Phase.Monitor();
    public Phase.Monitor boundaryMonitor = new Phase.Monitor();
    public Phase.Monitor iteratorFactoryMonitor = new Phase.Monitor();
//    public SimulationEventManager integratorMonitor = new SimulationEventManager();
    
    public static class Monitor extends Observable implements java.io.Serializable {
        
        public void notifyObservers() {
            this.notifyObservers(null);
        }
        public void notifyObservers(Object obj) {
            setChanged();
            super.notifyObservers(obj);
        }
        public void addObserver(Observer o) {if(o != null) super.addObserver(o);}
    }
     
    /**
     * Class for constructing linked lists of Phases.
     * Each Linker points to one Phase and another Linker, the next one in the list.
     */
    public static class Linker implements java.io.Serializable {
        private final Phase phase;
        private Phase.Linker next = null;
        //Constructors
        public Linker(Phase a) {phase = a;}
        public Linker(Phase a, Linker l) {phase = a; next = l;}
        //Access methods
        public final Phase phase() {return phase;}
        public final Phase.Linker next() {return next;}
        public final void setNext(Phase.Linker l) {next = l;}
    }//end of Phase.Linker
    
    /**
     * Makes an iterator that loops through all the (leaf) atoms present in this phase.
     */
    public AtomIterator makeAtomIterator() {
        return new AtomIteratorSequential(speciesMaster, true); //instead make an enumerated type to key for leaf iterator
    //        public boolean contains(Atom a) {return a != null && a.parentPhase() == Phase.this;}
    }
    
    /**
     * Makes an iterator that loops through all the molecules 
     * (children of the species agents) in this phase.
     */
    public AtomIterator makeMoleculeIterator() {
        return new AtomIteratorCompound(speciesMaster);
    }//end of makeMoleculeIterator
    
    /**
     * Makes an iterator that loops through all the (leaf) atoms 
     * derived from the given species in this phase.
     */
    public AtomIterator makeAtomIterator(Species s) {
        return new AtomIteratorSequential((SpeciesAgent)s.getAgent(this), true); //instead make an enumerated type to key for leaf iterator
//        return ((SpeciesAgent)s.getAgent(this)).new LeafAtomIterator();
    }

    /**
     * Makes an iterator that loops through all the child atoms (molecules) 
     * of the given species in this phase.
     */
    public AtomIterator makeMoleculeIterator(Species s) {
        return new AtomIteratorSequential((SpeciesAgent)s.getAgent(this));
 //       return ((SpeciesAgent)s.getAgent(this)).new ChildAtomIterator();
    }
//    public final AtomIterator atomIterator;
    public final AtomIterator moleculeIterator;
    
/*    public static void main(String[] args) {
        
        Simulation.instance = new etomica.graphics.SimulationGraphic(); 
        new Controller();
        Phase phase = new Phase();
        Species species = new SpeciesSpheresMono();
        new IntegratorHard();
        Potential2 potential = new P2HardSphere();
        potential.setSpecies(species,species);
        new etomica.graphics.DisplayPhase();
//        AtomIterator integrator = phase.makeMoleculeIterator();
        Simulation.instance.elementCoordinator.go();
        etomica.graphics.SimulationGraphic.makeAndDisplayFrame(Simulation.instance);
    }*/
    
} //end of Phase
        