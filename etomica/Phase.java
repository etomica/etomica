package etomica;

import java.util.Observable;
import java.util.Observer;

//Java2 imports
//import java.util.LinkedList;
//import java.util.Iterator;

//import etomica.utility.java2.LinkedList;
//import etomica.utility.java2.Iterator;

/* History of changes
 * 09/01/02 (DAK) setConfiguration sets new configuration so that it zeros total momentum
 *                when used.  Change made while modify behavior of momentum initiation
 *                in Configuration and randomizeMomentum methods in Space.Coord...
 * 01/21/04 (DAK) changed initializeCoordinate calls to take this phase as
 * argument.  As a result, Configuration will set its dimensions equal to that
 * of this phase before assigning coordinates.
 * 01/29/04 (DAK) added nearestAtom method.
 */

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
 * @see Space.Boundary
 * @see MeterAbstract
 */
 
 /* History of changes
  * 7/3/02  added reset method
  */
  
public class Phase extends SimulationElement {
        
    private Space.Boundary boundary;
    private PhaseAction.Inflate inflater;
    public final SpeciesMaster speciesMaster;
    private boolean lrcEnabled = true;
    public final SimulationEventManager boundaryEventManager = new SimulationEventManager();
    private static int nonSimCount = 0;//number of times instantiated without a parent simulation
    private PhaseAction.ImposePbc centralImageEnforcer;
    
    /**
     * Constructs phase and registers as part of the current Simulation.instance.
     */
    public Phase() {
        this(Simulation.instance);
    }
    
    /**
     * Constructs phase and registers it as part of the given simulation.
     */
    public Phase(SimulationElement parent) {
        super(parent, Phase.class);
        
        speciesMaster = new SpeciesMaster(this);
                
        inflater = new PhaseAction.Inflate(this);

//        setBoundary(Space.Boundary.DEFAULT);
        setBoundary(space.makeBoundary());

        if(space.D() < 3) 
            setConfiguration(new ConfigurationSequential(parent.simulation()));  //default configuration
        else
            setConfiguration(new ConfigurationFcc(parent.simulation()));
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
        return ((AtomTreeNodeGroup)s.node).getAtom(i-(sum-((AtomTreeNodeGroup)s.node).childAtomCount()));
    }
    
    /**
     * Finds and returns the atom nearest to the each of one or more given
     * positions, using the boundary associated with this phase.
     * @param r a positions in the phase.  Point may be outside the phase, as
     * minimum-image separations are used in accordance with phase's boundary.
     * @return Atom the (leaf) atoms in the phase nearest to each position. A
     * given atom can appear more than once in this list, if it is nearest to
     * more than one of the positions.
     */
    public Atom[] nearestAtom(Space.Vector[] r) {
    	AtomIterator iterator = makeAtomIterator();
    	Atom[] nearest = new Atom[r.length];
    	for(int i=0; i<r.length; i++) {
	    	double r2Min = Double.MAX_VALUE;
	    	iterator.reset();
	    	while(iterator.hasNext()) {
	    		Atom atom = iterator.nextAtom();
	    		double r2 = Space.r2(atom.coord.position(), r[i], boundary);
	    		if(r2 < r2Min) {
	    			r2Min = r2;
	    			nearest[i] = atom;
	    		}
	    	}
    	}
    	return nearest;
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
      * Resets phase by initializing coordinates according to current instance of configuration.
      * Fires PhaseEvent of type RESET after completing action.
      */
     public void reset() {
        getConfiguration().initializeCoordinates(this);
        fireEvent(new PhaseEvent(this, PhaseEvent.RESET));
     }

    /**
     * Returns the current boundary instance.
     * 
     * @return The current instance of the boundary class
     */
    public final Space.Boundary getBoundary() {return boundary;}
    /**
     * Same as getBoundary.
     */
    public final Space.Boundary boundary() {return boundary;}
    
    /**
     * Returns the agent of the given species in this phase.
     */
    public final SpeciesAgent getAgent(Species s) {return s.getAgent(this);}
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
        double scale = Math.pow(vNew/boundary.volume(), 1.0/space.D());
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
        
    public void setConfiguration(Configuration c) {
        configuration = c;
        c.setZeroTotalMomentum(true);
        configuration.initializeCoordinates(this);
    }
    
    public Configuration getConfiguration() {return configuration;}
    
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
        a.node.setParent((AtomTreeNodeGroup)s.node);
    }
    
    /**
     * Removes molecule from phase.  
     */
    public void removeMolecule(Atom a) {
        if(a == null) return;
        a.sendToReservoir();
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

    public synchronized void addListener(PhaseListener listener) {
        eventManager.addListener(listener);
    }

    public synchronized void removeListener(PhaseListener listener) {
        eventManager.removeListener(listener);
    }

    protected void fireEvent(PhaseEvent event) {
        eventManager.fireEvent(event);
    }    
                                 
    public Configuration configuration;
          
    private Integrator integrator;
              
//    public transient MeterEnergy energy;
        
    public Phase.Monitor integratorMonitor = new Phase.Monitor();
    public Phase.Monitor boundaryMonitor = new Phase.Monitor();
//    public SimulationEventManager integratorMonitor = new SimulationEventManager();

    //used to handle firing of reset events
    private SimulationEventManager eventManager = new SimulationEventManager();
    
    public static class Monitor extends Observable implements java.io.Serializable {
        
        public void notifyObservers() {
            this.notifyObservers(null);
        }
        public void notifyObservers(Object obj) {
            setChanged();
            super.notifyObservers(obj);
        }
        public void addObserver(Observer o) {if(o != null) super.addObserver(o);}
        public void deleteObserver(Observer o) {super.deleteObserver(o);}
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
        return new AtomIteratorList(speciesMaster.atomList);
    }
    
    /**
     * Makes an iterator that loops through all the molecules 
     * (children of the species agents) in this phase.
     */
    public AtomIterator makeMoleculeIterator() {
        return new AtomIteratorMolecule(this);
    }//end of makeMoleculeIterator
    
    /**
     * Makes an iterator that loops through all the (leaf) atoms 
     * derived from the given species in this phase.
     */
    public AtomIterator makeAtomIterator(Species s) {
        return new AtomIteratorTree(getAgent(s));
    }
    
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
    
	/**
	 * Returns the centralImageEnforcer.
	 * @return PhaseAction.ImposePbc
	 */
	public PhaseAction.ImposePbc getCentralImageEnforcer() {
		return centralImageEnforcer;
	}

	/**
	 * Sets the centralImageEnforcer.  Normally this is set by the
	 * IntegratorPhase mediator, and is also registered as an
	 * IntervalEventListener with the Integrator.
	 * @param centralImageEnforcer The centralImageEnforcer to set
	 */
	public void setCentralImageEnforcer(
		PhaseAction.ImposePbc centralImageEnforcer) {
		this.centralImageEnforcer = centralImageEnforcer;
	}

} //end of Phase
        