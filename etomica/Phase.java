package etomica;

import java.util.Observable;
import java.util.Observer;

import etomica.action.PhaseImposePbc;
import etomica.action.PhaseInflate;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.RectangularLattice;
import etomica.space.Boundary;
import etomica.space.BoundaryPeriodicSquare;
import etomica.space.Vector;
import etomica.utility.NameMaker;

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
 * A Phase collects all atoms that interact with one another; atoms in different
 * phases do not interact. These are the important features of a Phase:
 * <p>
 * <ol>
 * <li>It holds a SpeciesMaster instance, which provides the root of a
 * hierarchy of atoms that represent the physical object that interact.
 * <li>It holds a Boundary object, obtained from the governing Space, that
 * defines the volume of the phase and the behavio of atoms as they move into or
 * across the boundary of the phase.
 * <li>It has a Configuration object that determines the default initial
 * configuration of the atoms in the phase.
 * <li>It maintains a list of listeners that are informed when significant
 * events happen in the phase (such as a change in its boundary).
 * <li>Each Phase has a unique species index assigned when it is constructed.
 * The index assignment begins at 0 and is incremented after each Phase
 * construction. This index is useful when collecting things in reference to the
 * phase.
 * </ol>
 * A phase acted upon by an Integrator instance to move its atoms around and
 * generate configurations. Properties of a phase are measured by MeterAbstract
 * instances which are simply DataSource objects that require a phase to
 * generate their data. <br>
 * A simulation may involve more than one phase. All Phase instances are
 * normally registered with the default simulation upon their construction, and
 * may be accessed via the simulation's getPhaseList method.
 * 
 * @author David Kofke
 * @see Boundary
 */
 
 /* History of changes
  * 7/3/02  added reset method
  */
  
public class Phase {
        
    private Boundary boundary;
    private PhaseInflate inflater;
    public final SpeciesMaster speciesMaster;
    private boolean lrcEnabled = true;
    public final SimulationEventManager boundaryEventManager = new SimulationEventManager();
    private static int nonSimCount = 0;//number of times instantiated without a parent simulation
    private PhaseImposePbc centralImageEnforcer;
    private String name;
    public final int index;
    private final Space space;
    private static int instanceCount;
    
    /**
     * Constructs phase.
     */
    public Phase(Space space) {
        index = instanceCount++;
        speciesMaster = new SpeciesMaster(space, this);
        this.space = space;
        setName(NameMaker.makeName(this.getClass()));
        inflater = new PhaseInflate(this);

        setBoundary(new BoundaryPeriodicSquare(space));

        if(space.D() < 3) 
            setConfiguration(new ConfigurationSequential(space));  //default configuration
        else
            setConfiguration(new ConfigurationLattice(new LatticeCubicFcc()));

        if (Default.AUTO_REGISTER) {
            Simulation.getDefault().register(this);
        }
    }//end of constructor
   

    /**
     * Accessor method of the name of this phase
     * 
     * @return The given name of this phase
     */
    public final String getName() {return name;}
    /**
     * Method to set the name of this simulation element. The element's name
     * provides a convenient way to label output data that is associated with
     * it.  This method might be used, for example, to place a heading on a
     * column of data. Default name is the base class followed by the integer
     * index of this element.
     * 
     * @param name The name string to be associated with this element
     */
    public void setName(String name) {this.name = name;}

    /**
     * Overrides the Object class toString method to have it return the output of getName
     * 
     * @return The name given to the phase
     */
    public String toString() {return getName();}
    
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
    
    public final Space space() {return space;}
    
    public Vector randomPosition() {return boundary.randomPosition();}
        
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
     * Finds and returns the atom nearest to each of one or more given
     * positions, using the boundary associated with this phase.
     * @param r a positions in the phase.  Point may be outside the phase, as
     * minimum-image separations are used in accordance with phase's boundary.
     * @return Atom the (leaf) atoms in the phase nearest to each position. A
     * given atom can appear more than once in this list, if it is nearest to
     * more than one of the positions.
     */
    public Atom[] nearestAtom(Vector[] r) {
    	AtomIterator iterator = new AtomIteratorListSimple(speciesMaster.atomList);
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
     public void setBoundary(Boundary b) {
        boundaryMonitor.notifyObservers(b);
        boundary = b;
     }
     
     /**
      * Resets phase by initializing coordinates according to current instance of configuration.
      * Fires PhaseEvent of type RESET after completing action.
      */
     public void reset() {
         Configuration c = getConfiguration();
         if (c != null) {
             c.initializeCoordinates(this);
         }
         fireEvent(new PhaseEvent(this, PhaseEvent.RESET));
     }

    /**
     * Returns the current boundary instance.
     * 
     * @return The current instance of the boundary class
     */
    public final Boundary getBoundary() {return boundary;}
    /**
     * Same as getBoundary.
     */
    public final Boundary boundary() {return boundary;}
    
    /**
     * Returns the agent of the given species in this phase.
     */
    public final SpeciesAgent getAgent(Species s) {return s.getAgent(this);}
                       
    public final Vector dimensions() {return boundary.dimensions();}
    public final double volume() {return boundary.volume();}  //infinite volume unless using PBC
    
    public void setDensity(double rho) {
        double vNew = moleculeCount()/rho;
        double scale = Math.pow(vNew/boundary.volume(), 1.0/space.D());
        inflater.setScale(scale);
        inflater.actionPerformed();
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
        if (c != null) {
            c.setZeroTotalMomentum(true);
            c.initializeCoordinates(this);
        }
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
    
    public synchronized void addListener(PhaseListener listener) {
        eventManager.addListener(listener);
    }

    public synchronized void removeListener(PhaseListener listener) {
        eventManager.removeListener(listener);
    }

    protected void fireEvent(PhaseEvent event) {
        eventManager.fireEvent(event);
    }    
     
    
    /**
     * Temporary method for cell neighbor listing
     * @return Returns the lattice.
     */
    public RectangularLattice getLattice() {
        return lattice;
    }
    /**
     * Temporary method for cell neighbor listing
     * @param lattice The lattice to set.
     */
    public void setLattice(RectangularLattice lattice) {
        this.lattice = lattice;
    }
    //temporary construct for cell neighbor listing
    private RectangularLattice lattice;
    
    public Configuration configuration;
          
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
	 * @return PhaseActionAdapter.PhaseImposePbc
	 */
	public PhaseImposePbc getCentralImageEnforcer() {
		return centralImageEnforcer;
	}

	/**
	 * Sets the centralImageEnforcer.  Normally this is set by the
	 * IntegratorPhase mediator, and is also registered as an
	 * IntervalEventListener with the Integrator.
	 * @param centralImageEnforcer The centralImageEnforcer to set
	 */
	public void setCentralImageEnforcer(
		PhaseImposePbc centralImageEnforcer) {
		this.centralImageEnforcer = centralImageEnforcer;
	}

} //end of Phase
        