package etomica;

import java.util.Observable;
import java.util.Observer;

import etomica.action.PhaseInflate;
import etomica.atom.AtomLinker;
import etomica.atom.iterator.AtomIteratorListTabbed;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.RectangularLattice;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
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
 * defines the volume of the phase and the behavior of atoms as they move into or
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
  
public class Phase implements EtomicaElement {
        
    private Boundary boundary;
    public final SpeciesMaster speciesMaster;
    private boolean lrcEnabled = true;
    public final SimulationEventManager boundaryEventManager = new SimulationEventManager();
    private String name;
    protected final Space space;
    private PhaseCellManager cellManager;
    
    /**
     * Constructs phase with default rectangular periodic boundary.
     */
    public Phase(Simulation sim) {
        space = sim.space;
        speciesMaster = new SpeciesMaster(sim, this);
        speciesMaster.node.setParent(sim.speciesRoot);
        makeMolecules();
        setName(NameMaker.makeName(this.getClass()));

        setBoundary(new BoundaryRectangularPeriodic(space));

        if(space.D() < 3) 
            setConfiguration(new ConfigurationSequential(space));  //default configuration
        else
            setConfiguration(new ConfigurationLattice(new LatticeCubicFcc()));

    }//end of constructor

    public void setCellManager(PhaseCellManager manager) {
        cellManager = manager;
    }
    
    public PhaseCellManager getCellManager() {
        return cellManager;
    }
    
    public void makeMolecules() {
        AtomLinker agentHead = speciesMaster.node.childList.header;
        for (AtomLinker link=agentHead.next; link!=agentHead; link=link.next) {
            ((SpeciesAgent)link.atom).makeMolecules();
        }
    }
    
    public int getIndex() {
        return speciesMaster.node.getOrdinal();
    }
    
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
        return ((AtomTreeNodeGroup)s.node).childList.get(i-(sum-((AtomTreeNodeGroup)s.node).childAtomCount()));
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
    	AtomIterator iterator = new AtomIteratorListTabbed(speciesMaster.atomList);
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
    public final Boundary boundary() {return boundary;}
    
    public final void setDimensions(Vector d) {
        boundary.setDimensions(d);
        if (cellManager != null) {
            cellManager.getLattice().setDimensions(d);
        }
    }
    
    /**
     * Returns the agent of the given species in this phase.
     */
    public final SpeciesAgent getAgent(Species s) {return s.getAgent(this);}
                       
    public final double volume() {return boundary.volume();}  //infinite volume unless using PBC
    
    public void setDensity(double rho) {
        double vNew = moleculeCount()/rho;
        double scale = Math.pow(vNew/boundary.volume(), 1.0/space.D());
        PhaseInflate inflater = new PhaseInflate(this);
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
        
    public int moleculeCount() {return speciesMaster.moleculeCount();}
    
    public int atomCount() {return speciesMaster.node.leafAtomCount();}
        
    public void setConfiguration(Configuration c) {
        configuration = c;
        if (c != null) {
            c.initializeCoordinates(this);
        }
    }
    
    public Configuration getConfiguration() {return configuration;}
    
    /**
     * Adds the given molecule to this phase.
     * @param molecule the molecule to be added
     * @param s the species agent in this phase for the molecule's species.
     */
    public void addMolecule(Atom a, SpeciesAgent s) {
        a.node.setParent((AtomTreeNodeGroup)s.node);
    }
    
    /**
     * Removes molecule from phase.  
     */
    public void removeMolecule(Atom a) {
        if(a == null) return;
        a.node.dispose();
    }
    
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
     * @return Returns the lattice.
     */
    public RectangularLattice getLattice() {
        return cellManager.getLattice();
    }
    
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

} //end of Phase
        