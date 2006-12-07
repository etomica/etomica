package etomica.phase;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;

import etomica.EtomicaElement;
import etomica.action.PhaseInflate;
import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesAgent;
import etomica.atom.SpeciesMaster;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.ICoordinate;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.Species;
import etomica.species.SpeciesResolver;
import etomica.species.SpeciesSignature;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;
import etomica.units.Quantity;
import etomica.units.Volume;

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
 * hierarchy of atoms that represent the physical objects that interact.
 * <li>It holds a Boundary object, obtained from the governing Space, that
 * defines the volume of the phase and the behavior of atoms as they move into or
 * across the boundary of the phase.
 * <li>It maintains a list of listeners that are informed when significant
 * events happen in the phase (such as a change in its boundary).
 * <li>Each Phase has a unique index assigned when it is constructed.
 * The index assignment begins at 0 and is incremented after each Phase
 * construction. This index is useful when collecting things in reference to the
 * phase.
 * </ol>
 * A phase acted upon by an Integrator instance to move its atoms around and
 * generate configurations. Properties of a phase are measured by MeterAbstract
 * instances which are simply DataSource objects that require a phase to
 * generate their data. <br>
 * A simulation may involve more than one phase. All Phase instances are
 * registered with the simulation specified upon their construction, and
 * may be accessed via the simulation's getPhases method.
 * 
 * @author David Kofke
 * @see Boundary
 */
public class Phase implements EtomicaElement, java.io.Serializable {
        
    /**
     * Constructs phase with default rectangular periodic boundary.
     */
    public Phase(Simulation sim) {
        space = sim.space;
        eventManager = new PhaseEventManager();
        speciesMaster = new SpeciesMaster(sim, this, eventManager);
        setBoundary(new BoundaryRectangularPeriodic(sim));
        speciesMaster.getNode().setParent((AtomTreeNodeGroup)sim.speciesRoot.getNode());
        setName(null);

        inflateEvent = new PhaseInflateEvent(this);
    }

    public int getIndex() {
        return speciesMaster.getNode().getIndex();
    }
    
    /**
     * Accessor method of the name of this phase
     * 
     * @return The given name of this phase
     */
    public final String getName() {
        if (name == null) {
            return "Phase"+getIndex();
        }
        return name;
    }
    
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
        AtomArrayList agentList = ((AtomTreeNodeGroup)speciesMaster.getNode()).childList;
        for (int agentIndex=0; agentIndex<agentList.size(); agentIndex++) {
            AtomArrayList moleculeList = ((AtomTreeNodeGroup)agentList.get(agentIndex).getNode()).childList;
            int count = moleculeList.size();
            if (i < count) {
                return moleculeList.get(i);
            }
            i -= count;
        }
        throw new IllegalStateException("how can this be?!?!?!");
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
    	AtomIterator iterator = new AtomIteratorLeafAtoms(this);
    	Atom[] nearest = new Atom[r.length];
    	for(int i=0; i<r.length; i++) {
	    	double r2Min = Double.MAX_VALUE;
	    	iterator.reset();
	    	while(iterator.hasNext()) {
	    		AtomLeaf atom = (AtomLeaf)iterator.nextAtom();
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
        return molecule(Simulation.random.nextInt(moleculeCount()));
    }
      
    /**
     * Sets the boundary object of the phase.
     */
     public void setBoundary(Boundary b) {
        boundary = b;
     }
     
    /**
     * Returns the current boundary instance.
     * 
     * @return The current instance of the boundary class
     */
    public final Boundary getBoundary() {return boundary;}
    
    public final void setDimensions(Vector d) {
        boundary.setDimensions(d);
        eventManager.fireEvent(inflateEvent);
    }
    
    /**
     * Returns the agent of the given species in this phase.
     */
    public final SpeciesAgent getAgent(Species s) {
        return (SpeciesAgent)((AtomTreeNodeGroup)speciesMaster.getNode()).childList.get(s.getIndex());
    }
                       
    public final double volume() {return boundary.volume();}  //infinite volume unless using PBC
    
    public void setDensity(double rho) {
        double vNew = moleculeCount()/rho;
        double scale = Math.pow(vNew/boundary.volume(), 1.0/space.D());
        PhaseInflate inflater = new PhaseInflate(this);
        inflater.setScale(scale);
        inflater.actionPerformed();
    }
    
    public double getDensity() {return moleculeCount()/boundary.volume();}
    
    public Dimension getDensityDimension() {
        return new DimensionRatio("Density",Quantity.DIMENSION,Volume.DIMENSION);
    }

    /**
     * @return the first atom in the linked list of atoms in this Phase
     */
    public final Atom firstAtom() {
        return speciesMaster.getNode().firstLeafAtom();
    }
    
    /**
     * @return the last atom in the linked list of atoms in this Phase
     */
    public final Atom lastAtom() {
        return speciesMaster.getNode().lastLeafAtom();
    }

    /**
     * returns the number of molecules in the phase
     */
    public int moleculeCount() {return speciesMaster.moleculeCount();}
    
    /**
     * returns the number of leaf atoms in the phase
     */
    public int atomCount() {return speciesMaster.getNode().leafAtomCount();}
        
    /**
     * Adds the given molecule to this phase.
     * @param molecule the molecule to be added
     * @param s the species agent in this phase for the molecule's species.
     */
    public void addMolecule(Atom a, SpeciesAgent s) {
        a.getNode().setParent((AtomTreeNodeGroup)s.getNode());
    }
    
    /**
     * Removes molecule from phase.  
     */
    public void removeMolecule(Atom a) {
        if(a == null) return;
        a.getNode().dispose();
    }
    
    /**
     * @return Returns the speciesMaster.
     */
    public SpeciesMaster getSpeciesMaster() {
        return speciesMaster;
    }

    public void writePhase(ObjectOutputStream out) throws IOException {
        out.writeObject(boundary);
        AtomArrayList agents = ((AtomTreeNodeGroup)speciesMaster.getNode()).childList;
        out.writeInt(agents.size());
        for (int i=0; i<agents.size(); i++) {
            Species species = agents.get(i).getType().getSpecies();
            out.writeObject(species.getSpeciesSignature());
            out.writeInt(((SpeciesAgent)agents.get(i)).getNMolecules());
        }
        AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
        iterator.setPhase(this);
        iterator.reset();
        while (iterator.hasNext()) {
            out.writeObject(((AtomLeaf)iterator.nextAtom()).coord);
        }
    }
    
    public static Phase readPhase(ObjectInputStream in, Simulation sim, SpeciesResolver resolver) throws IOException, ClassNotFoundException {
        Boundary newBoundary = (Boundary)in.readObject();
        Phase newPhase = new Phase(sim);
        newPhase.setBoundary(newBoundary);
        
        Species[] mySpecies = sim.getSpecies();
        int numSpecies = in.readInt();
        for (int i = 0; i<numSpecies; i++) {
            SpeciesSignature speciesSignature = (SpeciesSignature)in.readObject();
            Species newSpecies = null;
            Species[] candidates = new Species[0];
            for (int j=0; j<mySpecies.length; j++) {
                Species candidate = mySpecies[j];
                if (speciesSignature.equals(candidate.getSpeciesSignature())) {
                    candidates = (Species[])etomica.util.Arrays.addObject(candidates,candidate);
                    break;
                }
                if (candidates.length > 0) {
                    newSpecies = resolver.whichOneDoYouLike(candidates,speciesSignature.name);
                }
            }
            if (newSpecies == null) {
                Constructor constructor = speciesSignature.constructor;
                Object[] parameters = new Object[speciesSignature.parameters.length+1];
                parameters[0] = sim;
                System.arraycopy(speciesSignature.parameters,0,parameters,1,parameters.length-1);
                try {
                    newSpecies = (Species)constructor.newInstance(parameters);
                }
                catch (IllegalAccessException e) {}
                catch (InstantiationException e) {}
                catch (InvocationTargetException e) {}
            }
            int nMolecules = in.readInt();
            newPhase.getAgent(newSpecies).setNMolecules(nMolecules);
        }

        AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
        iterator.setPhase(newPhase);
        iterator.reset();
        while (iterator.hasNext()) {
            ((AtomLeaf)iterator.nextAtom()).coord.E((ICoordinate)in.readObject());
        }
        return newPhase;
    }
    
    public PhaseEventManager getEventManager() {
        return eventManager;
    }

    private static final long serialVersionUID = 1L;
    private Boundary boundary;
    private SpeciesMaster speciesMaster;
    private boolean lrcEnabled = true;
    private String name;
    protected final Space space;
    private final PhaseEventManager eventManager;
    private final PhaseEvent inflateEvent;
} //end of Phase
        
