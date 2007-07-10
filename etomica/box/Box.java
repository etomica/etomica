package etomica.box;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;

import etomica.action.BoxInflate;
import etomica.atom.AtomManager;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.IAtomPositioned;
import etomica.atom.ISpeciesAgent;
import etomica.simulation.ISimulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.species.Species;
import etomica.species.SpeciesResolver;
import etomica.species.SpeciesSignature;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;
import etomica.units.Quantity;
import etomica.units.Volume;

/**
 * A Box collects all atoms that interact with one another; atoms in different
 * boxs do not interact. These are the important features of a Box:
 * <p>
 * <ol>
 * <li>It holds a SpeciesMaster instance, which provides the root of a
 * hierarchy of atoms that represent the physical objects that interact.
 * <li>It holds a Boundary object, obtained from the governing Space, that
 * defines the volume of the box and the behavior of atoms as they move into or
 * across the boundary of the box.
 * <li>It maintains a list of listeners that are informed when significant
 * events happen in the box (such as a change in its boundary).
 * <li>Each Box has a unique index assigned when it is constructed.
 * The index assignment begins at 0 and is incremented after each Box
 * construction. This index is useful when collecting things in reference to the
 * box.
 * </ol>
 * A box is acted upon by an Integrator instance to move its atoms around and
 * generate configurations. Properties of a box are measured by MeterAbstract
 * instances which are simply DataSource objects that require a box to
 * generate their data. <br>
 * A simulation may involve more than one box. All Box instances are
 * registered with the simulation specified upon their construction, and
 * may be accessed via the simulation's getBoxs method.
 * 
 * @author David Kofke, Andrew Schultz
 * @see Boundary
 */
public class Box implements java.io.Serializable {
        
    /**
     * Constructs box with default rectangular periodic boundary.
     */
    public Box(ISimulation sim) {
        this(new BoundaryRectangularPeriodic(sim));
    }
    
    /**
     * Constructs box with the given boundary
     */
    public Box(Boundary boundary) {
        space = boundary.getSpace();
        eventManager = new BoxEventManager();
        atomManager = new AtomManager(this, eventManager);
        setBoundary(boundary);
        
        inflateEvent = new BoxInflateEvent(this);
    }
    
    /**
     * Resets the Box's index.  This should only need to be called from the
     * Simulation class.
     * 
     * @param sim  The Simulation to which this Box was added.  Passing null
     *             notifies the Box that it was removed from the Simulation
     *             (the index is set to 0).
     */
    public void resetIndex(ISimulation sim) {
        if (sim == null) {
            // sim is notifying us that we got removed.
            index = 0;
            return;
        }
        Box[] boxs = sim.getBoxs();
        for (int i=0; i<boxs.length; i++) {
            if (boxs[i] == this) {
                index = i;
                return;
            }
        }
        // you really shouldn't be calling resetIndex unless you're a simulation!
        throw new IllegalArgumentException(sim+" does not contain me");
    }

    public int getIndex() {
        return index;
    }
    
    /**
     * Overrides the Object class toString method to have it return the output of getName
     * 
     * @return The name given to the box
     */
    public String toString() {
        return "Box"+getIndex();
    }
    
    public final Space getSpace() {return space;}
    
    public IAtom addNewMolecule(Species species) {
        IAtom aNew = species.getMoleculeFactory().makeAtom();
        getAgent(species).addChildAtom(aNew);
        return aNew;
    }
    
    public void addMolecule(IAtom molecule) {
        Species species = molecule.getType().getSpecies();
        getAgent(species).addChildAtom(molecule);
    }

    public void removeMolecule(IAtom molecule) {
        Species species = molecule.getType().getSpecies();
        getAgent(species).removeChildAtom(molecule);
    }
    
    /**
     * Sets the number of molecules for this species.  Molecules are either
     * added or removed until the given number is obtained.  Takes no action
     * at all if the new number of molecules equals the existing number.
     *
     * @param n  the new number of molecules for this species
     */
    public void setNMolecules(Species species, int n) {
        ISpeciesAgent agent = getAgent(species);
        int currentNMolecules = getAgent(species).getNMolecules();
        atomManager.notifyNewAtoms((n-currentNMolecules)*species.getMoleculeFactory().getNumTreeAtoms(),
                                     (n-currentNMolecules)*species.getMoleculeFactory().getNumLeafAtoms());
        if(n > agent.getChildList().getAtomCount()) {
            for(int i=agent.getChildList().getAtomCount(); i<n; i++) addNewMolecule(species);
        }
        else if(n < agent.getChildList().getAtomCount()) {
            if(n < 0) {
                throw new IllegalArgumentException("Number of molecules cannot be negative");
            }
            for (int i=agent.getChildList().getAtomCount(); i>n; i--) {
                removeMolecule(agent.getChildList().getAtom(i-1));
            }
        }
    }
    
    public int getNMolecules(Species species) {
        return getAgent(species).getNMolecules();
    }
    
    public AtomSet getMoleculeList(Species species) {
        return getAgent(species).getChildList();
    }
    
    /**
     * Returns the ith molecule in the linked list of molecules.
     * 0 returns the first molecule, and moleculeCount-1 returns the last.
     * An argument outside this range throws an IndexOutOfBoundsException
     */
    //could make more efficient by starting from first or last molecule, as appropriate
    public IAtom molecule(int i) {
        if(i >= moleculeCount() || i < 0) 
            throw new IndexOutOfBoundsException("Index: "+i+
                                                ", Number of molecules: "+moleculeCount());
        AtomSet agentList = atomManager.getAgentList();
        for (int agentIndex=0; agentIndex<agentList.getAtomCount(); agentIndex++) {
            AtomSet moleculeList = ((IAtomGroup)agentList.getAtom(agentIndex)).getChildList();
            int count = moleculeList.getAtomCount();
            if (i < count) {
                return moleculeList.getAtom(i);
            }
            i -= count;
        }
        throw new IllegalStateException("how can this be?!?!?!");
    }
    
    /**
     * Sets the boundary object of the box.
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
    
    public final void setDimensions(IVector d) {
        boundary.setDimensions(d);
        eventManager.fireEvent(inflateEvent);
    }
    
    /**
     * Returns the agent of the given species in this box.
     */
    public final ISpeciesAgent getAgent(Species s) {
        //brute force it
        AtomSet agentList = atomManager.getAgentList();
        for (int i=0; i<agentList.getAtomCount(); i++) {
            if (((ISpeciesAgent)agentList.getAtom(i)).getType().getSpecies() == s) {
                return (ISpeciesAgent)agentList.getAtom(i);
            }
        }
        return null;
    }

    public final double volume() {return boundary.volume();}  //infinite volume unless using PBC
    
    public void setDensity(double rho) {
        double vNew = moleculeCount()/rho;
        double scale = Math.pow(vNew/boundary.volume(), 1.0/space.D());
        BoxInflate inflater = new BoxInflate(this);
        inflater.setScale(scale);
        inflater.actionPerformed();
    }
    
    public double getDensity() {return moleculeCount()/boundary.volume();}
    
    public Dimension getDensityDimension() {
        return new DimensionRatio("Density",Quantity.DIMENSION,Volume.DIMENSION);
    }

    public AtomSet getLeafList() {
        return atomManager.getLeafList();
    }
    
    /**
     * returns the number of molecules in the box
     */
    public int moleculeCount() {return atomManager.moleculeCount();}
    
    /**
     * returns the number of leaf atoms in the box
     */
    public int atomCount() {return atomManager.getLeafList().getAtomCount();}

    /**
     * @return Returns the speciesMaster.
     */
    public AtomManager getSpeciesMaster() {
        return atomManager;
    }

    public void writeBox(ObjectOutputStream out) throws IOException {
        out.writeObject(boundary);
        AtomSet agents = atomManager.getAgentList();
        out.writeInt(agents.getAtomCount());
        for (int i=0; i<agents.getAtomCount(); i++) {
            Species species = agents.getAtom(i).getType().getSpecies();
            out.writeObject(species.getSpeciesSignature());
            out.writeInt(((ISpeciesAgent)agents.getAtom(i)).getNMolecules());
        }
        AtomSet leafList = atomManager.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomPositioned a = (IAtomPositioned)leafList.getAtom(iLeaf);
            out.writeObject(a.getPosition());
        }
    }
    
    public static Box readBox(ObjectInputStream in, ISimulation sim, SpeciesResolver resolver) throws IOException, ClassNotFoundException {
        Boundary newBoundary = (Boundary)in.readObject();
        Box newBox = new Box(sim);
        sim.addBox(newBox);
        newBox.setBoundary(newBoundary);
        
        Species[] mySpecies = sim.getSpeciesManager().getSpecies();
        int numSpecies = in.readInt();
        for (int i = 0; i<numSpecies; i++) {
            SpeciesSignature speciesSignature = (SpeciesSignature)in.readObject();
            Species newSpecies = null;
            Species[] candidates = new Species[0];
            for (int j=0; j<mySpecies.length; j++) {
                Species candidate = mySpecies[j];
                if (speciesSignature.equals(candidate.getSpeciesSignature())) {
                    candidates = (Species[])etomica.util.Arrays.addObject(candidates,candidate);
                }
            }
            if (candidates.length > 0) {
                newSpecies = resolver.whichOneDoYouLike(candidates);
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
            newBox.setNMolecules(newSpecies, nMolecules);
        }

            //XXX broken
        // loop over the atoms
//            ((AtomLeaf)iterator.nextAtom()).getCoord().E((ICoordinate)in.readObject());
        return newBox;
    }
    
    public BoxEventManager getEventManager() {
        return eventManager;
    }

    private static final long serialVersionUID = 2L;
    private Boundary boundary;
    private AtomManager atomManager;
    protected final Space space;
    private final BoxEventManager eventManager;
    private final BoxEvent inflateEvent;
    private int index;
}
        
