package etomica.box;

import etomica.action.BoxInflate;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomSet;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IBoxEventManager;
import etomica.api.IMolecule;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.atom.AtomArrayList;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.ISpace;
import etomica.util.Arrays;
import etomica.util.Debug;

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
public class Box implements java.io.Serializable, IBox {

    /**
     * Constructs box with default rectangular periodic boundary.
     */
    public Box(ISimulation sim, ISpace _space) {
        this(new BoundaryRectangularPeriodic(sim.getRandom(), _space), _space);
    }
    
    /**
     * Constructs box with the given boundary
     */
    public Box(IBoundary boundary, ISpace _space) {
    	this.space = _space;
        eventManager = new BoxEventManager();
        setBoundary(boundary);
        
        moleculeLists = new AtomArrayList[0];
        allMoleculeList = new AtomSetAllMolecules();
        allMoleculeList.setMoleculeLists(moleculeLists);
        leafList = new AtomArrayList();
        
        indexReservoir = new int[reservoirSize];
        maxIndex = -1;
        reservoirCount = 0;
    }
    
    public void setIndex(int newIndex) {
        index = newIndex;
    }

    public int getIndex() {
        return index;
    }
    
    public String toString() {
        return "Box"+getIndex();
    }
    
    public IMolecule addNewMolecule(ISpecies species) {
        IMolecule aNew = species.makeMolecule();
        addMolecule(aNew);
        return aNew;
    }

    public void addMolecule(IMolecule molecule) {
        int speciesIndex = molecule.getType().getIndex();
        if (Debug.ON) {
            for (int i=0; i<moleculeLists[speciesIndex].getAtomCount(); i++) {
                if (moleculeLists[speciesIndex].getAtom(i) == molecule) {
                    throw new RuntimeException("you bastard!");
                }
            }
        }
        molecule.setIndex(moleculeLists[speciesIndex].getAtomCount());
        moleculeLists[speciesIndex].add(molecule);
        allMoleculeList.setMoleculeLists(moleculeLists);

        IAtomSet childList = molecule.getChildList();
        int nLeafAtoms = leafList.getAtomCount();
        for (int iChild = 0; iChild < childList.getAtomCount(); iChild++) {
            IAtomLeaf childAtom = (IAtomLeaf)childList.getAtom(iChild);
            childAtom.setLeafIndex(nLeafAtoms++);
            leafList.add(childAtom);
        }
        eventManager.fireEvent(new BoxAtomAddedEvent(this, molecule));

        if (Debug.ON) {
            for (int i=0; i<moleculeLists[speciesIndex].getAtomCount(); i++) {
                if (moleculeLists[speciesIndex].getAtom(i).getIndex() != i) {
                    throw new RuntimeException("oops "+molecule+" "+moleculeLists[speciesIndex].getAtom(i)+" "+i);
                }
            }
        }
    }

    public void removeMolecule(IMolecule molecule) {
        int moleculeIndex = molecule.getIndex();
        AtomArrayList moleculeList = moleculeLists[molecule.getType().getIndex()];
        if (Debug.ON && moleculeList.getAtom(moleculeIndex) != molecule) {
            throw new IllegalArgumentException("can't find "+molecule);
        }
        if (moleculeIndex < moleculeList.getAtomCount()-1) {
            moleculeList.removeAndReplace(moleculeIndex);
            IMolecule replacingMolecule = (IMolecule)moleculeList.getAtom(moleculeIndex);
            replacingMolecule.setIndex(moleculeIndex);
            BoxMoleculeIndexChangedEvent event = new BoxMoleculeIndexChangedEvent(this, replacingMolecule, moleculeList.getAtomCount());
            eventManager.fireEvent(event);
        }
        else {
            moleculeList.remove(moleculeIndex);
        }
        allMoleculeList.setMoleculeLists(moleculeLists);

        eventManager.fireEvent(new BoxAtomRemovedEvent(this, molecule));
        IAtomSet childList = molecule.getChildList();
        for (int iChild = 0; iChild < childList.getAtomCount(); iChild++) {
            IAtomLeaf childAtom = (IAtomLeaf)childList.getAtom(iChild);
            int leafIndex = childAtom.getLeafIndex();
            leafList.removeAndReplace(leafIndex);
            if (leafList.getAtomCount() > leafIndex) {
                IAtomLeaf movedAtom = (IAtomLeaf)leafList.getAtom(leafIndex);
                int movedLeafIndex = movedAtom.getLeafIndex();
                BoxAtomLeafIndexChangedEvent event = new BoxAtomLeafIndexChangedEvent(this, movedAtom, movedLeafIndex);
                movedAtom.setLeafIndex(leafIndex);
                eventManager.fireEvent(event);
            }
        }
        leafList.maybeTrimToSize();
    }
    
    public void setNMolecules(ISpecies species, int n) {
        int speciesIndex = species.getIndex();
        AtomArrayList moleculeList = moleculeLists[speciesIndex];
        int currentNMolecules = moleculeList.getAtomCount();
        notifyNewMolecules(species, (n-currentNMolecules));
        if(n < 0) {
            throw new IllegalArgumentException("Number of molecules cannot be negative");
        }
        if (n > currentNMolecules) {
            moleculeLists[species.getIndex()].ensureCapacity(n);
            leafList.ensureCapacity(leafList.getAtomCount()+(n-currentNMolecules)*species.getNumLeafAtoms());
            for(int i=currentNMolecules; i<n; i++) {
                addMolecule(species.makeMolecule());
            }
        }
        else {
            for (int i=currentNMolecules; i>n; i--) {
                removeMolecule((IMolecule)moleculeList.getAtom(i-1));
            }
        }
    }
    
    public int getNMolecules(ISpecies species) {
        int speciesIndex = species.getIndex();
        return moleculeLists[speciesIndex].getAtomCount();
    }
    
    public IAtomSet getMoleculeList(ISpecies species) {
        return moleculeLists[species.getIndex()];
    }

    public IAtomSet getMoleculeList() {
        return allMoleculeList;
    }
  
    public void setBoundary(IBoundary b) {
        boundary = b;
        boundary.setBox(this);
     }
     
    public final IBoundary getBoundary() {return boundary;}
    
    public void setDensity(double rho) {
        double vNew = getMoleculeList().getAtomCount()/rho;
        double scale = Math.pow(vNew/boundary.volume(), 1.0/space.D());
        BoxInflate inflater = new BoxInflate(this, space);
        inflater.setScale(scale);
        inflater.actionPerformed();
    }

    public IBoxEventManager getEventManager() {
        return eventManager;
    }

    public void addSpeciesNotify(ISpecies species) {
        moleculeLists = (AtomArrayList[])Arrays.addObject(moleculeLists, new AtomArrayList());
        allMoleculeList.setMoleculeLists(moleculeLists);
    }
    
    public void removeSpeciesNotify(ISpecies species) {
        moleculeLists = (AtomArrayList[])Arrays.removeObject(moleculeLists, moleculeLists[species.getIndex()]);
        allMoleculeList.setMoleculeLists(moleculeLists);
    }

    public IAtomSet getLeafList() {
        return leafList;
    }
    
    /**
     * Notifies the SpeciesMaster that the given number of new Atoms will be
     * added to the system.  It's not required to call this method before
     * adding atoms, but if adding many Atoms, calling this will improve
     * performance.
     */
    protected void notifyNewMolecules(ISpecies species, int numNewMolecules) {
        if (numNewMolecules < 1) return;
        // has no actual effect within this object.  We just notify things to 
        // prepare for an increase in the max index.  If things assume that the
        // actual max index has already increased, there's no harm since
        // there's nothing that says the max index can't be too large.
        int numNewLeafAtoms = numNewMolecules * species.getNumLeafAtoms();
        if (numNewLeafAtoms + numNewMolecules > reservoirCount) {
            BoxGlobalAtomIndexEvent event = new BoxGlobalAtomIndexEvent(this, maxIndex + numNewMolecules + numNewLeafAtoms - reservoirCount);
            eventManager.fireEvent(event);
        }
        BoxNumMoleculesEvent event = new BoxNumMoleculesEvent(this, species, moleculeLists[species.getIndex()].getAtomCount() + numNewMolecules);
        eventManager.fireEvent(event);
        if (numNewLeafAtoms > 1) {
            BoxGlobalAtomLeafIndexEvent leafEvent = new BoxGlobalAtomLeafIndexEvent(this, leafList.getAtomCount() + numNewLeafAtoms);
            eventManager.fireEvent(leafEvent);
        }
    }
    
    public int getLeafIndex(IAtomLeaf atomLeaf) {
        return atomLeaf.getLeafIndex();
    }

    /**
     * List of leaf atoms in box
     */
    protected final AtomArrayList leafList;

    protected int[] indexReservoir;
    protected int reservoirSize = 50;
    protected int reservoirCount;
    protected int maxIndex;
    
    private static final long serialVersionUID = 2L;
    private IBoundary boundary;;
    private final IBoxEventManager eventManager;
    protected AtomArrayList[] moleculeLists;
    private int index;
    private final ISpace space;
    protected final AtomSetAllMolecules allMoleculeList;
}
