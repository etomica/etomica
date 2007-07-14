package etomica.atom.iterator;

import etomica.atom.AtomAddressManager;
import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.iterator.IteratorDirective.Direction;
import etomica.box.Box;
import etomica.species.Species;

/**
 * Iterator for the molecules of a single species in a box.  Can be targeted to
 * give the molecule containing a specific atom (if consistent with the species),
 * or all molecules of the species in a box.<br>
 * This class is used by PotentialMaster to iterate over molecules for single-body
 * potentials.
 */
public class AtomIteratorMolecule extends AtomIteratorAdapter implements
        AtomsetIteratorPDT, AtomIteratorBoxDependent {

    /**
     * @param species species for which molecules are returned as iterates.
     * species must not be null.
     */
    public AtomIteratorMolecule(Species species) {
        super(new AtomIteratorArrayListSimple());
        listIterator = (AtomIteratorArrayListSimple)iterator;
        this.species = species;
    }

    /**
     * Sets the box containing the molecules for iteration. A null
     * box conditions iterator to give no iterates.
     * @throws a NullPointerException if the Box is null
     */
    public void setBox(Box newBox) {
        box = newBox;
    }
    
    public void reset() {
        setList();
        super.reset();
    }

    /**
     * Specifies molecule returned by this iterator, as the one containing
     * the given target atom.  Only the first element of the array is relevant.
     * If argument is null, of zero length, or if targetAtom[0] is null, then
     * no target atom is specified and all molecules of the species in the
     * box will be given on iteration.
     * 
     * @throws NullPointerException
     *          if targetAtom is null
     * @throws IllegalArgumentException
     *          if targetAtom.count() is not 0 or 1
     */
    public void setTarget(IAtom newTargetAtom) {
        targetAtom = newTargetAtom;
    }

    /** 
     * Has no effect, but is included as part of the AtomsetIteratorPDT interface.
     */
    public void setDirection(Direction direction) {
        //ignore
    }

    /**
     * Configures the list iterator with a list appropriate to the specified
     * box and target.
     */
    private void setList() {
        if(targetAtom == null) {
            listIterator.setList(box.getMoleculeList(species));
        
        //target specified -- give it as only iterate if descended from species
        } else {
            IAtom molecule = targetAtom;
            if (molecule.getType().getDepth() > AtomAddressManager.MOLECULE_DEPTH) {
                molecule = molecule.getParentGroup();
            }
            if (molecule.getType().getSpecies() != species) {
                molecule = null;
            }
            littleList.clear();
            if(molecule != null) littleList.add(molecule);
            listIterator.setList(littleList);
        }
    }

    private static final long serialVersionUID = 1L;
    private final AtomIteratorArrayListSimple listIterator;
    private final Species species;
    private final AtomArrayList littleList = new AtomArrayList();
    private Box box;
    private IAtom targetAtom;
}
