package etomica.association;
import etomica.*;

/**
 * Class to define and track atom associations.  Constructed given an iterator
 * that defines the set of atoms that is managed, and an association definition
 * that is used to determine if atoms are associated.  Only pairwise associations
 * are considered.
 */
public class AssociationManager /*listener interface*/ {
    
    private AssociationDefinition associationDefinition;
    private final int index = Atom.requestAtomListIndex();
    private AtomIterator atomIterator;
    private final AtomList associatedAtoms = new AtomList();
    private int associatedAtomCount;
    
    public AssociationManager(AtomIterator iterator, AssociationDefinition definition) {
        associationDefinition = definition;
        atomIterator = iterator;
    }
    
    public void initialize() {
        AtomPairIterator pairIterator = new AtomPairIterator(
            Simulation.instance.space, atomIterator, new AtomIteratorList(atomIterator));
        pairIterator.reset();
        while(pairIterator.hasNext()) {
            AtomPair pair = pairIterator.next();
            if(associationDefinition.isAssociated(pair.atom1(), pair.atom2())) {
                pair.atom1().atomList[index].add(pair.atom2());
                pair.atom2().atomList[index].add(pair.atom1());
            }//end if
        }//end while
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom atom = atomIterator.next();
            if(atom.atomList[index] != null && atom.atomList[index].size() != 0) 
                associatedAtoms.add(atom);
        }
    }
    
//    public Atom randomAssociatedAtom() {return associatedAtoms.randomAtom();}
    
    public int associatedAtomCount() {return associatedAtoms.size();}
    
    /**
     * Updates the association status of the given atom.
     */
    public void atomMoveNotify(Atom atom) {
        //loop over all atoms and see if association status with given atom has changed
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom atomB = atomIterator.next();
            if(atomB == atom) continue;
            boolean wasAssociated = false;
            //see if B is on the list of associated atoms for A
            wasAssociated = atom.atomList[index].contains(atomB);
            boolean isAssociated = associationDefinition.isAssociated(atom, atomB);
            
            if(isAssociated != wasAssociated) {//something changed -- one is true and the other is false
                if(wasAssociated) { //but not associated now
                    atom.atomList[index].remove(atomB);
                    atomB.atomList[index].remove(atom);
                    if(atom.atomList[index].size() == 0) associatedAtoms.remove(atom);
                    if(atomB.atomList[index].size() == 0) associatedAtoms.remove(atomB);
                }
                else { //is associated now but wasn't before
                    if(atom.atomList[index].size() == 0) associatedAtoms.add(atom);
                    if(atomB.atomList[index].size() == 0) associatedAtoms.add(atomB);
                    atom.atomList[index].add(atomB);
                    atomB.atomList[index].add(atom);
                }
            }//end if(isAssociated != wasAssociated)
        }//end while
    }//end update
    
    /**
     * Returns the number of atoms on the list of associations of the given atom.
     */
    public int associationCount(Atom atom) {
       return atom.atomList[index].size();
    }
    
    public AtomIterator iterator() {return atomIterator;}
    
    /**
     * Returns the index of the atomList array where each atom holds references to the
     * other atoms that are currently associated with it.
     */
    public int associationListIndex() {return index;}
    
    /**
     * Class to replace the default EnergySum used to compute the total energy.
     * Replaces the default pair-energy calculation method with one that first
     * updates the association status of each iterated pair before calling the energy method.
     */
    public class EnergySum extends PotentialCalculationEnergySum {
        public void calculate(AtomPairIterator iterator, Potential2 potential) {
            while(iterator.hasNext()) {
                AtomPair pair = iterator.next();
                Atom atom = pair.atom1;
                Atom atomB = pair.atom2;

                boolean wasAssociated = atom.atomList[index].contains(atomB);
                boolean isAssociated = associationDefinition.isAssociated(atom, atomB);
                
                if(isAssociated != wasAssociated) {//something changed -- one is true and the other is false
                    if(wasAssociated) { //but not associated now
                        atom.atomList[index].remove(atomB);
                        atomB.atomList[index].remove(atom);
                        if(atom.atomList[index].size() == 0) associatedAtoms.remove(atom);
                        if(atomB.atomList[index].size() == 0) associatedAtoms.remove(atomB);
                    }
                    else { //is associated now but wasn't before
                        if(atom.atomList[index].size() == 0) associatedAtoms.add(atom);
                        if(atomB.atomList[index].size() == 0) associatedAtoms.add(atomB);
                        atom.atomList[index].add(atomB);
                        atomB.atomList[index].add(atom);
                    }
                }//end if(isAssociated != wasAssociated)
                
                sum += potential.energy(pair);
                if(sum >= Double.MAX_VALUE) return;
            }//end while
        }//end of calculate
    }//end of EnergySum inner class
}
    