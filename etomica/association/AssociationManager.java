package etomica.association;
import etomica.*;

/**
 * Class to define and track atom associations.  Constructed given an iterator
 * that defines the set of atoms that is managed, and an association definition
 * that is used to determine if atoms are associated.  Only pairwise associations
 * are considered.
 */
public class AssociationManager implements MCMoveEventListener {
    
    private AssociationDefinition associationDefinition;
    private final int index = Atom.requestAtomListIndex();
    private final int indexOld = Atom.requestAtomListIndex();
    private AtomIterator atomIterator;
    private final AtomIteratorList listIterator0 = new AtomIteratorList();
    private final AtomIteratorList listIterator1 = new AtomIteratorList();
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
    
    //need also to handle associatedAtoms list
    public void mcMoveAction(MCMoveEvent evt) {
        int fromIndex, toIndex;
        if(evt.acceptedMove) {
            fromIndex = index;
            toIndex = indexOld;
        }
        else {//move was rejected
            fromIndex = indexOld;
            toIndex = index;
        }
        AtomIterator iterator = evt.mcMove.affectedAtoms();
        iterator.reset();
        while(iterator.hasNext()) {
            Atom atom = iterator.next();
            //restore lists for atoms on atom's new and old lists
            //something may not work correctly if one of these atoms shows up as "atom" later in this loop
            listIterator1.setBasis(atom.atomList[index]);
            while(listIterator1.hasNext()) {
                Atom a = listIterator1.next();
                copyList(a.atomList[fromIndex], a.atomList[toIndex]);
            }
            listIterator1.setBasis(atom.atomList[indexOld]);
            while(listIterator1.hasNext()) {
                Atom a = listIterator1.next();
                copyList(a.atomList[fromIndex], a.atomList[toIndex]);
            }
            //restore atom's list
            copyList(atom.atomList[fromIndex], atom.atomList[toIndex]);
        }
    }
    
    /**
     * Replaces elements in second list with those in the first.
     */
    private void copyList(AtomList fromList, AtomList toList) {
        toList.clear();
        listIterator0.setBasis(fromList);
        toList.addAll(listIterator0);
    }
    
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
                
                //these two lines are the only ones present in the superclass version of the loop
                sum += potential.energy(pair);
                if(sum >= Double.MAX_VALUE) return;
            }//end while
        }//end of calculate
    }//end of EnergySum inner class
}
    