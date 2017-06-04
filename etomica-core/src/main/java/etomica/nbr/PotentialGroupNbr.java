/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomType;
import etomica.api.IMolecule;
import etomica.api.IPotentialAtomic;
import etomica.atom.MoleculeSetSinglet;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.atom.iterator.AtomsetIteratorDirectable;
import etomica.atom.iterator.IteratorDirective;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialGroup;

public class PotentialGroupNbr extends PotentialGroup {

    protected PotentialGroupNbr(int nBody) {
        super(nBody);
        atomSetSinglet = new MoleculeSetSinglet();
    }

    /**
     * Performs the specified calculation over the iterates given by the iterator,
     * using the directive to set up the iterators for the sub-potentials of this group.
     */
    //TODO consider what to do with sub-potentials after target atoms are reached
    public void calculateRangeIndependent(IMolecule atom, IteratorDirective.Direction direction, IAtom targetAtom, PotentialCalculation pc) {
        //loop over sub-potentials
        //TODO consider separate loops for targetable and directable
        for (PotentialLinker link=firstRangeIndependent; link!= null; link=link.next) {
            link.iterator.setTarget(targetAtom);
            // are all iterators with basis size=1 directable and all
            // iterators with basis size=2 not-directable
            if (link.iterator instanceof AtomsetIteratorDirectable) {
                ((AtomsetIteratorDirectable)link.iterator).setDirection(direction);
            }
        }
        for (PotentialLinker link=firstRangeIndependent; link!= null; link=link.next) {
            if(!link.enabled) continue;
            atomSetSinglet.atom = atom;
            AtomsetIteratorBasisDependent atomIterator = link.iterator;
            atomIterator.setBasis(atomSetSinglet);
            atomIterator.reset();
            final IPotentialAtomic potential = link.potential;
            for (IAtomList atoms = atomIterator.next(); atoms != null;
                 atoms = atomIterator.next()) {
                pc.doCalculation(atoms, potential);
            }
        }
    }
    
    protected void addPotential(IPotentialAtomic potential, AtomsetIteratorBasisDependent iterator, IAtomType[] types) {
        super.addPotential(potential, iterator, types);
        if (potential.getRange() == Double.POSITIVE_INFINITY) {
            firstRangeIndependent = new PotentialLinker(potential, iterator, types, firstRangeIndependent);
        }
    }
    
    public boolean removePotential(IPotentialAtomic potential) {
        super.removePotential(potential);
        
        PotentialLinker previous = null;
        for(PotentialLinker link=firstRangeIndependent; link!=null; link=link.next) {
            if(link.potential == potential) {
                //found it
                if(previous == null) firstRangeIndependent = link.next;  //it's the first one
                else previous.next = link.next;          //it's not the first one
                return true;
            }
            previous = link;
        }
        return false;
    }
    
    private static final long serialVersionUID = 1L;
    protected PotentialLinker firstRangeIndependent;
    protected final MoleculeSetSinglet atomSetSinglet;
}
