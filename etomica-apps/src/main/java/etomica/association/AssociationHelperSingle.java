/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.AtomArrayList;
import etomica.util.Debug;

/**
 * AssociationHelperSingle is capable of populating a list of atoms in an smer 
 * within a simulation with single-bonding (one bonding site per atom).
 * AssociationHelperSingle verifies that each atom has (at most) one bond, and
 * that the reverse bonding info is consistent (bond[A] = B, bond[B] = A).
 * Each smer might be a monomer or a dimer.
 *
 * @author Andrew Schultz
 */
public class AssociationHelperSingle implements IAssociationHelper {

    protected final AssociationManager associationManager;

    public AssociationHelperSingle(AssociationManager associationManager) {
        this.associationManager = associationManager;
    }

    public boolean populateList(AtomArrayList mySmerList, IAtom atom, boolean mightBeBroken) {
        mySmerList.clear();
        mySmerList.add(atom);
        IAtomList bondList = associationManager.getAssociatedAtoms(atom);
        if (bondList.size() > 1){
            if (mightBeBroken) {
                return true;
            }
            System.out.println(atom+" has too many bonds");
            System.out.println(bondList);
            throw new RuntimeException();
        }
        if (bondList.size() == 0){
            return false;
        }
        mySmerList.add(bondList.get(0));
        IAtom bondedAtom = bondList.get(0);
        mySmerList.add(bondedAtom);
        bondList = associationManager.getAssociatedAtoms(bondedAtom);
        if (bondList.size() != 1 || bondList.get(0) != atom) {
            if (mightBeBroken && bondList.size() == 1) {
                return true;
            }
            System.out.println("invalid bonding encountered");
            if (Debug.ON) {
                System.out.println("step "+Debug.stepCount);
            }
            System.out.println(atom+" bonded to "+bondedAtom);
            System.out.println(bondedAtom+" bonded to "+bondList);
            throw new RuntimeException();
        }
        return false;
    }
}
