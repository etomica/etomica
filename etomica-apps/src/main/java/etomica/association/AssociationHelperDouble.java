/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.api.IBoundary;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.atom.AtomArrayList;
import etomica.space.Space;

/**
 * AssociationHelperDouble is capable of populating a list of atoms in an smer 
 * within a simulation with double-bonding (two bonding sites per atom).
 * AssociationHelperDouble verifies that each atom has (at most) two bond, and
 * that bonding information is consistent -- if B is listed as a bonding
 * partner of A, then A is listed as a bonding parter of B.
 * AssociationHelperDouble also verifies that all atoms bonded to a single atom
 * have a minimum distance from each other, which ensures that none of the
 * atoms are bonded to the same site of the central atom.
 *
 * @author Andrew Schultz
 */
public class AssociationHelperDouble implements IAssociationHelper {
    
    protected final AssociationManager associationManager;
    protected final IBoundary boundary;
    protected double minR2;
    protected final Vector dr;

    public AssociationHelperDouble(Space space, Box box, AssociationManager associationManager) {
        this.associationManager = associationManager;
        boundary = box.getBoundary();
        dr = space.makeVector();
    }

    public void setMinBondedAtomDistance(double newDistance) {
        minR2 = newDistance*newDistance;
    }

    public double getMinBondedAtomDistance() {
        return Math.sqrt(minR2);
    }

    public boolean populateList(AtomArrayList smerList, IAtom atom, boolean mightBeBroken) {
        smerList.clear();
        smerList.add(atom);
        IAtomList bondList = associationManager.getAssociatedAtoms(atom);
        if (!validateBondList(atom, bondList, mightBeBroken)) {
            return true;
        }
        if (bondList.getAtomCount() == 0){
            return false;
        }
        int rv = populate(smerList, atom, bondList.getAtom(0), mightBeBroken);
        if (rv == 1) {
            // inappropriate bonding
            return true;
        }
        if (rv == 2) {
            // ring
            return false;
        }
        if (rv == 1 && bondList.getAtomCount() == 2) {
            // need to handle the other atom bonded to "atom"
            if (populate(smerList, atom, bondList.getAtom(1), mightBeBroken) == 1) {
                // inappropriate bonding
                return true;
            }
        }
        // all is well
        return false;
    }

    // return value indicates status
    // 0: encountered only valid bonds ending in a terminal atom
    // 1: encountered invalid bonds (multiple atoms on one site)
    // 2: encountered only valid bonds forming a ring
    protected int populate(AtomArrayList smerList, IAtom firstAtom, IAtom secondAtom, boolean mightBeBroken) {
        IAtom previousAtom = firstAtom;
        IAtom nextAtom = secondAtom;
        while (true) {
            smerList.add(nextAtom);
            IAtomList bondList = associationManager.getAssociatedAtoms(nextAtom);
            if (bondList.getAtomCount() == 1) {
                if (bondList.getAtom(0) != previousAtom) {
                    System.out.println("invalid bonding");
                    System.out.println(previousAtom+" bonded to "+nextAtom);
                    System.out.println(nextAtom+" bonded to "+bondList);
                    throw new RuntimeException();
                }
                return 0;
            }
            if (bondList.getAtomCount() == 0) {
                System.out.println("invalid bonding");
                System.out.println(previousAtom+" bonded to "+nextAtom);
                System.out.println(nextAtom+" bonded to nothing");
                throw new RuntimeException();
            }                
            if (!validateBondList(nextAtom, bondList, mightBeBroken)) {
                return 1;
            }
            IAtom previousPreviousAtom = previousAtom;
            previousAtom = nextAtom;
            nextAtom = bondList.getAtom(0);
            if (nextAtom == previousPreviousAtom){
                // this is the atom we just encountered (going backwards).
                // go forwards instead
                nextAtom = bondList.getAtom(1);
            } 
            else if (bondList.getAtom(1) != previousPreviousAtom) {
                // our next atom was bonded to 2 atoms, but neither were the previous atom
                System.out.println("invalid bonding");
                System.out.println(previousPreviousAtom+" bonded to "+previousAtom);
                System.out.println(previousAtom+" bonded to "+bondList);
                throw new RuntimeException();
            }                
            if (nextAtom == firstAtom){
                // we circled back (ring).  we don't need to check atom's second
                // bond
                return 2;
            }
        }
    }

    protected boolean validateBondList(IAtom atom, IAtomList bondList, boolean mightBeBroken) {
        if (bondList.getAtomCount() > 2){
            if (mightBeBroken) {
                return false;
            }
            System.out.println(atom+" has too many bonds");
            System.out.println(bondList);
            throw new RuntimeException();
        }
        if (minR2 > 0 && bondList.getAtomCount() > 1){
            for (int i=0; i<bondList.getAtomCount()-1; i++) {
                IAtom atomi = bondList.getAtom(i);
                for (int j=i+1; j<bondList.getAtomCount(); j++) {
                    IAtom atomj = bondList.getAtom(j);
                    dr.Ev1Mv2(atomi.getPosition(), atomj.getPosition());
                    boundary.nearestImage(dr);
                    if (dr.squared() < minR2) {
                        if (mightBeBroken) {
                            return true;
                        }
                        System.out.println(atom+" has multiple bonds at one site");
                        System.out.println(bondList);
                        System.out.println(atomi+" "+atomj+" are too close, |dr|="+Math.sqrt(dr.squared()));
                        throw new RuntimeException();
                    }
                }
            }
        }
        return true;
    }        
}
