/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.api.*;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.atom.AtomArrayList;
import etomica.space.Vector;
import etomica.space.Space;

/**
 * AssociationHelperBranched is capable of populating a list of atoms in an smer 
 * within a simulation with multiple association sites per atom.  The number of
 * sites per atom can be set, but ought to be at least 3
 * (AssociationManagerSingle or AssociationManagerDouble are more efficient for
 * single-site and double-site bonding).  With at least 3 sites, branching can
 * occur and is handled by this class.  AssociationHelperBranched verifies that
 * each atom has (at most) maxBonds bonds and that all atoms bonded to a single
 * atom have a minimum distance from each other, which ensures that none of the
 * atoms are bonded to the same site of the central atom.
 * AssociationHelperBranched does not attempt to verify that bonding information
 * is consistent.
 *
 * @author Andrew Schultz
 */
public class AssociationHelperBranched implements IAssociationHelper {
    
    protected final AssociationManager associationManager;
    protected int maxBonds;
    protected final IBoundary boundary;
    protected double minR2;
    protected final Vector dr;

    public AssociationHelperBranched(Space space, Box box, AssociationManager associationManager) {
        this.associationManager = associationManager;
        boundary = box.getBoundary();
        dr = space.makeVector();
        // default value.
        maxBonds = 3;
    }
    
    public void setMaxBonds(int newMaxBonds) {
        if (maxBonds < 3) {
            System.out.println("using AssociationHelperBranched with maxBonds<3 is overkill.");
        }
        maxBonds = newMaxBonds;
    }
    
    public int getMaxBonds() {
        return maxBonds;
    }

    public void setMinBondedAtomDistance(double newDistance) {
        minR2 = newDistance*newDistance;
    }

    public double getMinBondedAtomDistance() {
        return Math.sqrt(minR2);
    }

    public boolean populateList(AtomArrayList smerList, IAtom atom, boolean mightBeBroken){
        smerList.clear();
        if (populate(smerList, atom, mightBeBroken)) {
            return true;
        }
        // all is well
        return false;
    }

    // return true indicates invalid bonding encountered
    // return false indicates no invalid bonding encountered
    protected boolean populate(AtomArrayList smerList, IAtom atom, boolean mightBeBroken) {
        smerList.add(atom);
        IAtomList bondList = associationManager.getAssociatedAtoms(atom);
        if (!validateBondList(atom, bondList, mightBeBroken)) {
            return true;
        }
        for (int i=0; i<bondList.getAtomCount(); i++) {
            if (smerList.indexOf(bondList.getAtom(i)) == -1) {
                if (populate(smerList, bondList.getAtom(i), mightBeBroken)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Returns false if invalid bonding is recognized
     */
    protected boolean validateBondList(IAtom atom, IAtomList bondList, boolean mightBeBroken) {
        if (bondList.getAtomCount() > maxBonds){
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
                            return false;
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
