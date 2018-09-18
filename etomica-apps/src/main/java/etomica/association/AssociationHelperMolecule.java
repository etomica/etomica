/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.box.Box;
import etomica.models.OPLS.SpeciesAceticAcid;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * AssociationHelperBranched is capable of populating a list of molecules in an smer 
 * within a simulation with multiple association sites per molecule.  The number of
 * sites per molecule can be set, but ought to be at least 3
 * (AssociationManagerSingle or AssociationManagerDouble are more efficient for
 * single-site and double-site bonding).  With at least 3 sites, branching can
 * occur and is handled by this class.  AssociationHelperBranched verifies that
 * each atom has (at most) maxBonds bonds and that all atoms bonded to a single
 * atom have a minimum distance from each other, which ensures that none of the
 * atoms are bonded to the same site of the central atom.
 * AssociationHelperBranched does not attempt to verify that bonding information
 * is consistent.
 *
 * @author Hye Min Kim
 */
public class AssociationHelperMolecule implements IAssociationHelperMolecule {
    
    protected final AssociationManagerMolecule associationManager;
    protected int maxBonds;
    protected final Boundary boundary;
    protected double minR2;
    protected final Vector dr;

    public AssociationHelperMolecule(Space space, Box box, AssociationManagerMolecule associationManager) {
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

    public void setMinBondedMoleculeDistance(double newDistance) {
        minR2 = newDistance*newDistance;
    }

    public double getMinBondedMoleculeDistance() {
        return Math.sqrt(minR2);
    }

    public boolean populateList(MoleculeArrayList smerList, IMolecule molecule, boolean mightBeBroken){
        smerList.clear();
        if (populate(smerList, molecule, mightBeBroken)) {
            return true;
        }
        // all is well
        return false;
    }

    // return true indicates invalid bonding encountered
    // return false indicates no invalid bonding encountered
    protected boolean populate(MoleculeArrayList smerList, IMolecule molecule, boolean mightBeBroken) {
        smerList.add(molecule);
        IMoleculeList bondList = associationManager.getAssociatedMolecules(molecule);
        if (!validateBondList(molecule, bondList, mightBeBroken)) {
            return true;
        }
        for (int i = 0; i<bondList.size(); i++) {
            if (smerList.indexOf(bondList.get(i)) == -1) {
                if (populate(smerList, bondList.get(i), mightBeBroken)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Returns false if invalid bonding is recognized
     */
    protected boolean validateBondList(IMolecule molecule, IMoleculeList bondList, boolean mightBeBroken) {
        if (bondList.size() > maxBonds){
            if (mightBeBroken) {
                return false;
            }
            System.out.println(molecule+" has too many bonds");
            System.out.println(bondList);
            throw new RuntimeException();
        }
        if (minR2 > 0 && bondList.size() > 1){
            for (int i = 0; i<bondList.size()-1; i++) {
                IMolecule moleculei = bondList.get(i);
                for (int j = i+1; j<bondList.size(); j++) {
                    IMolecule moleculej = bondList.get(j);
                    dr.Ev1Mv2(moleculei.getChildList().get(SpeciesAceticAcid.indexC).getPosition(), moleculej.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
                    boundary.nearestImage(dr);
                    if (dr.squared() < minR2) {
                        if (mightBeBroken) {
                            return false;
                        }
                        System.out.println(molecule+" has multiple bonds at one site");
                        System.out.println(bondList);
                        System.out.println(moleculei+" "+moleculej+" are too close, |dr|="+Math.sqrt(dr.squared()));
                        throw new RuntimeException();
                    }
                }
            }
        }
        return true;
    }        
}
