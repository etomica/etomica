/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.api.IPotentialAtomic;
import etomica.api.ISpecies;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.space.Space;

/**
 * Zero-body potential implementing the long-range correction 
 * that compensates for truncation of the potential.  Normally, the
 * concrete instance of this class is defined as an inner class to
 * a PotentialTruncation class.  The PotentialTruncation class defines
 * the scheme used to truncate the potential, and its inner Potential0Lrc
 * subclass implements the formulas needed to integrate the potential
 * and its derivatives over the range affected by the truncation.<br>
 *
 * The components and procedures related to the long-range correction 
 * are as follows.  The Potential2 constructor defines a PotentialTruncation
 * that determines whether and how the pair potential may be truncated.
 * If the PotentialTruncation is not passed as a parameter to the constructor,
 * a default is used as indicated by the Default.TRUNCATE_POTENTIALS boolean
 * flag.  By default, all pair potentials are truncated using the 
 * PotentialTruncationSimple scheme; if TRUNCATE_POTENTIALS is set to
 * false, the PotentialTruncation.Null truncation is applied to all new pair potentials.  The
 * PotentialTruncation for a potential cannot be changed after the potential
 * is instantiated.<br>
 *
 * Each PotentialTruncation class defines an inner Potential0Lrc subclass that 
 * provides the appropriate methods for computing the long-range correction 
 * to the energy and its first two derivatives.  This class is instantiated
 * by the PotentialTruncation class in its constructor.  Upon its instantiation,
 * the Potential0Lrc class is added to the group of long-range correction potentials
 * that is kept by a single Potential0GroupLrc instance in the PotentialMaster.<br>
 *
 * Before the calculate method of PotentialMaster is called to compute something,
 * its set(Box) method must have been previously called, which identifies to
 * all potentials (including Potential0GroupLrc) which box is subject to the 
 * ensuing calculation.  Potential0Group ignores this notification if the
 * given box is the same as the one specified in the previous call; otherwise
 * it passes the identified box to all the set(Box) methods (inherited from Potential0)
 * of the Potential0Lrc classes it holds.<br>
 *
 * Then when the calculate(IteratorDirective, PotentialCalculation) method of
 * PotentialMaster is invoked, it passes the call on to the calculate methods of
 * its child potentials, including the Potential0GroupLrc instance if present.
 * Potential0GroupLrc checks that a box has been specified, that its
 * enableLrc flag is <b>true</b> (the default), and that the given iteratorDirective's
 * includeP0Lrc flag is also <b>true</b> (default is <b>false</b>).  If so, it 
 * calls the calculate methods of all child Potential0Lrc classes.
 *
 * The Potential0Lrc class will use the volume from the specified box and the
 * size method of the iterator of its associated potential to determine the
 * pair density.  The Potential0Lrc methods are called if the PotentialCalculation
 * implements Potential0Calculation.
 *
 * @author David Kofke
 */
public abstract class Potential0Lrc extends Potential0 implements PotentialSoft, IPotential0Lrc {

    protected final AtomType[] types;
    protected final boolean interType;
    protected final IPotentialAtomic truncatedPotential;
    protected final int[] lrcAtomsPerMolecule = new int[2];
    protected double divisor;
    protected Box box;

    public Potential0Lrc(Space space, AtomType[] types, IPotentialAtomic truncatedPotential) {
        super(space);
        this.types = types.clone();
        if(types.length != 2) {
            throw new IllegalArgumentException("LRC developed only for two-body potentials; must give two species to constructor");
        }
        interType = (types[0] != types[1]);
        this.truncatedPotential = truncatedPotential;
        divisor = 1;
    }
    
    /**
     * Returns the potential whose truncation this lrcPotential exists to correct.
     */
    public IPotentialAtomic getTruncatedPotential() {
        return truncatedPotential;
    }

    public void setBox(Box b) {
        for (int i=0; i<2; i++) {
            if (lrcAtomsPerMolecule[i] != 0) {
                continue;
            }
            IMolecule mol = types[i].getSpecies().makeMolecule();
            IAtomList childList = mol.getChildList();
            for (int j=0; j<childList.getAtomCount(); j++) {
                if (childList.getAtom(j).getType() == types[i]) {
                    lrcAtomsPerMolecule[i]++;
                }
            }
        }
        box = b;
    }

    public void setTargetAtom(IAtom targetAtom) {
        if (targetAtom == null) {
            divisor = 1;
            return;
        }
        int typeIndex = 1;
        if (types[0] == targetAtom.getType()) {
            typeIndex = 0;
        }
        else if (types[1] != targetAtom.getType()) {
            divisor = 0;
            return;
        }
        ISpecies species = types[typeIndex].getSpecies();
        divisor = box.getNMolecules(species) * lrcAtomsPerMolecule[typeIndex];
        if (!interType) {
            divisor = divisor/2.0;
        }
    }

    public void setTargetMolecule(IMolecule targetAtom) {
        if (targetAtom == null) {
            divisor = 1;
            return;
        }
        int typeIndex = 1;
        if (types[0].getSpecies() == targetAtom.getType()) {
            typeIndex = 0;
        }
        else if (types[1].getSpecies() != targetAtom.getType()) {
            divisor = 0;
            return;
        }
        ISpecies species = types[typeIndex].getSpecies();
        divisor = box.getNMolecules(species) * lrcAtomsPerMolecule[typeIndex];
        if (!interType) {
            divisor = divisor/2.0;
        }
    }

    /**
     * Returns the number of pairs formed from molecules of the current
     * species, in the given box.
     * @return int
     */
    protected long nPairs() {
        if(divisor == 0) return 0;
        long nPairs = 0;
        ISpecies species0 = types[0].getSpecies();
        long n0 = box.getNMolecules(species0)*lrcAtomsPerMolecule[0];
        if(interType) {
            ISpecies species1 = types[1].getSpecies();
            int n1 = box.getNMolecules(species1);
            nPairs = n0*n1*lrcAtomsPerMolecule[1];
        }
        else {
            nPairs = n0*(n0-1)/2L;
        }
        return Math.round(nPairs/divisor);
    }

}
