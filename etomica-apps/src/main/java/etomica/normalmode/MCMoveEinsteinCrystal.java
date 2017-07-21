/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * MC move whose purpose in life is to sample an  Einstein crystal.
 * Since the energy is harmonic, each configuration can be
 * independent.
 *
 * @author Andrew Schultz
 */
public class MCMoveEinsteinCrystal extends MCMoveBox {

    public MCMoveEinsteinCrystal(Space space, IRandom random) {
        super(null);
        this.random = random;
        fixedCOM = true;
        dr = space.makeVector();
    }
    
    public void setFixedCOM(boolean newFixedCOM) {
        this.fixedCOM = newFixedCOM;
    }
    
    public boolean getFixedCOM() {
        return fixedCOM;
    }

    /**
     * Sets the Einstein spring constant
     * u = alpha * r^2
     */
    public void setAlphaEin(double newAlpha) {
        alpha = newAlpha;
    }
    
    public double getAlpha() {
        return alpha;
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public void setCoordinateDefinition(CoordinateDefinition newCoordinateDefinition) {
        coordinateDefinition = newCoordinateDefinition;
    }

    public boolean doTrial() {
        IAtomList atomList = box.getLeafList();
        double einFac = Math.sqrt(temperature/(2*alpha));
        int end = atomList.getAtomCount();
        if (fixedCOM) {
            dr.E(0);
        }
        for (int i=0; i<end; i++) {
            IAtom a = atomList.getAtom(i);
            Vector p = a.getPosition();
            Vector site = coordinateDefinition.getLatticePosition(a);
            for (int k=0; k<p.getD(); k++) {
                p.setX(k, einFac * random.nextGaussian() );
            }
            if (fixedCOM) {
                dr.PE(p);
            }
            p.PE(site);
        }
        if (fixedCOM) {
            dr.TE(-1.0/end);
            for (int i=0; i<end; i++) {
                IAtom a = atomList.getAtom(i);
                Vector p = a.getPosition();
                p.PE(dr);
            }
        }
        return true;
    }

    public AtomIterator affectedAtoms() {return null;}

    public double energyChange() {return 0;}

    public double getChi(double temperature) {
        return 1;
    }

    public void acceptNotify() {}

    public void rejectNotify() {}

    protected double alpha;
    protected double temperature;
    protected final IRandom random;
    protected CoordinateDefinition coordinateDefinition;
    protected boolean fixedCOM;
    protected final Vector dr;
}
