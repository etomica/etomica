/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.api.IRandom;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.Space;

/**
 * Extension of MCMoveClusterAtom that moves only the second atom and only
 * along the x axis.  The atom is moved in discrete steps such that its
 * position is always an integer multiple of the discrete step size.
 * 
 * @author Andrew Schultz
 */
public class MCMoveClusterAtomDiscrete extends MCMoveAtom {

    public MCMoveClusterAtomDiscrete(IRandom random, Space _space, double dr) {
        super(random, null, _space);
        setStepSize(1.2);
        setStepSizeMin(2*dr);
        this.dr = dr;
        rPow = 0;
	}

    public void setBox(Box p) {
        super.setBox(p);
        if (translationVectors == null) {
            translationVectors = new Vector[box.getLeafList().getAtomCount()-1];
            for (int i=0; i<translationVectors.length; i++) {
                translationVectors[i] = space.makeVector();
            }
        }
    }

    public void setRPow(double newRPow) {
        rPow = newRPow;
    }

    public double getRPow() {
        return rPow;
    }

    public boolean doTrial() {
        uOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        int imax = (int)Math.ceil(stepSize / dr);
        int idr = random.nextInt(2*imax) - imax;
        if (idr >= 0) idr++;
        Vector p1 = box.getLeafList().getAtom(1).getPosition();
        oldR = p1.getX(0);
        int iOldR = (int)Math.round(oldR/dr);
        newR = (iOldR + idr)*dr;
        p1.setX(0, newR);

        IAtomList leafAtoms = box.getLeafList();
        for(int i=2; i<leafAtoms.getAtomCount(); i++) {
            translationVectors[i-1].setRandomCube(random);
            translationVectors[i-1].TE(0.5*stepSize);
            leafAtoms.getAtom(i).getPosition().PE(translationVectors[i-1]);
        }

        ((BoxCluster)box).trialNotify();
        uNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }
    
    public double getA() {
        if (uOld == 0) return Double.POSITIVE_INFINITY;
        double ratio = uNew/uOld;
        // if rPow is given, use it (if r=0, we have trouble)
        if (rPow != 0) ratio *= Math.pow(Math.abs(newR/oldR), rPow);
        // we need to give r=0 double weight since we visit r=-1 and r=+1
        else if (oldR == 0) ratio *= 0.5;
        else if (newR == 0) ratio *= 2;
        return ratio;
    }

    public double getB() {
        return 0;
    }
    
    public void rejectNotify() {
        Vector p1 = box.getLeafList().getAtom(1).getPosition();
        p1.setX(0, oldR);

        IAtomList leafAtoms = box.getLeafList();
        for(int i=2; i<leafAtoms.getAtomCount(); i++) {
            leafAtoms.getAtom(i).getPosition().ME(translationVectors[i-1]);
        }

        ((BoxCluster)box).rejectNotify();
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }

    protected Vector[] translationVectors;
    protected double dr, oldR, newR, rPow;
}
