/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.space.Vector;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

/**
 * Extension of MCMoveClusterAtom that moves only the second atom and only
 * along the x axis.  The atom is moved in discrete steps such that its
 * position is always an integer multiple of the discrete step size.
 * 
 * @author Andrew Schultz
 */
public class MCMoveClusterAtomDiscrete extends MCMoveBoxStep {

    protected double wOld, wNew;
    protected final IRandom random;
    protected Vector[] translationVectors;
    protected double dr, oldR, newR, rPow;

    public MCMoveClusterAtomDiscrete(IRandom random, Box box, double dr) {
        super();
        this.random = random;
        setBox(box);
        setStepSize(1.2);
        setStepSizeMin(2*dr);
        this.dr = dr;
        rPow = 0;
	}

    public void setBox(Box p) {
        super.setBox(p);
        if (translationVectors == null) {
            translationVectors = new Vector[box.getLeafList().size()-1];
            for (int i=0; i<translationVectors.length; i++) {
                translationVectors[i] = box.getSpace().makeVector();
            }
        }
    }

    @Override
    public AtomIterator affectedAtoms() {
        return null;
    }

    @Override
    public double energyChange() {
        return 0;
    }

    public void setRPow(double newRPow) {
        rPow = newRPow;
    }

    public double getRPow() {
        return rPow;
    }

    public boolean doTrial() {
        wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        int imax = (int)Math.ceil(stepSize / dr);
        int idr = random.nextInt(2*imax) - imax;
        if (idr >= 0) idr++;
        Vector p1 = box.getLeafList().get(1).getPosition();
        oldR = p1.getX(0);
        int iOldR = (int)Math.round(oldR/dr);
        newR = (iOldR + idr)*dr;
        p1.setX(0, newR);

        IAtomList leafAtoms = box.getLeafList();
        for(int i = 2; i<leafAtoms.size(); i++) {
            translationVectors[i-1].setRandomCube(random);
            translationVectors[i-1].TE(0.5*stepSize);
            leafAtoms.get(i).getPosition().PE(translationVectors[i-1]);
        }

        ((BoxCluster)box).trialNotify();
        wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }
    
    public double getA() {
        if (wOld == 0) return Double.POSITIVE_INFINITY;
        double ratio = wNew/wOld;
        // if rPow is given, use it (if r=0, we have trouble)
        if (rPow != 0) ratio *= Math.pow(Math.abs(newR/oldR), rPow);
        // we need to give r=0 double weight since we visit r=-1 and r=+1
        else if (oldR == 0) ratio *= 0.5;
        else if (newR == 0) ratio *= 2;
        return ratio;
    }

    public double getChi(double temperature) {
        return 1;
    }
    
    public void rejectNotify() {
        Vector p1 = box.getLeafList().get(1).getPosition();
        p1.setX(0, oldR);

        IAtomList leafAtoms = box.getLeafList();
        for(int i = 2; i<leafAtoms.size(); i++) {
            leafAtoms.get(i).getPosition().ME(translationVectors[i-1]);
        }

        ((BoxCluster)box).rejectNotify();
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }
}
