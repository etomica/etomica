/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.util.Debug;

/**
 * Scales the momenta of all the leaf atoms in a Box such that the kinetic
 * temperature matches some value.  The net momentum is also subtracted off so
 * that there is no net momentum.
 * 
 * @author Andrew Schultz
 */
public class BoxScaleMomenta implements IAction {

    public BoxScaleMomenta(Box box, Space space) {
        this.box = box;
        momentum = space.makeVector();
    }

    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public double getTemperature() {
        return temperature;
    }
    
    public Box getBox() {
        return box;
    }
    
    public void actionPerformed() {
        momentum.E(0);
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        if (nLeaf == 0) return;
        if (nLeaf > 1) {
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtom a = leafList.getAtom(iLeaf);
                double mass = a.getType().getMass();
                if (mass != Double.POSITIVE_INFINITY) {
                    momentum.PEa1Tv1(mass,((IAtomKinetic)a).getVelocity());
                }
            }
            momentum.TE(1.0/nLeaf);
            //set net momentum to 0
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
                double rm = a.getType().rm();
                if (rm != 0) {
                    a.getVelocity().PEa1Tv1(-rm,momentum);
                }
            }
            if (Debug.ON) {
                momentum.E(0);
                for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
                    double mass = a.getType().getMass();
                    if (mass != Double.POSITIVE_INFINITY) {
                        momentum.PEa1Tv1(mass,a.getVelocity());
                    }
                }
                momentum.TE(1.0/nLeaf);
                if (Math.sqrt(momentum.squared()) > 1.e-10) {
                    System.out.println("Net momentum per leaf atom is "+momentum+" but I expected it to be 0");
                }
            }
            momentum.E(0);
        }
        
        // calculate current kinetic temperature.
        for (int i = 0; i < momentum.getD(); i++) {
            // scale independently in each dimension
            double sum = 0.0;
            for (int iAtom = 0; iAtom<nLeaf; iAtom++) {
                IAtomKinetic atom = (IAtomKinetic)leafList.getAtom(iAtom);
                double mass = atom.getType().getMass();
                if(mass == Double.POSITIVE_INFINITY) continue;
                double v = atom.getVelocity().getX(i);
                sum += mass*v*v;
            }
            if (sum == 0 && temperature != 0) {
                // wonky.  possible if you try to scale up velocities after T=0.
                // but then, you called scaleMomenta, so you're probably a bad
                // person and deserve this.
                throw new RuntimeException("atoms have no velocity component in "+i+" dimension");
            }
            double s = Math.sqrt(temperature / (sum / nLeaf));
            if (s == 1) continue;
            for (int iAtom = 0; iAtom<nLeaf; iAtom++) {
                IAtomKinetic atom = (IAtomKinetic)leafList.getAtom(iAtom);
                Vector vel = atom.getVelocity();
                vel.setX(i, vel.getX(i)*s); //scale momentum
            }
        }
    }

    protected final Box box;
    protected final Vector momentum;
    protected double temperature;
}
