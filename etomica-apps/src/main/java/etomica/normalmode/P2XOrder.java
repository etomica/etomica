/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.potential.Potential2;
import etomica.potential.Potential2HardSpherical;
import etomica.potential.Potential2Spherical;
import etomica.potential.PotentialHard;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * Hard potential that enforces ordering of the x-coordinates of the
 * pairs.  Returns infinite energy if the difference in atom indexes
 * and difference in x coordinates are of opposite sign; returns
 * zero otherwise.  Designed for use in 1D simulations.
 * 
 * @author David Kofke
 * @author Jhumpa Adhikari
 */
public class P2XOrder extends Potential2 implements Potential2Spherical, PotentialHard {
    
    private static final long serialVersionUID = 1L;
    protected final Vector dr;
    protected Box box;
    protected Potential2HardSpherical wrappedPotential;
    protected boolean hasPBC;
    
    public P2XOrder(Space space, Potential2HardSpherical wrappedPotential) {
        super(space);
        dr = space.makeVector();
        this.wrappedPotential = wrappedPotential;
    }

    /**
     * Interaction energy of the pair.
     * Zero if x coordinates are ordered differently from atom indexes.
     */
    public double energy(IAtomList pair) {
        IAtom atom0 = pair.get(0);
        IAtom atom1 = pair.get(1);
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        int dI = atom1.getParentGroup().getIndex() - atom0.getParentGroup().getIndex();
        // assume 1 species
        if (hasPBC && Math.abs(dI) == box.getMoleculeList().getMoleculeCount()-1) {
            if (box.getMoleculeList().getMoleculeCount() == 2) {
                // order is more complicated for 2 atoms since an atom is on both sides of the other
                if (dr.getX(0)*dI < 0 || Math.abs(dr.getX(0)) > box.getBoundary().getBoxSize().getX(0)) {
                    // out of order
                    return Double.POSITIVE_INFINITY;
                }

                double drnew = Math.abs(dr.getX(0) - dI*box.getBoundary().getBoxSize().getX(0));
                double mindr;
                if (drnew > Math.abs(dr.getX(0))) {
                    mindr = Math.abs(dr.getX(0));
                }
                else {
                    mindr = drnew;
                }
                return wrappedPotential.u(mindr*mindr);
            }
            dr.PEa1Tv1(dI > 0 ? -1 : 1, box.getBoundary().getBoxSize());
            return (dr.getX(0) * dI > 0.0) ? Double.POSITIVE_INFINITY : wrappedPotential.u(dr.squared());
        }
        else if (dI == 1 || dI == -1) {
            return (dr.getX(0) * dI < 0.0) ? Double.POSITIVE_INFINITY : wrappedPotential.u(dr.squared());
        }
        else {
            return 0;
        }
    }
    
    /**
     * Returns infinity.
     */
    public double getRange() {
        return wrappedPotential.getRange();
    }
    
    public Potential2Spherical getWrappedPotential() {
        return wrappedPotential;
    }

    public void setBox(Box newBox) {
        box = newBox;
        wrappedPotential.setBox(newBox);
        dr.E(box.getBoundary().getBoxSize());
        box.getBoundary().nearestImage(dr);
        hasPBC = dr.getX(0) < 0.5*box.getBoundary().getBoxSize().getX(0);
    }

    public void bump(IAtomList atom, double falseTime) {
        wrappedPotential.bump(atom, falseTime);
    }

    public double collisionTime(IAtomList atom, double falseTime) {
        return wrappedPotential.collisionTime(atom, falseTime);
    }

    public double energyChange() {
        return wrappedPotential.energyChange();
    }

    public double lastCollisionVirial() {
        return wrappedPotential.lastCollisionVirial();
    }

    public Tensor lastCollisionVirialTensor() {
        return wrappedPotential.lastCollisionVirialTensor();
    }
    public double u(double r2) {
        return wrappedPotential.u(r2);
    }
}
