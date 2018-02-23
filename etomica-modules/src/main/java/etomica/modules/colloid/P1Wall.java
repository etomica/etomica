/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.colloid;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.atom.AtomLeafAgentManager;
import etomica.potential.PotentialHard;
import etomica.space.Space;
import etomica.space.Tensor;

public class P1Wall implements PotentialHard {

    public P1Wall(Space space, AtomLeafAgentManager monomerMonomerBondManager) {
        lastVirialTensor = space.makeTensor();
        this.monomerMonomerBondManager = monomerMonomerBondManager;
    }

    public void setSigma(double newSigma) {
        sigma = newSigma;
    }

    public double getSigma() {
        return sigma;
    }

    public void setRange(double newRange) {
        range = newRange;
    }
    
    public void setEpsilon(double newEpsilon) {
        epsilon = newEpsilon;
    }

    public double getEpsilon() {
        return epsilon;
    }

    public double energy(IAtomList atoms) {
        Vector p = atoms.get(0).getPosition();
        double y = Math.abs(p.getX(1));
        double Ly = boundary.getBoxSize().getX(1);
        if (y > 0.5*Ly-0.5*sigma) {
            return Double.POSITIVE_INFINITY;
        }
        if (y > 0.5*Ly-range) {
            return -epsilon;
        }
        return 0;
    }

    public double getRange() {
        return range;
    }

    public int nBody() {
        return 1;
    }

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    public void bump(IAtomList a, double falseTime) {
        IAtomKinetic atom = (IAtomKinetic)a.get(0);
        Vector v = atom.getVelocity();
        double vy = v.getX(1);
        // dv = 2*NewVelocity
        double y = atom.getPosition().getX(1);
        y += vy*falseTime;
        int idx = atom.getLeafIndex();
        lastEnergyChange = 0;
        if (monomerMonomerBondManager != null) {
            IAtomList bondAtoms = (IAtomList)monomerMonomerBondManager.getAgent(atom);
            for (int i = 0; i<bondAtoms.size(); i++) {
                if (bondAtoms.get(i).getLeafIndex() < idx) {
                    IAtomKinetic bondAtom = (IAtomKinetic)bondAtoms.get(i);
                    double by = bondAtom.getPosition().getX(1);
                    by += falseTime*bondAtom.getVelocity().getX(1);
                    if (y*by < 0) {
                        // opposite sides of the box.  ignore this collision so that this atom can go back
                        lastVirial = 0;
                        double Ly = boundary.getBoxSize().getX(1);
                        if (y > 0 && Math.abs(y-(0.5*Ly-0.5*sigma)) < 1e-10) {
                            atom.getPosition().setX(1, atom.getPosition().getX(1)+1e-10);
                        }
                        else if (y < 0 && Math.abs(y+(0.5*Ly-0.5*sigma)) < 1e-10) {
                            atom.getPosition().setX(1, atom.getPosition().getX(1)-1e-10);
                        }
                        return;
                    }
                }
            }
        }
        double m = atom.getType().getMass();
        double vNew = -vy;
        double nudgeEps = 0;
        if (Math.abs(y) < 0.5*boundary.getBoxSize().getX(1)-0.51*sigma) {
            // not colliding with wall, must be well
            if (y*vy > 0) {
                // moving out, must be capture
                double ke = 0.5*vy*vy*atom.getType().getMass();
                if (ke > -epsilon) {
                    // capture
                    lastEnergyChange = -epsilon;
                    vNew = Math.sqrt(vy*vy+2*epsilon/m);
                    nudgeEps = 1e-10;
                    if (vy < 0) {
                        vNew = -vNew;
                        nudgeEps = -nudgeEps;
                    }
                }
                else {
                    nudgeEps = 1e-10;
                    if (vNew < 0) {
                        nudgeEps = -nudgeEps;
                    }
                }
            }
            else {
                // moving out, must be escape attempt
                double ke = 0.5*vy*vy*atom.getType().getMass();
                if (ke > epsilon) {
                    vNew = Math.sqrt(vy*vy-2*epsilon/m);
                    // escape
                    lastEnergyChange = epsilon;
                    nudgeEps = 1e-10;
                    if (vy < 0) {
                        vNew = -vNew;
                        nudgeEps = -nudgeEps;
                    }
                }
                else {
                    nudgeEps = 1e-10;
                    if (vNew < 0) {
                        nudgeEps = -nudgeEps;
                    }
                }
            }
        }

        double newP = atom.getPosition().getX(1) + falseTime*(vy-vNew) + nudgeEps;
        atom.getPosition().setX(1,newP);
        lastVirial = 2.0*m*(vNew-vy);
        v.setX(1,vNew);
    }

    public double collisionTime(IAtomList a, double falseTime) {
        IAtomKinetic atom = (IAtomKinetic)a.get(0);
        double y = atom.getPosition().getX(1);
        double v = atom.getVelocity().getX(1);
        if (v == 0) return Double.POSITIVE_INFINITY;
        y += falseTime*v;
        double Ly = boundary.getBoxSize().getX(1);
        double tmin = Double.POSITIVE_INFINITY;
        if (y < 0) {
            // pretend.  it works out in the end
            y = -y;
            v = -v;
        }
        if (v < 0) {
            if (y < 0.5*Ly-range) {
                // in the middle, moving toward the middle
                return Double.POSITIVE_INFINITY;
            }
            // in the well, moving toward the middle (hit the well)
            tmin = (0.5*Ly-range-y)/v;
        }
        else if (y < 0.5*Ly-range) {
            // in the middle, moving toward the wall (hit the well)
            tmin = (0.5*Ly-range-y)/v;
        }
        else if (y > 0.5*Ly-0.5*sigma) {
            // inside the wall, collide immediately
            tmin = 0.01*sigma/Math.abs(v);
        }
        else {
            // in the well, moving toward the wall
            tmin = (0.5*Ly-0.5*sigma-y)/v;
        }
        return tmin + falseTime;
    }

    public double energyChange() {
        return lastEnergyChange;
    }

    public double lastCollisionVirial() {
        // return 0 because the wall is not a molecule!
        return 0;
    }
    
    public Tensor lastCollisionVirialTensor() {
        // let's hope people only call this on purpose.  It should really be 0.
        // we're really returning the change in momentum 
        lastVirialTensor.E(0);
        lastVirialTensor.setComponent(1, 1, lastVirial);
        return lastVirialTensor;
    }

    public double lastWallVirial() {
        double area = 1.0;
        final Vector dimensions = boundary.getBoxSize();
        area *= dimensions.getX(0)*dimensions.getX(2);
        double s = lastVirial / area;
        return s;
    }


    protected Boundary boundary;
    protected final AtomLeafAgentManager monomerMonomerBondManager;
    protected double range;
    protected double sigma;
    protected double epsilon;
    protected double lastVirial;
    protected double lastEnergyChange;
    protected final Tensor lastVirialTensor;
}
