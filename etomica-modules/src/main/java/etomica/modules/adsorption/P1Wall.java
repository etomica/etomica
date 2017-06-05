/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.adsorption;

import etomica.action.AtomActionRandomizeVelocity;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.api.IRandom;
import etomica.space.Vector;
import etomica.integrator.IntegratorBox;
import etomica.potential.PotentialHard;
import etomica.space.Space;
import etomica.space.Tensor;

public class P1Wall implements PotentialHard {

    public P1Wall(Space space) {
        lastVirialTensor = space.makeTensor();
        vOld = space.makeVector();
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
        Vector p = atoms.getAtom(0).getPosition();
        double y = p.getX(1);
        double Ly = boundary.getBoxSize().getX(1);
        if (Math.abs(y) > 0.5*Ly-0.5*sigma) {
            return Double.POSITIVE_INFINITY;
        }
        if (y < -0.5*Ly+0.5*sigma+range) {
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
        IAtomKinetic atom = (IAtomKinetic)a.getAtom(0);
        Vector v = atom.getVelocity();
        double vy = v.getX(1);
        // dv = 2*NewVelocity
        double y = atom.getPosition().getX(1);
        y += vy*falseTime;
        lastEnergyChange = 0;
        double m = atom.getType().getMass();
        double vNew = -vy;
        double nudgeEps = 0;
        boolean doThermalize = false;
        if (Math.abs(y) < 0.5*boundary.getBoxSize().getX(1)-0.51*sigma) {
            // not colliding with wall, must be well
            if (y*vy > 0) {
                // moving out, must be capture
                double ke = 0.5*vy*vy*atom.getType().getMass();
                if (ke > -epsilon) {
                    // capture
                    lastEnergyChange = -epsilon;
                    vNew = -Math.sqrt(vy*vy+2*epsilon/m);
                    nudgeEps = -1e-10;
                }
                else {
                    nudgeEps = 1e-10;
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
                }
                else {
                    nudgeEps = -1e-10;
                }
            }
        }
        else {
            doThermalize = y < 0 && pThermalize > 0 && random.nextDouble() < pThermalize;
        }

        if (doThermalize) {
            vOld.E(v);
            randomizer.setTemperature(integrator.getTemperature());
            do {
                randomizer.actionPerformed(atom);
            } while (v.getX(1) < 0);
            vOld.ME(v);
            atom.getPosition().PEa1Tv1(falseTime, vOld);
        }
        else {
            double newP = atom.getPosition().getX(1) + falseTime*(vy-vNew) + nudgeEps;
            atom.getPosition().setX(1,newP);
            lastVirial = 2.0*m*(vNew-vy);
            v.setX(1,vNew);
        }
        
    }

    public double collisionTime(IAtomList a, double falseTime) {
        IAtomKinetic atom = (IAtomKinetic)a.getAtom(0);
        double y = atom.getPosition().getX(1);
        double v = atom.getVelocity().getX(1);
        if (v == 0) return Double.POSITIVE_INFINITY;
        y += falseTime*v;
        double Ly = boundary.getBoxSize().getX(1);
        double ttop = (0.5*Ly-0.5*sigma - y)/v;
        double twell = (-0.5*Ly+0.5*sigma+range - y)/v;
        double tbottom = (-0.5*Ly+0.5*sigma - y)/v;
        double tmin = Double.POSITIVE_INFINITY;
        if (ttop > 0 && ttop < tmin && v>0) tmin = ttop;
        if (twell > 0 && twell < tmin) tmin = twell;
        if (tbottom> 0 && tbottom < tmin && v<0) tmin = tbottom;
        if (tmin==Double.POSITIVE_INFINITY && v!=0) {
            // already overlapping the wall and going farther in
            // collide immediately
            tmin = 0.01*sigma/Math.abs(v);
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

    public void setThermalize(IntegratorBox integrator, double pThermalize, IRandom random) {
        this.integrator = integrator;
        this.pThermalize = pThermalize;
        this.random = random;
        randomizer = new AtomActionRandomizeVelocity(integrator.getTemperature(), random);
    }

    protected Boundary boundary;
    protected double range;
    protected double sigma;
    protected double epsilon;
    protected double lastVirial;
    protected double lastEnergyChange;
    protected final Tensor lastVirialTensor;
    protected double pThermalize;
    protected IntegratorBox integrator;
    protected IRandom random;
    protected AtomActionRandomizeVelocity randomizer;
    protected final Vector vOld;
}
