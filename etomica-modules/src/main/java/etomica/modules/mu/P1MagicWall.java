/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.mu;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.IPotential;
import etomica.potential.Potential1;
import etomica.potential.PotentialHard;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 */
 
public class P1MagicWall extends Potential1 implements PotentialHard {
    
    private static final long serialVersionUID = 1L;
    protected final PotentialMasterList potentialMaster;
    protected NeighborListManager neighborManager;
    protected final Vector dr, dv;
    protected double lastDeltaU;
    
    public P1MagicWall(Space space, PotentialMasterList potentialMaster) {
        super(space);
        this.potentialMaster = potentialMaster;
        dr = space.makeVector();
        dv = space.makeVector();
    }
    
    public double energy(IAtomList a) {
        double e = 0.0;
        // this is probably wrong, unfortunately
        return e;
    }
    
    public void setBox(Box newBox) {
        super.setBox(newBox);
        neighborManager = potentialMaster.getNeighborManager(newBox);
    }

     
    public double collisionTime(IAtomList a, double falseTime) {
        IAtomKinetic atom = (IAtomKinetic)a.get(0);
        Vector r = atom.getPosition();
        Vector v = atom.getVelocity();
        double vx = v.getX(0);
        double rx = r.getX(0) + vx * falseTime;
        double t = - rx / vx;
        if (t < 0) {
            // moving away from the wall
            t = Double.POSITIVE_INFINITY;
        }
        return t+falseTime;
    }

    public void bump(IAtomList a, double falseTime) {
        IAtomKinetic atom = (IAtomKinetic)a.get(0);
        Vector v = atom.getVelocity();
        Vector p = atom.getPosition();
        double x = p.getX(0);
        double de = getDeltaU(atom, falseTime, x<0, true);
        if (x<0) {
            // ideal gas trying to become a SQW atom
            if (de < Double.POSITIVE_INFINITY) {
                double v0 = v.getX(0);
                double m = atom.getType().getMass();
                double ke = 0.5 * m * v0*v0;
                if (ke > de) {
                    // we have enough kinetic energy to go through the wall
                    double newV = Math.sqrt(v0*v0 - 2*de/m);
                    v.setX(0, newV);
                    p.setX(0, p.getX(0) + falseTime*(v0-newV) + 1e-10);
                    lastDeltaU = de;
//                    System.out.println(atom+" IG => SQW "+de+" "+(p.getX(0)+falseTime*v.getX(0)));
                    return;
                }
            }
            else {
                // bounce
                lastDeltaU = 0;
                v.setX(0, -v.getX(0));
                double newP = atom.getPosition().getX(0) - falseTime*v.getX(0)*2.0 - 1E-10;
                atom.getPosition().setX(0,newP);
//                System.out.println(atom+" IG => IG "+(p.getX(0)+falseTime*v.getX(0)));
            }
        }
        else {
            de = -de;
            // SQW trying to become an ideal gas atom
            double v0 = v.getX(0);
            double m = atom.getType().getMass();
            double ke = 0.5 * m * v0*v0;
            if (ke > de) {
                // we have enough kinetic energy to go through the wall
                lastDeltaU = de;
                double newV = -Math.sqrt(v0*v0 - 2*de/m);
                v.setX(0, newV);
                p.setX(0, p.getX(0) + falseTime*(v0-newV) - 1e-10);
//                System.out.println(atom+" SQW => IG "+de+" "+(p.getX(0)+falseTime*newV));
            }
            else {
                // bounce
                lastDeltaU = 0;
                v.setX(0, -v.getX(0));
                double newP = atom.getPosition().getX(0) - falseTime*v.getX(0)*2.0;
                atom.getPosition().setX(0,newP + 1e-10);
//                System.out.println(atom+" SQW => SQW "+(p.getX(0)+falseTime*v.getX(0)));
            }
        }
    }

    protected double getDeltaU(IAtomKinetic atom, double falseTime, boolean isIG2SQW, boolean countHigh) {
        Vector v = atom.getVelocity();
        Vector p = atom.getPosition();
        IAtomList[] upList = neighborManager.getUpList(atom);
        IAtomList[] downList = neighborManager.getDownList(atom);
        IPotential[] potentials = potentialMaster.getRangedPotentials(atom.getType());
        double de = 0;
        for (int ip=0; ip<upList.length; ip++) {
            if (potentials[ip].nBody() != 2) continue;
            double epsilon = ((P2SquareWellOneSide)potentials[ip]).getEpsilon();
            double sigmaSq = ((P2SquareWellOneSide)potentials[ip]).getCoreDiameter();
            sigmaSq *= sigmaSq;
            double wellSigmaSq = ((P2SquareWellOneSide)potentials[ip]).getLambda();
            wellSigmaSq *= wellSigmaSq;
            wellSigmaSq *= sigmaSq;
            for (int i = 0; i<upList[ip].size(); i++) {
                IAtomKinetic atom2 = ((IAtomKinetic)upList[ip].get(i));
                Vector pos2 = atom2.getPosition();
                Vector vel2 = atom2.getVelocity();
                double x2 = pos2.getX(0) + vel2.getX(0)*falseTime;
                if (x2 < 0 == countHigh) {
                    // we're already interacting with this atom
                    continue;
                }
                dv.Ev1Mv2(v, vel2);
                
                dr.Ev1Mv2(p, pos2);
                dr.PEa1Tv1(falseTime,dv);
                boundary.nearestImage(dr);
    
                double r2 = dr.squared();
                if (r2 < sigmaSq && isIG2SQW) {
                    return Double.POSITIVE_INFINITY;
                }
                if (r2 < wellSigmaSq) {
                    de -= epsilon;
                }
            }
            for (int i = 0; i<downList[ip].size(); i++) {
                IAtomKinetic atom2 = ((IAtomKinetic)downList[ip].get(i));
                Vector pos2 = atom2.getPosition();
                Vector vel2 = atom2.getVelocity();
                double x2 = pos2.getX(0) + vel2.getX(0)*falseTime;
                if (x2 < 0 == countHigh) {
                    // we're already interacting with this atom
                    continue;
                }
                dv.Ev1Mv2(v, vel2);
                
                dr.Ev1Mv2(p, pos2);
                dr.PEa1Tv1(falseTime,dv);
                boundary.nearestImage(dr);
    
                double r2 = dr.squared();
                if (r2 < sigmaSq && isIG2SQW) {
                    return Double.POSITIVE_INFINITY;
                }
                if (r2 < wellSigmaSq) {
                    de -= epsilon;
                }
            }
        }
        return de;
    }

    public double energyChange() {
        return lastDeltaU;
    }
    
    /**
     * not yet implemented
     */
    public double lastCollisionVirial() {return Double.NaN;}
    
    /**
     * not yet implemented.
     */
    public Tensor lastCollisionVirialTensor() {return null;}
}
   
