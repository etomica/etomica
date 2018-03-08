/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.colloid;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.potential.P2SquareWell;
import etomica.space.Space;
import etomica.util.Debug;

public class P2SquareWellMonomer extends P2SquareWell {
    public P2SquareWellMonomer(Space space, AtomLeafAgentManager bondManager) {
        super(space, 1.0, 2.0, 1.0, true);
        this.bondManager = bondManager;
    }
    
    public double energy(IAtomList pair) {
        IAtom atom0 = pair.get(0);
        IAtom atom1 = pair.get(1);

        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        if (r2 < coreDiameterSquared*(1+1e-9)) {
            IAtomList bondList = (IAtomList)bondManager.getAgent(pair.get(0));
            boolean bonded = false;
            for (int i = 0; !bonded && i<bondList.size(); i++) {
                bonded = bondList.get(i) == pair.get(1);
            }
            if (bonded) {
                return 0;
            }
        }
        return u(r2);
    }

    /**
     * Computes next time of collision of two square-well atoms, assuming free-flight kinematics.
     * Collision may occur when cores collides, or when wells first encounter each other on
     * approach, or when they edge of the wells are reached as atoms diverge.
     */
    public double collisionTime(IAtomList pair, double falseTime) {
        IAtomKinetic coord0 = (IAtomKinetic)pair.get(0);
        IAtomKinetic coord1 = (IAtomKinetic)pair.get(1);
        dv.Ev1Mv2(coord1.getVelocity(), coord0.getVelocity());
        
        dr.Ev1Mv2(coord1.getPosition(), coord0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double v2 = dv.squared();
//        System.out.println("in CT "+r2+" "+bij);
        double time = Double.POSITIVE_INFINITY;

        if(r2 < wellDiameterSquared) {  // Already inside wells

            if(bij < 0.0) {    // Check for hard-core collision

                double discr = bij*bij - v2 * ( r2 - coreDiameterSquared );
                if(discr > 0) {  // Hard cores collide next
                    if (r2 < coreDiameterSquared+1e-9) {
                        // check for bonding
                        IAtomList bondList = (IAtomList)bondManager.getAgent(pair.get(0));
                        boolean bonded = false;
                        for (int i = 0; !bonded && i<bondList.size(); i++) {
                            bonded = bondList.get(i) == pair.get(1);
                        }
                        if (bonded) {
                            double discrBonding = discr + v2 * coreDiameterSquared*(bondFac*bondFac-1);
                            if (discrBonding > 0) {
                                // bond repulsion
                                time = (-bij - Math.sqrt(discrBonding))/v2;
                            }
                            else {
                                // bonded repulsion misses, take bond stretch time
                                time = (-bij + Math.sqrt(discr))/v2;
                            }
                        }
                        else {
                            // overlapped but unbonded
                            if (ignoreOverlap) {
                                // collide soon
                                return falseTime+0.001*Math.sqrt(dr.squared())/Math.sqrt(v2);
                            }
                            throw new RuntimeException("Overlap encountered, please fix me");
                        }
                    }
                    else {
                        // not overlapped, so not bonded, 
                        time = (-bij - Math.sqrt(discr))/v2;
                    }
                }
                else {           // Moving toward each other, but wells collide next
                    discr = bij*bij - v2 * ( r2 - wellDiameterSquared );
                    time = (-bij + Math.sqrt(discr))/v2;
                }
            }
            else {           // Moving away from each other, wells collide next
                if (r2 < coreDiameterSquared+1e-9) {
                    // check for bonding
                    IAtomList bondList = (IAtomList)bondManager.getAgent(pair.get(0));
                    boolean bonded = false;
                    for (int i = 0; !bonded && i<bondList.size(); i++) {
                        bonded = bondList.get(i) == pair.get(1);
                    }
                    if (bonded) {
                        // take bond stretch time
                        double discr = bij*bij - v2 * ( r2 - coreDiameterSquared );
                        time = (-bij + Math.sqrt(discr))/v2;
                    }
                    else {
                        // well collision
                        double discr = bij*bij - v2 * ( r2 - wellDiameterSquared );
                        time = (-bij + Math.sqrt(discr))/v2;
                    }
                }
                else {
                    double discr = bij*bij - v2 * ( r2 - wellDiameterSquared );  // This is always > 0
                    time = (-bij + Math.sqrt(discr))/v2;
                }
            }
        }
        else {              // Outside wells; look for collision at well
            if(bij < 0.0) {
                double discr = bij*bij - v2 * ( r2 - wellDiameterSquared );
                if(discr > 0) {
                    time = (-bij - Math.sqrt(discr))/v2;
                }
            }
        }
        if (Debug.ON && Debug.DEBUG_NOW && ((Debug.LEVEL > 1 && Debug.allAtoms(pair)) || time < 0)) {
            System.out.println(pair+" r2 "+r2+" bij "+bij+" time "+(time+falseTime));
        }
        return time + falseTime;
    }

    public void bump(IAtomList pair, double falseTime) {
        IAtomKinetic atom0 = (IAtomKinetic) pair.get(0);
        IAtomKinetic atom1 = (IAtomKinetic) pair.get(1);
        dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());

        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        dr.PEa1Tv1(falseTime, dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double eps = 1.0e-10;
        double rm0 = atom0.getType().rm();
        double rm1 = atom1.getType().rm();
        double reduced_m = 1.0 / (rm0 + rm1);
        double nudge = 0;
        if (2 * r2 < (coreDiameterSquared + wellDiameterSquared) || lambda == 1) {   // Hard-core collision
            if (Debug.ON && !ignoreOverlap && Math.abs(r2 - coreDiameterSquared) / coreDiameterSquared > 1.e-9) {
                throw new RuntimeException("atoms " + pair + " not at the right distance " + r2 + " " + coreDiameterSquared);
            }
            // check for bonding
            IAtomList bondList = (IAtomList) bondManager.getAgent(pair.get(0));
            boolean bonded = false;
            for (int i = 0; !bonded && i < bondList.size(); i++) {
                bonded = bondList.get(i) == pair.get(1);
            }
            if (bonded) {
                // bonded
                lastCollisionVirial = 2.0 * reduced_m * bij;
            } else {
                // unbonded
                if (bij > 0) {
                    lastCollisionVirial = 0;
                    nudge = eps;
                } else {
                    lastCollisionVirial = 2.0 * reduced_m * bij;
                }
            }
            lastEnergyChange = 0.0;
        } else {    // Well collision
            if (Debug.ON && Math.abs(r2 - wellDiameterSquared) / wellDiameterSquared > 1.e-9) {
                throw new RuntimeException("atoms " + pair + " not at the right distance " + r2 + " " + wellDiameterSquared);
            }
            // ke is kinetic energy due to components of velocity
            double ke = bij * bij * reduced_m / (2.0 * r2);
            if (bij > 0.0) {         // Separating
                if (ke < epsilon) {     // Not enough kinetic energy to escape
                    lastCollisionVirial = 2.0 * reduced_m * bij;
                    nudge = -eps;
                    lastEnergyChange = 0.0;
                } else {                 // Escape
                    lastCollisionVirial = reduced_m * (bij - Math.sqrt(bij * bij - 2.0 * r2 * epsilon / reduced_m));
                    nudge = eps;
                    lastEnergyChange = epsilon;
                }
            } else if (ke > -epsilon) {   // Approach/capture
                lastCollisionVirial = reduced_m * (bij + Math.sqrt(bij * bij + 2.0 * r2 * epsilon / reduced_m));
                nudge = -eps;
                lastEnergyChange = -epsilon;
            } else {                     // Not enough kinetic energy to overcome square-shoulder
                lastCollisionVirial = 2.0 * reduced_m * bij;
                nudge = eps;
                lastEnergyChange = 0.0;
            }
        }
        lastCollisionVirialr2 = lastCollisionVirial / r2;
        dv.Ea1Tv1(lastCollisionVirialr2, dr);
        atom0.getVelocity().PEa1Tv1(rm0, dv);
        atom1.getVelocity().PEa1Tv1(-rm1, dv);
        atom0.getPosition().PEa1Tv1(-falseTime * rm0, dv);
        atom1.getPosition().PEa1Tv1(falseTime * rm1, dv);
        if (nudge != 0) {
            if (rm0 > 0) {
                atom0.getPosition().PEa1Tv1(-nudge, dr);
            }
            if (rm1 > 0) {
                atom1.getPosition().PEa1Tv1(nudge, dr);
            }
        }
    }

    public void setBondFac(double newBondFac) {
        bondFac = newBondFac;
    }
    
    public double getBondFac() {
        return bondFac;
    }

    protected final AtomLeafAgentManager bondManager;
    protected double bondFac;
}
