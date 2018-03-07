/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.colloid;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.atom.AtomLeafAgentManager;
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

    public void setBondFac(double newBondFac) {
        bondFac = newBondFac;
    }
    
    public double getBondFac() {
        return bondFac;
    }

    protected final AtomLeafAgentManager bondManager;
    protected double bondFac;
}
