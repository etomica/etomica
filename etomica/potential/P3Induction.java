/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IPotentialAtomic;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.AtomTypeAgentManager;
import etomica.atom.IAtomOriented;
import etomica.space.IOrientation;
import etomica.space.ISpace;
import etomica.space3d.OrientationFull3D;

/**
 * 3-body induction potential based on form used by Oakley and Wheatley.
 * 
 * http://dx.doi.org/10.1063/1.3059008
 *
 * @author Andrew Schultz
 */
public class P3Induction implements IPotentialAtomic {

    protected final AtomTypeAgentManager paramsManager;
    protected final ISpace space;
    protected final double[] I = new double[3];
    protected final double[] alpha = new double[3];
    protected final IVectorMutable dr1, dr2;
    protected final IVectorMutable ri, rj, rk;
    protected final IVectorMutable rij, rik;
    protected final IVectorMutable or3;
    protected IBoundary boundary;
    
    public P3Induction(ISpace space, AtomTypeAgentManager paramsManager) {
        this.space = space;
        this.paramsManager = paramsManager;
        dr1 = space.makeVector();
        dr2 = space.makeVector();
        ri = space.makeVector();
        rj = space.makeVector();
        rk = space.makeVector();
        rij = space.makeVector();
        rik = space.makeVector();
        or3 = space.makeVector();
    }
    
    public double energy(IAtomList atoms) {
        double sum = 0;
        for (int i=0; i<3; i++) {
            
            int j = (i+1)%3;
            int k = (i+2)%3;
            IAtomOriented atomi = (IAtomOriented)atoms.getAtom(i);
            IOrientation ori = atomi.getOrientation();
            IAtomOriented atomj = (IAtomOriented)atoms.getAtom(j);
            IOrientation orj = atomj.getOrientation();
            IAtomOriented atomk = (IAtomOriented)atoms.getAtom(k);
            IOrientation ork = atomk.getOrientation();
            MyAgent agi = (MyAgent)paramsManager.getAgent(atoms.getAtom(i).getType());
            MyAgent agj = (MyAgent)paramsManager.getAgent(atoms.getAtom(j).getType());
            MyAgent agk = (MyAgent)paramsManager.getAgent(atoms.getAtom(k).getType());
            for (int ip=0; ip<agi.alpha.length; ip++) {
                ri.E(atomi.getPosition());
                ri.PEa1Tv1(agi.polSite[ip].getX(0), ori.getDirection());
                if (ori instanceof OrientationFull3D) {
                    IVector or2 = ((OrientationFull3D)ori).getSecondaryDirection();
                    ri.PEa1Tv1(agi.polSite[ip].getX(1), or2);
                    or3.E(ori.getDirection());
                    or3.XE(or2);
                    ri.PEa1Tv1(agi.polSite[ip].getX(2), or3);
                }
                for (int jq=0; jq<agj.q.length; jq++) {
                    rj.E(atomj.getPosition());
                    rj.PEa1Tv1(agj.qSite[jq].getX(0), orj.getDirection());
                    if (orj instanceof OrientationFull3D) {
                        IVector or2 = ((OrientationFull3D)orj).getSecondaryDirection();
                        rj.PEa1Tv1(agj.qSite[jq].getX(1), or2);
                        or3.E(ori.getDirection());
                        or3.XE(or2);
                        rj.PEa1Tv1(agj.qSite[jq].getX(2), or3);
                    }
                    rij.Ev1Mv2(rj, ri);
                    double rij2 = rij.squared();
                    // normalize rij
                    rij.TE(1/Math.sqrt(rij2));
                    for (int kq=0; kq<agk.q.length; kq++) {
                        rk.E(atomk.getPosition());
                        rk.PEa1Tv1(agk.qSite[kq].getX(0), orj.getDirection());
                        if (ork instanceof OrientationFull3D) {
                            IVector or2 = ((OrientationFull3D)ork).getSecondaryDirection();
                            rk.PEa1Tv1(agk.qSite[kq].getX(1), or2);
                            or3.E(ori.getDirection());
                            or3.XE(or2);
                            rk.PEa1Tv1(agk.qSite[kq].getX(2), or3);
                        }
                        rik.Ev1Mv2(rk, ri);
                        double rik2 = rik.squared();
                        rik.TE(1/Math.sqrt(rik2));
                        
                        sum += agi.alpha[ip]*agj.q[jq]*agk.q[kq]/(rij2*rik2)*rij.dot(rik); 
                    }
                }
            }
        }
        
        return -0.5*sum;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setBox(IBox box) {
        boundary = box.getBoundary();
    }

    public int nBody() {
        return 3;
    }

    public static class MyAgent {
        public final double[] alpha, q;
        public final IVector[] polSite, qSite;
        public MyAgent(double[] alpha, IVector[] polSite, double[] q, IVector[] qSite) {
            this.alpha = alpha;
            this.polSite = polSite;
            this.q = q;
            this.qSite = qSite;
        }
    }
}
