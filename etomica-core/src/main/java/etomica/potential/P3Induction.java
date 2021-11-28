/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomOriented;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.OrientationFull3D;

import java.util.Map;

/**
 * 3-body induction potential based on form used by Oakley and Wheatley.
 * 
 * http://dx.doi.org/10.1063/1.3059008
 *
 * @author Andrew Schultz
 */
public class P3Induction implements Potential3Soft {

    protected final Map<AtomType, MyAgent> paramsManager;
    protected final Space space;
    protected final double[] I = new double[3];
    protected final double[] alpha = new double[3];
    protected final Vector dr1, dr2;
    protected final Vector ri, rj, rk;
    protected final Vector rij, rik;
    protected final Vector or3;

    public P3Induction(Space space, Map<AtomType, MyAgent> paramsManager) {
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

    @Override
    public double u(Vector dr12, Vector dr13, Vector dr23, IAtom atom1, IAtom atom2, IAtom atom3) {
        IAtomOriented a1 = (IAtomOriented) atom1;
        IAtomOriented a2 = (IAtomOriented) atom2;
        IAtomOriented a3 = (IAtomOriented) atom3;
        return u123(a1, a2, a3) + u123(a2, a3, a1) + u123(a3, a1, a2);
    }

    public double u123(IAtomOriented atomi, IAtomOriented atomj, IAtomOriented atomk) {
        IOrientation ori = atomi.getOrientation();
        IOrientation orj = atomj.getOrientation();
        IOrientation ork = atomk.getOrientation();
        MyAgent agi = paramsManager.get(atomi.getType());
        MyAgent agj = paramsManager.get(atomj.getType());
        MyAgent agk = paramsManager.get(atomk.getType());
        double sum = 0;
        for (int ip=0; ip<agi.alpha.length; ip++) {
            ri.E(atomi.getPosition());
            ri.PEa1Tv1(agi.polSite[ip].getX(0), ori.getDirection());
            if (ori instanceof OrientationFull3D) {
                Vector or2 = ((OrientationFull3D)ori).getSecondaryDirection();
                ri.PEa1Tv1(agi.polSite[ip].getX(1), or2);
                or3.E(ori.getDirection());
                or3.XE(or2);
                ri.PEa1Tv1(agi.polSite[ip].getX(2), or3);
            }
            for (int jq=0; jq<agj.q.length; jq++) {
                rj.E(atomj.getPosition());
                rj.PEa1Tv1(agj.qSite[jq].getX(0), orj.getDirection());
                if (orj instanceof OrientationFull3D) {
                    Vector or2 = ((OrientationFull3D)orj).getSecondaryDirection();
                    rj.PEa1Tv1(agj.qSite[jq].getX(1), or2);
                    or3.E(orj.getDirection());
                    or3.XE(or2);
                    rj.PEa1Tv1(agj.qSite[jq].getX(2), or3);
                }
                rij.Ev1Mv2(rj, ri);
                double rij2 = rij.squared();
                // normalize rij
                rij.TE(1/Math.sqrt(rij2));
                for (int kq=0; kq<agk.q.length; kq++) {
                    rk.E(atomk.getPosition());
                    rk.PEa1Tv1(agk.qSite[kq].getX(0), ork.getDirection());
                    if (ork instanceof OrientationFull3D) {
                        Vector or2 = ((OrientationFull3D)ork).getSecondaryDirection();
                        rk.PEa1Tv1(agk.qSite[kq].getX(1), or2);
                        or3.E(ork.getDirection());
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

        return -sum;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public static class MyAgent {
        public final double[] alpha, q;
        public final Vector[] polSite, qSite;
        public MyAgent(double[] alpha, Vector[] polSite, double[] q, Vector[] qSite) {
            this.alpha = alpha;
            this.polSite = polSite;
            this.q = q;
            this.qSite = qSite;
        }
    }
}
