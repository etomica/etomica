/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.freeenergy.npath;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Created by andrew on 4/11/17.
 */
public class P1ImageHarmonic extends Potential1 implements PotentialSoft {

    protected final Vector offset;
    protected final Vector dr;
    protected final Vector[] gradient;
    protected double w;
    protected Boundary boundary;
    protected IAtomList allAtoms;
    protected boolean zeroF = true;
    protected boolean do21;
    protected int[] partners;

    public P1ImageHarmonic(Space space, Vector offset, double w, boolean do21) {
        super(space);
        this.offset = offset;
        this.w = w;
        this.do21 = do21;
        dr = space.makeVector();
        gradient = new Vector[1];
        gradient[0] = space.makeVector();
    }
    
    public void findNOffset(Box box) {
        IAtomList atoms = box.getLeafList();
        partners = new int[atoms.getAtomCount()];
        Vector p0 = atoms.getAtom(0).getPosition();
        Boundary boundary = box.getBoundary();
        int nOffset = 0;
        for (int i=1; i<atoms.getAtomCount(); i++) {
            Vector p = atoms.getAtom(i).getPosition();
            dr.Ev1Mv2(p, p0);
            dr.ME(offset);
            boundary.nearestImage(dr);
            if (dr.squared() < 0.1) {
                nOffset = i;
            }
        }
        if (nOffset == 0) throw new RuntimeException("could not find N offset");
        for (int i = 0; i < atoms.getAtomCount(); i++) {
            if (i % (nOffset * 2) >= nOffset) {
                partners[i] = i - nOffset;
            } else {
                partners[i] = i + nOffset;
            }
        }
    }
    
    public int getPartner(int idx0) {
        return partners[idx0];
    }
    
    public void setPartner(int idx0, int idx1) {
        partners[idx0] = idx1;
        partners[idx1] = idx0;
    }

    public Vector getOffset() {
        return offset;
    }

    public double getW() {
        return w;
    }
    
    public void setW(double newW) {
        w = newW;
    }

    public double getDUDW(IAtomList atoms) {
        IAtom atom0 = atoms.getAtom(0);
        int idx0 = atom0.getLeafIndex();
        if (partners[idx0] < idx0) return 0;
        IAtom atom1 = allAtoms.getAtom(partners[idx0]);
        Vector p0 = atom0.getPosition();
        Vector p1 = atom1.getPosition();
        dr.Ev1Mv2(p1,p0);
        dr.ME(offset);
        boundary.nearestImage(dr);
        // return the full contribution for this pair.  our partner will be skipped
        double r2 = dr.squared();
        if (!do21 || w==0) return r2;
        if (r2 == 0) return 0;
        double wr2 = w*r2;
        double sqrtwr2 = Math.sqrt(wr2);
        double foo = wr2+Math.sqrt(wr2);
        return wr2*r2*(2+sqrtwr2)/(2*foo*foo);
    }

    @Override
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public void setBox(Box box) {
        allAtoms = box.getLeafList();
        boundary = box.getBoundary();
    }

    @Override
    public double energy(IAtomList atoms) {
        IAtom atom0 = atoms.getAtom(0);
        int idx0 = atom0.getLeafIndex();
        IAtom atom1 = allAtoms.getAtom(partners[idx0]);
        Vector p0 = atom0.getPosition();
        Vector p1 = atom1.getPosition();
        dr.Ev1Mv2(p1,p0);
        dr.ME(offset);
        boundary.nearestImage(dr);
        // half the energy for this pair.  energy will be called again for our partner
        double wr2 = w*dr.squared();
        if (!do21) return 0.5*wr2;
        if (wr2==0) return 0;
        return 0.5/(1/wr2+1/Math.sqrt(wr2));
    }

    @Override
    public double virial(IAtomList atoms) {
        throw new RuntimeException("Implement me (please don't)");
    }

    public void setZeroForce(boolean doZeroForce) {
        this.zeroF = doZeroForce;
    }

    @Override
    public Vector[] gradient(IAtomList atoms) {
        if (zeroF) {
            gradient[0].E(0);
            return gradient;
        }

        IAtom atom0 = atoms.getAtom(0);
        int idx0 = atom0.getLeafIndex();
        IAtom atom1 = allAtoms.getAtom(partners[idx0]);
        Vector p0 = atom0.getPosition();
        Vector p1 = atom1.getPosition();
        dr.Ev1Mv2(p1,p0);
        dr.ME(offset);
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        if (!do21) {
            gradient[0].Ea1Tv1(-2*w, dr);
            return gradient;
        }
        if (r2 == 0) {
            gradient[0].E(0);
            return gradient;
        }
        double wr2 = w*r2;
        double sqrtwr2 = Math.sqrt(wr2);
        // full gradient on this atom
        gradient[0].Ea1Tv1(-w*(2+sqrtwr2)/(1+wr2+2*sqrtwr2), dr);
        return gradient;
    }

    @Override
    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        throw new RuntimeException("Implement me (just kidding.  call a different method instead)");
    }
}
