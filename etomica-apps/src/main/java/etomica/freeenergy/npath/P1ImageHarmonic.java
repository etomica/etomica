/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.freeenergy.npath;

import etomica.api.*;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.ISpace;
import etomica.space.Tensor;

/**
 * Created by andrew on 4/11/17.
 */
public class P1ImageHarmonic extends Potential1 implements PotentialSoft {

    protected double w;
    protected IBoundary boundary;
    protected IAtomList allAtoms;
    protected final IVector offset;
    protected final IVectorMutable dr;
    protected final IVectorMutable[] gradient;
    protected int nOffset;
    protected boolean zeroF = true;

    public P1ImageHarmonic(ISpace space, IVector offset, double w) {
        super(space);
        this.offset = offset;
        this.w = w;
        dr = space.makeVector();
        gradient = new IVectorMutable[1];
        gradient[0] = space.makeVector();
    }
    
    public void findNOffset(IBox box) {
        IAtomList atoms = box.getLeafList();
        IVector p0 = atoms.getAtom(0).getPosition();
        IBoundary boundary = box.getBoundary();
        for (int i=1; i<atoms.getAtomCount(); i++) {
            IVector p = atoms.getAtom(i).getPosition();
            dr.Ev1Mv2(p, p0);
            dr.ME(offset);
            boundary.nearestImage(dr);
            if (dr.squared() < 0.1) {
                nOffset = i;
                return;
            }
        }
        throw new RuntimeException("could not find N offset");
    }

    public int getNOffset() {
        return nOffset;
    }

    public IVector getOffset() {
        return offset;
    }

    public void setW(double newW) {
        w = newW;
    }

    public double getW() {
        return w;
    }

    public double getDUDW(IAtomList atoms) {
        IAtom atom0 = atoms.getAtom(0);
        int idx0 = atom0.getLeafIndex();
        IAtom atom1 = null;
        if (idx0%(nOffset*2) >= nOffset) return 0;
        atom1 = allAtoms.getAtom(idx0+nOffset);
        IVector p0 = atom0.getPosition();
        IVector p1 = atom1.getPosition();
        dr.Ev1Mv2(p1,p0);
        dr.ME(offset);
        boundary.nearestImage(dr);
        // return the full contribution for this pair.  our partner will be skipped
        double r2 = dr.squared();
        if (r2 == 0) return 0;
        double wr2 = w*r2;
        double sqrtwr2 = Math.sqrt(wr2);
        double foo = wr2+Math.sqrt(wr2);
        double u = wr2*r2*(2+sqrtwr2)/(2*foo*foo);
        return u;
    }

    @Override
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public void setBox(IBox box) {
        allAtoms = box.getLeafList();
        boundary = box.getBoundary();
    }

    @Override
    public double energy(IAtomList atoms) {
        IAtom atom0 = atoms.getAtom(0);
        int idx0 = atom0.getLeafIndex();
        IAtom atom1 = null;
        if (idx0%(nOffset*2) >= nOffset) {
            atom1 = atom0;
            atom0 = allAtoms.getAtom(idx0-nOffset);
        }
        else {
            atom1 = allAtoms.getAtom(idx0+nOffset);
        }
        IVector p0 = atom0.getPosition();
        IVector p1 = atom1.getPosition();
        dr.Ev1Mv2(p1,p0);
        dr.ME(offset);
        boundary.nearestImage(dr);
        // half the energy for this pair.  energy will be called again for our partner
        double wr2 = w*dr.squared();
        double u = 0.5/(1/wr2+1/Math.sqrt(wr2));
        return u;
    }

    @Override
    public double virial(IAtomList atoms) {
        throw new RuntimeException("Implement me (please don't)");
    }

    public void setZeroForce(boolean doZeroForce) {
        this.zeroF = doZeroForce;
    }

    @Override
    public IVector[] gradient(IAtomList atoms) {
        if (zeroF) {
            gradient[0].E(0);
            return gradient;
        }

        IAtom atom0 = atoms.getAtom(0);
        int idx0 = atom0.getLeafIndex();
        int iOffset = nOffset;
        if (idx0%(nOffset*2) >= nOffset) iOffset = -nOffset;
        IAtom atom1 = allAtoms.getAtom(idx0+iOffset);
        IVector p0 = atom0.getPosition();
        IVector p1 = atom1.getPosition();
        dr.Ev1Mv2(p1,p0);
        dr.ME(offset);
        boundary.nearestImage(dr);
        double r2 = dr.squared();
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
    public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        throw new RuntimeException("Implement me (just kidding.  call a different method instead)");
    }
}
