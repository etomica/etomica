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

    public P1ImageHarmonic(ISpace space, IVector offset, double w) {
        super(space);
        this.offset = offset;
        this.w = w;
        dr = space.makeVector();
        gradient = new IVectorMutable[1];
        gradient[0] = space.makeVector();
    }

    public void setW(double newW) {
        w = newW;
    }

    public double getW() {
        return w;
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
        int n = allAtoms.getAtomCount();
        IAtom atom0 = atoms.getAtom(0);
        int idx0 = atom0.getLeafIndex();
        IAtom atom1 = null;
        if (idx0 >= n/2) {
            atom1 = atom0;
            atom0 = allAtoms.getAtom(idx0-n/2);
        }
        else {
            atom1 = allAtoms.getAtom(idx0+n/2);
        }
        IVector p0 = atom0.getPosition();
        IVector p1 = atom1.getPosition();
        dr.Ev1Mv2(p1,p0);
        dr.ME(offset);
        double u = 0.5*w*dr.squared();
        return u;
    }

    @Override
    public double virial(IAtomList atoms) {
        throw new RuntimeException("Implement me (please don't)");
    }

    @Override
    public IVector[] gradient(IAtomList atoms) {
        int n = allAtoms.getAtomCount();
        IAtom atom0 = atoms.getAtom(0);
        int idx0 = atom0.getLeafIndex();
        IAtom atom1 = null;
        boolean swapped = false;
        if (idx0 >= n/2) {
            swapped = true;
            atom1 = atom0;
            atom0 = allAtoms.getAtom(idx0-n/2);
        }
        else {
            atom1 = allAtoms.getAtom(idx0+n/2);
        }
        IVector p0 = atom0.getPosition();
        IVector p1 = atom1.getPosition();
        dr.Ev1Mv2(p1,p0);
        dr.ME(offset);
        boundary.nearestImage(dr);
        dr.DE(boundary.getBoxSize());
        gradient[0].Ea1Tv1((swapped?1:-1)*2*w, dr);
        return gradient;
    }

    @Override
    public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        throw new RuntimeException("Implement me (just kidding.  call a different method instead)");
    }
}
