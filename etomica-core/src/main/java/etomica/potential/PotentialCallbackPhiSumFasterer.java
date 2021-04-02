/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.molecule.DipoleSourceMolecular;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePositionCOMPBC;
import etomica.potential.compute.PotentialCallback;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

public class PotentialCallbackPhiSumFasterer implements PotentialCallback {

    protected final Box box;
    protected Vector ei, ej, dr;
    protected double secondDerivativeSum= 0;
    protected DipoleSourceMolecular dipoleSource;
    protected final Tensor iT, t, ao1o1, ao1o2, ao2o2;
    protected Vector[] com;
    protected final Tensor Ri, Rj;

    public PotentialCallbackPhiSumFasterer(Box box, DipoleSourceMolecular dipoleSource) {
        this.box = box;
        this.dipoleSource = dipoleSource;
        com = new Vector[0];
        Space space = box.getSpace();
        dr = space.makeVector();
        ei = space.makeVector();
        ej = space.makeVector();
        Ri = space.makeTensor();
        Rj = space.makeTensor();
        Vector[] a = new Vector[3];
        a[0] = space.makeVector();
        a[1] = space.makeVector();
        a[2] = space.makeVector();
        iT = space.makeTensor();
        double [] xD = {1,0,0};
        double [] yD = {0,1,0};
        double [] zD = {0,0,1};
        a[0].E(xD);
        a[1].E(yD);
        a[2].E(zD);
        iT.E(a);
        t = space.makeTensor();
        ao1o1 = space.makeTensor();
        ao1o2 = space.makeTensor();
        ao2o2 = space.makeTensor();
    }

    private Vector makeDr(IAtom atom) {
        Vector ri = atom.getPosition();
        int m = atom.getParentGroup().getIndex();
        dr.Ev1Mv2(ri, com[m]);
        box.getBoundary().nearestImage(dr);
        return dr;
    }

    protected void setR(Vector dr, Tensor R) {
        // assume 3D
        R.setComponent(0, 1, +dr.getX(2));
        R.setComponent(1, 0, -dr.getX(2));
        R.setComponent(0, 2, -dr.getX(1));
        R.setComponent(2, 0, +dr.getX(1));
        R.setComponent(1, 2, +dr.getX(0));
        R.setComponent(2, 1, -dr.getX(0));
    }

    @Override
    public void pairCompute(int i, int j, Vector dr, double[] u012) {
        // atom interacting with its own image doesn't contribute
        if (i==j) return;
        // dr = rj-ri
        double r2 = dr.squared();
        t.Ev1v2(dr, dr);
        t.TE((u012[1] - u012[2]) / (r2*r2));
        t.PEa1Tt1(-u012[1]/r2, iT);

        IAtom iAtom = box.getLeafList().get(i);
        IAtom jAtom = box.getLeafList().get(i);
        setR(makeDr(iAtom), Ri);
        setR(makeDr(jAtom), Rj);

        // RT
        ao1o2.E(Ri);
        ao1o2.TE(t);
        // RR
        ao1o2.TE(Rj);

        ao1o1.E(Ri);
        ao1o1.TE(t);
        ao1o1.TE(Ri);
        ao1o1.TE(-1);

        ao2o2.E(Rj);
        ao2o2.TE(t);
        ao2o2.TE(Rj);
        ao2o2.TE(-1);

        computeSum(iAtom.getParentGroup(), jAtom.getParentGroup(), ao1o1, ao1o2, ao2o2);
    }

    @Override
    public void pairComputeGeneral(Potential2Soft pij, IAtom atom1, IAtom atom2, Vector drij, Vector fij, Vector tij, Vector tji) {

        IPotentialPair.Hessian h = pij.d2u(drij, atom1, atom2);
        computeSum(atom1.getParentGroup(), atom2.getParentGroup(), h.o1o1, h.o1o2, h.o2o2);
    }

    protected void computeSum(IMolecule mol1, IMolecule mol2, Tensor o1o1, Tensor o1o2, Tensor o2o2) {

        ei.E(dipoleSource.getDipole(mol1));
        ej.E(dipoleSource.getDipole(mol2));
        ei.normalize();
        ej.normalize();

        double traceij = o1o2.trace();
        double traceii = o1o1.trace();
        double tracejj = o2o2.trace();


        o1o2.transpose();
        o1o2.TE(-1);
        o1o1.transpose();
        o1o1.TE(-1);
        o2o2.transpose();
        o2o2.TE(-1);


        o1o2.PEa1Tt1(traceij, iT);
        o1o1.PEa1Tt1(traceii, iT);
        o2o2.PEa1Tt1(tracejj, iT);

        dr.E(ej);
        o1o2.transform(dr);
        secondDerivativeSum += 2*ei.dot(dr);//ij

        dr.E(ei);
        o1o1.transform(dr);

        secondDerivativeSum += ei.dot(dr);//ii

        dr.E(ej);
        o2o2.transform(dr);

        secondDerivativeSum += ej.dot(dr);//jj
    }

    public void zeroSum() {
        if (com.length != box.getMoleculeList().size()) {
            com = new Vector[box.getMoleculeList().size()];
            for (int i=0; i<com.length; i++) {
                com[i] = box.getSpace().makeVector();
            }
        }
        for (int i=0; i<box.getMoleculeList().size(); i++) {
            com[i].E(MoleculePositionCOMPBC.com(box.getBoundary(), box.getMoleculeList().get(i)));
        }
        secondDerivativeSum = 0.0;
    }

    /**
     * Returns the current value of the energy sum.
     */
    public double getSum(Vector[] forces) {
        IMoleculeList molecules = box.getMoleculeList();
        for (int i=0; i<molecules.size(); i++) {
            IMolecule m = molecules.get(i);
            if (m.getChildList().size() == 1) continue;
            for (int j=0; j<m.getChildList().size(); j++) {
                IAtom a = m.getChildList().get(j);
                Vector dri = makeDr(a);
                // intramolecular RR correction
                double xdotf = dri.dot(forces[a.getLeafIndex()]);
                ao1o1.E(iT);
                ao1o1.TE(xdotf);
                ao1o1.MEv1v2(dri, forces[i]);
                ao1o2.E(0);
                ao2o2.E(0);
                computeSum(m, m, ao1o1, ao1o2, ao2o2);
            }

        }
        return secondDerivativeSum;
    }


}
