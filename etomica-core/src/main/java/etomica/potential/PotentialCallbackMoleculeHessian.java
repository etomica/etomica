/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePositionCOMPBC;
import etomica.potential.compute.PotentialCallback;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.species.SpeciesManager;

/**
 * Computes molecular Hessian (translation and rotation).
 */
public class PotentialCallbackMoleculeHessian implements PotentialCallback {

    protected final Box box;
    protected Vector[] com;
    protected final int molD;
    protected final Tensor id; // identity tensor
    protected final HessianConsumer callback;

    public PotentialCallbackMoleculeHessian(SpeciesManager sm, Box box, HessianConsumer callback) {
        this.box = box;
        molD = sm.getSpecies(0).getLeafAtomCount() > 1 ? 2 : 1;
        id = box.getSpace().makeTensor();
        for (int i=0; i<id.D(); i++) id.setComponent(i,i,1);
        this.callback = callback;
    }

    @Override
    public boolean wantsHessian() {
        return true;
    }

    public void reset() {
        int numMolecules = box.getMoleculeList().size();
        if (com == null || com.length != molD*numMolecules) {
            com = new Vector[numMolecules];
            for (int i=0; i<numMolecules; i++) com[i] = box.getSpace().makeVector();
        }
        IMoleculeList molecules = box.getMoleculeList();
        for (int i=0; i<numMolecules; i++) {
            com[i].E(MoleculePositionCOMPBC.com(box.getBoundary(), molecules.get(i)));
        }
    }

    @Override
    public void pairCompute(int i, int j, Vector dr, double[] u012) {
        // atom interacting with its own image doesn't contribute
        if (i == j) return;
        // dr = rj-ri
        double r2 = dr.squared();
        Tensor t = box.getSpace().makeTensor();
        t.Ev1v2(dr, dr);
        double f1 = (u012[1] - u012[2]) / (r2 * r2);
        double f2 = -u012[1] / r2;
        t.TE(f1);
        Tensor tt = box.getSpace().makeTensor();
        tt.E(t);
        tt.PEa1Tt1(f2, id);

        pairComputeHessian(i, j, tt);
//        pairComputeHessian(j, i, tt);
    }

    @Override
    public void pairComputeHessian(int i, int j, Tensor tt) {

        IAtom iAtom = box.getLeafList().get(i);
        IMolecule imol = iAtom.getParentGroup();
        IAtom jAtom = box.getLeafList().get(j);
        IMolecule jmol = jAtom.getParentGroup();
        int im = box.getMoleculeGlobalIndex(imol);
        int jm = box.getMoleculeGlobalIndex(jmol);

        if (molD > 1) {
            Tensor Ri = box.getSpace().makeTensor();
            setR(makeDr(iAtom, im), Ri);
            Ri.transpose();

            Tensor Rj = box.getSpace().makeTensor();
            setR(makeDr(jAtom, jm), Rj);
            // TR
            Tensor tr = box.getSpace().makeTensor();
            tr.E(tt);
            tr.TE(Rj);

            // RT
            Tensor rt = box.getSpace().makeTensor();
            rt.E(Ri);
            rt.TE(tt);

            // RR
            Tensor rr = box.getSpace().makeTensor();
            rr.E(rt);
            rr.TE(Rj);

            callback.takeHessian(im, jm, tt, tr, rt, rr);

            tr.transpose();
            rt.transpose();
            rr.transpose();

            callback.takeHessian(jm, im, tt, rt, tr, rr);

            tt.TE(-1);

            // RT
            rt.E(Ri);
            rt.TE(tt);

            // TR
            tr.E(tt);
            Ri.transpose();
            tr.TE(Ri); // Ri not transpose

            // RR
            rr.E(rt);
            rr.TE(Ri); // Ri not transpose

            callback.takeHessian(im, im, tt, tr, rt, rr);

            // TR
            tr.E(tt);
            tr.TE(Rj);

            // RT
            Rj.transpose();
            rt.E(Rj);  // Rj transpose
            rt.TE(tt);

            // RR
            rr.E(rt);
            Rj.transpose();
            rr.TE(Rj); // Rj not transpose

            callback.takeHessian(jm, jm, tt, tr, rt, rr);
        }

        else {
            callback.takeHessian(im, jm, tt, null, null, null);
            callback.takeHessian(jm, im, tt, null, null, null);
            tt.TE(-1);
            callback.takeHessian(im, im, tt, null, null, null);
            callback.takeHessian(jm, jm, tt, null, null, null);
        }
    }

    private Vector makeDr(IAtom atom, int m) {
        Vector ri = atom.getPosition();
        Vector dr = box.getSpace().makeVector();
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

    public void intramolecularCorrection(Vector[] forces) {
        if (molD == 1) return;
        Tensor tmpmat = box.getSpace().makeTensor();
        IMoleculeList molecules = box.getMoleculeList();
        Tensor zero = box.getSpace().makeTensor();
        for (int im=0; im<molecules.size(); im++) {
            IMolecule mi = molecules.get(im);
            IAtomList atomsi = mi.getChildList();
            for (int ia=0; ia<atomsi.size(); ia++) {
                IAtom iAtom = atomsi.get(ia);
                int i = iAtom.getLeafIndex();
                Vector dri = makeDr(iAtom, im);

                // intramolecular RR correction
                double xdotf = dri.dot(forces[i]);
                tmpmat.E(id);
                tmpmat.TE(xdotf);
                // force dr F to be symmetric.... ?????
                Tensor foo = box.getSpace().makeTensor();
                foo.Ev1v2(dri, forces[i]);
                foo.TE(0.5);
                tmpmat.ME(foo);
                foo.transpose();
                tmpmat.ME(foo);
                callback.takeHessian(im, im, zero, zero, zero, tmpmat);
            }
        }
    }

    public interface HessianConsumer {
        void takeHessian(int i, int j, Tensor tt, Tensor tr, Tensor rt, Tensor rr);
    }
}
