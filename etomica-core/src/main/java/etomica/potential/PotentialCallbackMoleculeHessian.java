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
    protected Tensor[][] atomPhiTotal;
    protected Tensor[][] moleculePhiTotal;
    protected Vector[] com;
    protected final int molD;
    protected final Tensor id; // identity tensor

    public PotentialCallbackMoleculeHessian(SpeciesManager sm, Box box) {
        this.box = box;
        atomPhiTotal = new Tensor[0][0];
        moleculePhiTotal = new Tensor[0][0];
        molD = sm.getSpecies(0).getLeafAtomCount() > 1 ? 2 : 1;
        id = box.getSpace().makeTensor();
        for (int i=0; i<id.D(); i++) id.setComponent(i,i,1);
    }

    @Override
    public boolean wantsHessian() {
        return true;
    }

    public void reset() {
        int numMolecules = box.getMoleculeList().size();
        if (moleculePhiTotal.length != molD*numMolecules) {
            moleculePhiTotal = new Tensor[molD*numMolecules][molD*numMolecules];
            for (int i=0; i<molD*numMolecules; i++) {
                for (int j=0; j<molD*numMolecules; j++) {
                    moleculePhiTotal[i][j] = box.getSpace().makeTensor();
                }
            }
            com = new Vector[numMolecules];
            for (int i=0; i<numMolecules; i++) com[i] = box.getSpace().makeVector();
        }
        IMoleculeList molecules = box.getMoleculeList();
        for (int i=0; i<numMolecules; i++) {
            com[i].E(MoleculePositionCOMPBC.com(box.getBoundary(), molecules.get(i)));
        }
        int numAtoms = box.getLeafList().size();
        if (atomPhiTotal.length != numAtoms) {
            atomPhiTotal = new Tensor[numAtoms][numAtoms];
            for (int i=0; i<numAtoms; i++) {
                for (int j=0; j<numAtoms; j++) {
                    atomPhiTotal[i][j] = box.getSpace().makeTensor();
                }
            }
        }
        else {
            for (int i=0; i<numAtoms; i++) {
                for (int j=0; j<numAtoms; j++) {
                    atomPhiTotal[i][j].E(0);
                }
            }
        }
    }

    @Override
    public void pairCompute(int i, int j, Vector dr, double[] u012) {
        // atom interacting with its own image doesn't contribute
        if (i==j) return;
        // dr = rj-ri
        double r2 = dr.squared();
        Tensor t = box.getSpace().makeTensor();
        t.Ev1v2(dr, dr);
        double f1 = (u012[1] - u012[2]) / (r2*r2);
        double f2 = -u012[1]/r2;
        t.TE(f1);
        atomPhiTotal[i][j].PE(t);
        atomPhiTotal[i][j].PEa1Tt1(f2, id);
        atomPhiTotal[j][i].PE(t);
        atomPhiTotal[j][i].PEa1Tt1(f2, id);
    }

    @Override
    public void pairComputeHessian(int i, int j, Tensor phi) {
        atomPhiTotal[i][j].PE(phi);
        phi.transpose();
        atomPhiTotal[j][i].PE(phi);
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

    public Tensor[][] getMoleculePhi(Vector[] forces) {
        int numAtoms = box.getLeafList().size();
        for (int i=0; i<numAtoms; i++) {
            // compute self terms
            atomPhiTotal[i][i].E(0);
            for (int j=0; j<numAtoms; j++) {
                if (j==i) continue;
                atomPhiTotal[i][i].ME(atomPhiTotal[i][j]);
            }
        }
        for (int i=0; i<moleculePhiTotal.length; i++) {
            for (int j=0; j<moleculePhiTotal.length; j++) {
                moleculePhiTotal[i][j].E(0);
            }
        }

        Tensor Ri = box.getSpace().makeTensor();
        Tensor Rj = box.getSpace().makeTensor();
        Tensor tmpmat = box.getSpace().makeTensor();
        IMoleculeList molecules = box.getMoleculeList();
        for (int im=0; im<molecules.size(); im++) {
            IMolecule mi = molecules.get(im);
            IAtomList atomsi = mi.getChildList();
            for (int ia=0; ia<atomsi.size(); ia++) {
                IAtom iAtom = atomsi.get(ia);
                int i = iAtom.getLeafIndex();
                Vector dri = makeDr(iAtom, im);
                setR(dri, Ri);
                Ri.transpose();

                for (int jm=0; jm<molecules.size(); jm++) {
                    IMolecule mj = molecules.get(jm);
                    IAtomList atomsj = mj.getChildList();
                    for (int ja=0; ja<atomsj.size(); ja++) {
                        IAtom jAtom = atomsj.get(ja);
                        int j = jAtom.getLeafIndex();
                        // TT
                        moleculePhiTotal[molD*im][molD*jm].PE(atomPhiTotal[i][j]);
                        if (molD > 1) {
                            setR(makeDr(jAtom, jm), Rj);
                            // TR
                            tmpmat.E(atomPhiTotal[i][j]);
                            tmpmat.TE(Rj);
                            moleculePhiTotal[molD*im][molD*jm+1].PE(tmpmat);
                            // RT
                            tmpmat.E(Ri);
                            tmpmat.TE(atomPhiTotal[i][j]);
                            moleculePhiTotal[molD*im+1][molD*jm].PE(tmpmat);
                            // RR
                            tmpmat.TE(Rj);
                            moleculePhiTotal[molD*im+1][molD*jm+1].PE(tmpmat);
                        }
                    }
                }
                if (molD > 1) {
                    // intramolecular RR correction
                    double xdotf = dri.dot(forces[i]);
                    moleculePhiTotal[molD*im+1][molD*im+1].PEa1Tt1(xdotf, id);
                    // force dr F to be symmetric.... ?????
                    Tensor foo = box.getSpace().makeTensor();
                    foo.Ev1v2(dri, forces[i]);
                    foo.TE(0.5);
                    moleculePhiTotal[molD*im+1][molD*im+1].ME(foo);
                    foo.transpose();
                    moleculePhiTotal[molD*im+1][molD*im+1].ME(foo);
                }
            }
        }
        return moleculePhiTotal;
    }
}
