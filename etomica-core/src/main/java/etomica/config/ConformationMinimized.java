/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.config;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.math.function.FunctionMultiDimensionalDifferentiable;
import etomica.math.numerical.SteepestDescent;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMolecule;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * Conformation class that will find the appropriate conformation for a molecule
 * through minimization of the intramolecular energy.  The initial guess for the
 * molecule has appropriate bond lengths and random bond angles (this class
 * cannot be used with ring molecules).  The minimize method must be called
 * before any molecules are added to a box.
 */
public class ConformationMinimized implements IConformation {

    protected boolean minimizing;
    protected Vector[] coords;
    protected final List<BondInfo> bonds;
    protected final IRandom random;

    public ConformationMinimized(List<BondInfo> bonds, IRandom random) {
        bonds.sort(new Comparator<BondInfo>() {
            @Override
            public int compare(BondInfo o1, BondInfo o2) {
                int rv = Integer.compare(o1.idx1, o2.idx1);
                if (rv != 0) return rv;
                return Integer.compare(o1.idx2, o2.idx2);
            }
        });
        this.bonds = bonds;
        this.random = random;
    }

    public void minimize(Box box, PotentialCompute pc, ISpecies species, double tol) {
        minimizing = true;
        if (box.getMoleculeList().size() > 0) {
            throw new RuntimeException("Conformation must be minimized before any molecules are in the box");
        }
        coords = null;
        box.setNMolecules(species, 1);
        IMolecule m = box.getMoleculeList(species).get(0);
        IAtomList atoms = m.getChildList();
        SteepestDescent sd = new SteepestDescent(new FunctionMultiDimensionalDifferentiable() {

            @Override
            public double df(int[] d, double[] x) {
                throw new RuntimeException("This is not the method you're looking for");
            }

            protected void assignCoords(double[] x) {
                for (int i=1; i<atoms.size(); i++) {
                    Vector p = atoms.get(i).getPosition();
                    for (int j=0; j<3; j++) {
                        p.setX(j, x[3*(i-1)+j]);
                    }
                }
            }

            @Override
            public double[] gradf(double[] x) {
                assignCoords(x);
                pc.computeAll(true);
                Vector[] forces = pc.getForces();
                double[] rv = new double[x.length];
                for (int i=1; i<atoms.size(); i++) {
                    for (int j=0; j<3; j++) {
                        rv[3 * (i - 1) + j] = -forces[i].getX(j);
                    }
                }

                return rv;
            }

            @Override
            public double f(double[] x) {
                assignCoords(x);
                double u = pc.computeAll(false);
                return u;
            }

            @Override
            public int getDimension() {
                return 3*(species.getLeafAtomCount()-1);
            }
        });
        double[] guess = new double[3*(atoms.size()-1)];
        double[] xStep = new double[guess.length];
        for (int i=1; i<atoms.size(); i++) {
            Vector p = atoms.get(i).getPosition();
            for (int j=0; j<3; j++) {
                guess[3*(i-1)+j] = p.getX(j);
                xStep[3*(i-1)+j] = 0.01;
            }
        }
        sd.minimize(guess, xStep, tol, 10000);
        Vector com = CenterOfMass.position(box, m);
        coords = new Vector[atoms.size()];
        for (int i = 0; i<coords.length; i++) {
            coords[i] = box.getSpace().makeVector();
            coords[i].Ev1Mv2(atoms.get(i).getPosition(), com);
        }
        box.setNMolecules(species, 0);
        minimizing = false;
    }

    @Override
    public void initializePositions(IAtomList atomList) {
        if (coords != null) {
            for (int i = 0; i<atomList.size(); i++) {
                atomList.get(i).getPosition().E(coords[i]);
            }
            return;
        }
        if (!minimizing) {
            throw new RuntimeException("Must call minimize before adding molecules to the box");
        }

        atomList.get(0).getPosition().E(0);
        boolean[] done = new boolean[atomList.size()];
        done[0] = true;
        List<BondInfo> pendingBonds = bonds;
        Vector dr = Space.getInstance(atomList.get(0).getPosition().getD()).makeVector();
        while (!pendingBonds.isEmpty()) {
            List<BondInfo> nowBonds = pendingBonds;
            pendingBonds = new ArrayList<>();
            for (ConformationMinimized.BondInfo bi : nowBonds) {
                int idx1 = bi.idx1;
                int idx2 = bi.idx2;
                if (!done[bi.idx1]) {
                    if (done[bi.idx2]) {
                        idx1 = bi.idx2;
                        idx2 = bi.idx1;
                    } else {
                        pendingBonds.add(bi);
                        continue;
                    }
                }
                Vector p1 = atomList.get(idx1).getPosition();
                Vector p2 = atomList.get(idx2).getPosition();
                dr.setRandomSphere(random);
                dr.TE(bi.bondLength);
                p2.Ev1Pv2(p1, dr);
                done[idx2] = true;
            }
        }
    }

    public static class BondInfo {
        public final int idx1, idx2;
        public final double bondLength;

        public BondInfo(int idx1, int idx2, double bondLength) {

            this.idx1 = Math.min(idx1, idx2);
            this.idx2 = Math.max(idx1, idx2);
            this.bondLength = bondLength;
        }
    }
}
