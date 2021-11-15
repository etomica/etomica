/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.IPotentialField;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

import java.util.HashMap;
import java.util.Map;

public class P1ConstraintNbr implements IPotentialField {

    protected final Vector drj, drk;
    protected int[][] neighborAtoms;
    protected double neighborRadiusSq;
    protected Boundary boundary;
    protected IAtomList leafList;
    private final Map<Box, int[][]> boxManager;
    protected int boxIndex;
    protected double constraint;

    // this could take a NeighborListManager to try to speed up finding neighbors
    public P1ConstraintNbr(Space space, double neighborDistance) {
        this(space, neighborDistance, 3.0);
    }

    public P1ConstraintNbr(Space space, double neighborDistance, double constraint) {
        boxManager = new HashMap<>();

        neighborRadiusSq = neighborDistance * neighborDistance;

        //Check for neighboring sites
        drj = space.makeVector();
        drk = space.makeVector();
        this.constraint = constraint;
    }

    /**
     * Assigns & keeps track of a set of neighbors for each box
     *
     * @param box
     */
    public void initBox(Box box) {
        //Does boxAgentManager already know what to do with this box?
        if (boxManager.containsKey(box)) {
            return;
        }

        boundary = box.getBoundary();
        IAtomList list = box.getLeafList();
        neighborAtoms = new int[list.size()][12];
        AtomArrayList tmpList = new AtomArrayList(12);

        for (int i = 0; i < list.size(); i++) {
            IAtom atomi = list.get(i);
            tmpList.clear();
            for (int j = 0; j < list.size(); j++) {
                if (i == j) continue;
                IAtom atomj = list.get(j);
                drj.Ev1Mv2(atomi.getPosition(), atomj.getPosition());
                boundary.nearestImage(drj);
                if (drj.squared() < neighborRadiusSq * 1.01) {
                    tmpList.add(atomj);
                }
            }

            for (int j = 0; j < 6; j++) {
                IAtom atomj = tmpList.get(0);
                drj.Ev1Mv2(atomi.getPosition(), atomj.getPosition());
                boundary.nearestImage(drj);
                drj.TE(1.0 / neighborRadiusSq);
                neighborAtoms[i][j * 2] = atomj.getLeafIndex();
                tmpList.remove(0);
                int indexj2 = findOppositeAtomIndex(atomi, tmpList);
                neighborAtoms[i][j * 2 + 1] = tmpList.get(indexj2).getLeafIndex();
                tmpList.remove(indexj2);
            }
        }

        boxManager.put(box, neighborAtoms);
    }

    public int[][] getNbrAtoms(Box box) {
        return boxManager.get(box);
    }

    @Override
    public double udu(IAtom atom, Vector f) {
        return u(atom);
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setBox(Box box) {
        boundary = box.getBoundary();
        leafList = box.getLeafList();
        neighborAtoms = boxManager.get(box);
        boxIndex = box.getIndex();
    }

    /**
     * find atom in list that is opposite atomi from atomj (as defined by drj)
     */
    protected int findOppositeAtomIndex(IAtom atomi, IAtomList list) {
        for (int k = 0; k < list.size(); k++) {
            IAtom atomk = list.get(k);
            drk.Ev1Mv2(atomi.getPosition(), atomk.getPosition());
            boundary.nearestImage(drk);
            double dot = drj.dot(drk);
            if (Math.abs(dot + 1.0) < 1e-10) {
                return k;
            }
        }
        throw new RuntimeException("couldn't find opposite atom");
    }

    /**
     * Returns sum of energy for all triplets containing the given atom
     */
    public double energy(IAtomList atoms) {
        IAtom atom = atoms.get(0);
        double u = energyi(atom);
        if (u == Double.POSITIVE_INFINITY) {
            return u;
        }

        int atomIndex = atom.getLeafIndex();
        int[] list = neighborAtoms[atomIndex];
        for (int i = 0; i < 12; i++) {
            u += energyij(leafList.get(list[i]), atom);
            if (u == Double.POSITIVE_INFINITY) {
                return u;
            }
        }
        return u;
    }

    /**
     * Returns the energy for atom i due to requirements for i and its
     * neighbors (but not for i as a neighbor)
     */
    public double energyi(IAtom atom) {
        Vector posAtom = atom.getPosition();

        int atomIndex = atom.getLeafIndex();
        int[] list = neighborAtoms[atomIndex];
        for (int i = 0; i < 12; i += 2) {
            IAtom atomj = leafList.get(list[i]);
            IAtom atomk = leafList.get(list[i + 1]);
            drj.Ev1Mv2(posAtom, atomj.getPosition());
            boundary.nearestImage(drj);
            // second-nearest neighbor should be at 2.66, 3rd nearest around 3.5
            // 3 seems to work OK
            if (drj.squared() > neighborRadiusSq * constraint) {
                return Double.POSITIVE_INFINITY;
            }
            drk.Ev1Mv2(posAtom, atomk.getPosition());
            boundary.nearestImage(drk);
            if (drk.squared() > neighborRadiusSq * constraint) {
                return Double.POSITIVE_INFINITY;
            }
            if (drj.dot(drk) > 0) {
                return Double.POSITIVE_INFINITY;
            }
        }
        return 0;
    }

    /**
     * Return the energy for i due to requirement that it be between j and the
     * atom opposite from j.
     */
    protected double energyij(IAtom atomi, IAtom atomj) {
        Vector posAtom = atomi.getPosition();

        int atomIndex = atomi.getLeafIndex();
        int[] list = neighborAtoms[atomIndex];
        for (int i = 0; i < 12; i += 2) {
            IAtom atomk = leafList.get(list[i]);
            if (atomk == atomj) {
                atomk = leafList.get(list[i + 1]);
            } else {
                IAtom atomkk = leafList.get(list[i + 1]);
                if (atomkk != atomj) {
                    continue;
                }
                // we found atomj, which means list[i] was atomk as we assumed!
            }
            drj.Ev1Mv2(posAtom, atomj.getPosition());
            boundary.nearestImage(drj);
            drk.Ev1Mv2(posAtom, atomk.getPosition());
            boundary.nearestImage(drk);
            if (drj.dot(drk) > 0) {
                return Double.POSITIVE_INFINITY;
            }
            // we tested the pair we want, so we're done
            return 0;
        }
        throw new RuntimeException("couldn't find " + atomj + " in " + atomi + " neighbors");
    }
}
