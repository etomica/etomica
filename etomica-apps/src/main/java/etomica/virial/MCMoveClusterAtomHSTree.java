/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Class that samples molecule positions based on a tree of hard spheres of
 * diameter sigma.  The tree structure is chosen randomly as a Prüfer sequence.
 * <p>
 * https://en.wikipedia.org/wiki/Pr%C3%BCfer_sequence#Algorithm_to_convert_a_Pr.C3.BCfer_sequence_into_a_tree
 *
 * @author Andrew
 */
public class MCMoveClusterAtomHSTree extends MCMoveAtom {

    public MCMoveClusterAtomHSTree(IRandom random, Space _space, double sigma) {
        super(random, null, _space);
        this.sigma = sigma;
    }

    public void setBox(Box box) {
        super.setBox(box);
        int n = box.getLeafList().getAtomCount();
        degree = new int[n];
        a = new int[n-2];
        inserted = new int[n];
        bonds = new int[n*(n-1)/2][2];
    }

    public boolean doTrial() {

        IAtomList leafAtoms = box.getLeafList();
        int n = leafAtoms.getAtomCount();
        for (int i=0; i<n; i++) {
            degree[i] = 1;
        }
        for (int i=0; i<n-2; i++) {
            a[i] = random.nextInt(n);
            degree[a[i]]++;
        }

        int numBonds = 0;
        for (int i=0; i<n-2; i++) {
            int ii = a[i];
            for (int j=0; j<n; j++) {
                if (degree[j] == 1) {
                    bonds[numBonds][0] = ii;
                    bonds[numBonds][1] = j;
                    numBonds++;
                    degree[ii]--;
                    degree[j]--;
                    break;
                }
            }
        }
        int u = -1, v = -1;
        for (int i=0; i<n; i++) {
            if (degree[i] == 1) {
                if (u==-1) {
                    u = i;
                } else {
                    v = i;
                }
            }
        }
        bonds[numBonds][0] = u;
        bonds[numBonds][1] = v;
        numBonds++;

        leafAtoms.getAtom(0).getPosition().E(0);
        inserted[0] = 0;
        int numInserted = 1;
        // inserted is a list of points that have inserted, but not coordinated
        // coordinated is a list of points that have been inserted and coordinated
        int coordinatedMask = 0;
        while (numInserted > 0) {
            int nbr = inserted[numInserted-1];
            numInserted--;
            coordinatedMask |= 1<<nbr;
            for (int iBond = 0; iBond<numBonds; iBond++) {
                int[] b = bonds[iBond];
                int nbr2 = -1;
                if (b[0] == nbr) {
                    nbr2 = b[1];
                } else if (b[1] == nbr) {
                    nbr2 = b[0];
                } else {
                    continue;
                }
                if ((coordinatedMask & (1<<nbr2)) != 0) {
                    // already inserted nbr2, move along
                    continue;
                }
                // insert nbr2 around nbr
                Vector pos = leafAtoms.getAtom(nbr2).getPosition();

                pos.setRandomInSphere(random);
                pos.TE(getSigma(nbr, nbr2));
                pos.PE(leafAtoms.getAtom(nbr).getPosition());
                inserted[numInserted] = nbr2;
                numInserted++;
            }
        }

        ((BoxCluster)box).trialNotify();
        return true;
    }

    // override this to do a mixture
    protected double getSigma(int i, int j) {
        return sigma;
    }

    public double getChi(double temperature) {
        return 1;
    }

    public void rejectNotify() {
        throw new RuntimeException("nope");
    }

    public void acceptNotify() {
        ((BoxCluster)box).acceptNotify();
    }

    protected final double sigma;
    protected int[][] bonds;
    protected int[] degree, a;
    protected int[] inserted;
}
