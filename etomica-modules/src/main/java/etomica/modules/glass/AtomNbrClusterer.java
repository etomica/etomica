/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Space;
import etomica.space.Vector;

public class AtomNbrClusterer {

    protected final Box box;
    protected final Vector dr;
    protected double nbrMax2 = 1.5 * 1.5;
    protected final AtomTestDeviation atomTestDeviation;
    protected final int[] clusters;
    protected final int[][] nbrList;
    protected final int[] nextAtom, firstAtoms, lastAtoms;
    protected final int[] skip;
    protected final boolean doNbrs;

    public AtomNbrClusterer(Box box, AtomTestDeviation atomTest){
        this(box, atomTest, false);
    }

    public AtomNbrClusterer(Box box, AtomTestDeviation atomTest, boolean doNbrs) {
        this.box = box;
        atomTestDeviation = atomTest;
        this.doNbrs = doNbrs;
        Space space = box.getSpace();
        dr = space.makeVector();
        int n = box.getLeafList().size();
        clusters = new int[n];
        firstAtoms = new int[n];
        lastAtoms = new int[n];
        nextAtom = new int[n];
        nbrList = new int[n][n];
        skip = new int[n];
    }

    public int[] getClusters() {
        return clusters;
    }

    public int[] getNextAtom() { return nextAtom; }

    public int[] getFirstAtom() {
        return firstAtoms;
    }

    public void setNbrMax(double nbrMax) {
        this.nbrMax2 = nbrMax * nbrMax;
    }

    public double getNbrMax() {
        return Math.sqrt(this.nbrMax2);
    }

    public void findClusters() {
        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < atoms.size(); i++) {
            skip[i] = nextAtom[i] = clusters[i] = -1;
            for(int j=0; j<atoms.size(); j++){
                nbrList[i][j] = -1;
            }
        }
        int nClusters = 0;

        firstAtoms[0] = -1;
        //i
        for (int i = 0; i < atoms.size(); i++) {
            IAtom a = atoms.get(i);
            if (skip[i] == -1) {
                // check if i is immobile
                skip[i] = atomTestDeviation.test(a) ? 0 : 1;
            }
            if (skip[i] > 0) continue;
            int iCluster;
            if (clusters[i] == -1) {
                // new cluster for i
                clusters[i] = nClusters;
                lastAtoms[nClusters] = firstAtoms[nClusters] = i;
                if (nClusters < firstAtoms.length - 1) firstAtoms[nClusters + 1] = -1;
                iCluster = nClusters;
                nClusters++;
            } else {
                // i in existing cluster
                iCluster = clusters[i];
            }
            Vector ri = a.getPosition();

            //j
            for (int j = i + 1; j < atoms.size(); j++) {
                if ((!doNbrs && clusters[j] == clusters[i]) || skip[j] == 1) continue;
                IAtom aj = atoms.get(j);
                if (skip[j] == -1) {
                    // check if j is immobile
                    skip[j] = atomTestDeviation.test(aj) ? 0 : 1;
                    if (skip[j] == 1) continue;
                }
                Vector rj = aj.getPosition();
                dr.Ev1Mv2(rj, ri);
                box.getBoundary().nearestImage(dr);
                if (dr.squared() > nbrMax2) continue;

                if(doNbrs) {
                    for (int k = 0; k < nbrList[i].length; k++) {
                        if (nbrList[i][k] == -1) {
                            nbrList[i][k] = j;
                            break;
                        }
                    }
                    for (int kk = 0; kk < nbrList[j].length; kk++) {
                        if (nbrList[j][kk] == -1) {
                            nbrList[j][kk] = i;
                            break;
                        }
                    }


                    if (clusters[j] == clusters[i]) continue;
                }

                if (clusters[j] == -1) {
                    // add j to iCluster
                    nextAtom[j] = firstAtoms[iCluster];
                    firstAtoms[iCluster] = j;
                    clusters[j] = iCluster;
                    continue;
                }
                int jCluster = clusters[j];
                if (iCluster < jCluster) {
                    iCluster += jCluster;
                    jCluster = iCluster - jCluster;
                    iCluster = iCluster - jCluster;
                }
                // fold iCluster into jCluster
                nextAtom[lastAtoms[jCluster]] = firstAtoms[iCluster];
                lastAtoms[jCluster] = lastAtoms[iCluster];
                for (int ii = firstAtoms[iCluster]; ii != -1; ii = nextAtom[ii]) {
                    clusters[ii] = jCluster;
                }
                lastAtoms[iCluster] = firstAtoms[iCluster] = -1;

                int lastCluster = nClusters - 1;
                if (lastCluster > iCluster) {
                    // fill in our hole
                    firstAtoms[iCluster] = firstAtoms[lastCluster];
                    lastAtoms[iCluster] = lastAtoms[lastCluster];
                    for (int ii = firstAtoms[iCluster]; ii != -1; ii = nextAtom[ii]) {
                        clusters[ii] = iCluster;
                    }
                }
                // iCluster now needs to be i's cluster
                iCluster = jCluster;
                firstAtoms[lastCluster] = lastAtoms[lastCluster] = -1;
                nClusters--;
            }
        }
    }
}





