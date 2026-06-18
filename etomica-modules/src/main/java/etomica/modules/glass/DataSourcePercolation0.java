/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;
import etomica.util.random.IRandom;

import java.util.Arrays;

/**
 * This class checks for percolation of neighbors from a random set of particles.
 */
public class DataSourcePercolation0 implements IDataSource, DataSourceIndependent {

    protected final Box box;
    protected double[] immFracs;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected final DataTag tTag, tag;

    protected final AtomTestRandom atomTest;
    protected final AtomNbrClusterer clusterer;
    protected final int[][] clusterSize;
    protected final int numAtoms;
    protected final int[] clusterStack;
    protected final boolean[] isVisited;
    protected final Space space;
    protected final Vector[] r;


    public DataSourcePercolation0(Box box, IRandom random) {
        this.box = box;
        atomTest = new AtomTestRandom(random);
        clusterer = new AtomNbrClusterer(box, atomTest, true);
        numAtoms = box.getLeafList().size();
        clusterSize = new int[numAtoms][2];
        clusterStack = new int[numAtoms];
        isVisited = new boolean[numAtoms];
        space = box.getSpace();
        tag = new DataTag();
        tTag = new DataTag();
        r = space.makeVectorArray(numAtoms);
        immFracs = new double[0];
        reset();
    }

    public void setNbrMax(double nbrMax) {
        clusterer.setNbrMax(nbrMax);
    }

    public double getNbrMax() {
        return clusterer.getNbrMax();
    }

    public void setImmFracs(double[] newImmFracs) {
        this.immFracs = newImmFracs;
        data = null;
        reset();
    }

    public void reset() {
        int n = immFracs.length;
        if (data != null && n == data.getLength()) return;
        data = new DataFunction(new int[]{n});
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("immFrac", Null.DIMENSION, new int[]{n});
        tDataInfo.addTag(tTag);
        dataInfo = new DataFunction.DataInfoFunction("percProb", new CompoundDimension(new Dimension[]{Length.DIMENSION}, new double[]{2}), this);
        dataInfo.addTag(tag);

        double[] t = tData.getData();
        for (int i = 0; i < t.length; i++) {
            t[i] = immFracs[i];
        }
    }

    @Override
    public IData getData() {
        reset(); // reallocates if needed
        double[] percP = data.getData();
        for (int i = 0; i < percP.length; i++) {
            percP[i] = 0;
        }
        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < immFracs.length; i++) {
            atomTest.setFraction(immFracs[i]);
            for (int j = 0; j < numAtoms; j++) {
                isVisited[j] = false;
            }

            clusterer.findClusters();
            int[] firstAtom = clusterer.getFirstAtom();
            int[] nextAtom = clusterer.getNextAtom();
            int nClusters = 0;
            for (int j = 0; j < firstAtom.length; j++) {
                if (firstAtom[j] == -1) break;
                nClusters++;
                clusterSize[j][0] = j;
                clusterSize[j][1] = 0;
                for (int ii = firstAtom[j]; ii != -1; ii = nextAtom[ii]) {
                    clusterSize[j][1]++;
                }
            }
            Arrays.sort(clusterSize, 0, nClusters, new java.util.Comparator<int[]>() {
                public int compare(int[] a, int[] b) {
                    return Integer.compare(b[1], a[1]);
                }
            });

            outer:
            for (int j = 0; j < nClusters; j++) {
                int c = clusterSize[j][0]; //cluster No.
                if (clusterSize[j][1] == 1) break;
                int a = firstAtom[c]; // a is 1st atom in cluster
                clusterStack[0] = a; //push a
                isVisited[a] = true; // a is visited
                r[a].E(atoms.get(a).getPosition());
                Vector tmp = space.makeVector();
                int k_top = 0;
                while (k_top > -1) {//BFS
                    a = clusterStack[k_top];
                    k_top--;
                    int[] nbrs = clusterer.nbrList[a];
                    for (int m = 0; m < nbrs.length && nbrs[m] != -1; m++) {
                        int b = nbrs[m];
                        tmp.Ev1Mv2(atoms.get(b).getPosition(), atoms.get(a).getPosition());
                        box.getBoundary().nearestImage(tmp);
                        tmp.PE(r[a]);
                        if (isVisited[b] && r[b].Mv1Squared(tmp) > 1e-8) {//percolation
                            percP[i] = 1;
                            break outer;//
                        }
                        if (!isVisited[b]) {//New nbrs
                            k_top++;
                            clusterStack[k_top] = b; //push b
                            isVisited[b] = true; // b is visited
                        }
                        r[b].E(tmp);
                    }
                }
            }
        }
        return data;
    }


    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    @Override
    public DataDoubleArray getIndependentData(int i) {
        return tData;
    }

    @Override
    public DataDoubleArray.DataInfoDoubleArray getIndependentDataInfo(int i) {
        return tDataInfo;
    }

    @Override
    public int getIndependentArrayDimension() {
        return 1;
    }

    @Override
    public DataTag getIndependentTag() {
        return tTag;
    }
}







