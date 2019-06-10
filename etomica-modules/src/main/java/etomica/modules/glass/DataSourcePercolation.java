/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.atom.IAtomList;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Time;

import java.util.Arrays;

public class DataSourcePercolation implements IDataSource, ConfigurationStorage.ConfigurationStorageListener, DataSourceIndependent {

    protected final ConfigurationStorage configStorage;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction data, errData;
    protected DataFunction.DataInfoFunction dataInfo;
    protected double[] percP;
    protected final DataTag tTag, tag;
    protected long[] nSamples;
    protected final AtomNbrClusterer clusterer;
    protected AtomTestDeviation atomTest;
    protected final int[][] clusterSize;
    protected final int numAtoms;
    protected final int[] clusterStack;
    protected final boolean[] isVisited;
    protected final Vector[] dr;
    protected final Space space;
    protected final Vector L;


    public DataSourcePercolation(ConfigurationStorage configStorage, AtomTestDeviation atomTest) {
        this.configStorage = configStorage;

        clusterer = new AtomNbrClusterer(configStorage.getBox(), atomTest,true);
        this.atomTest = atomTest;
        numAtoms = configStorage.getBox().getLeafList().size();
        clusterSize = new int[numAtoms][2];
        clusterStack = new int[numAtoms];
        dr = new Vector[numAtoms];
        isVisited= new boolean[numAtoms];
        space = configStorage.box.getSpace();
        L = configStorage.box.getBoundary().getBoxSize();
        percP = new double[0];
        nSamples = new long[0];
        tag = new DataTag();
        tTag = new DataTag();
        reset();
    }

    public void setNbrMax(double nbrMax) { clusterer.setNbrMax(nbrMax);}

    public double getNbrMax() {
        return clusterer.getNbrMax();
    }


    public void reset() {
        int n = configStorage.getLastConfigIndex();
        if (n + 1 == percP.length && data != null) return;
        if (n < 1) n = 0;
        percP = Arrays.copyOf(percP, n);
        nSamples = Arrays.copyOf(nSamples, n);
        data = new DataFunction(new int[]{n});
        errData = new DataFunction(new int[]{n});
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        tDataInfo.addTag(tTag);
        dataInfo = new DataFunction.DataInfoFunction("MSD", new CompoundDimension(new Dimension[]{Length.DIMENSION}, new double[]{2}), this);
        dataInfo.addTag(tag);
        double[] t = tData.getData();
        if (t.length > 0) {
            double[] savedTimes = configStorage.getSavedTimes();
            double dt = savedTimes[0] - savedTimes[1];
            for (int i = 0; i < t.length; i++) {
                t[i] = dt * (1L << i);
            }
        }
        for(int i=0; i<numAtoms; i++){
            clusterStack[i] = -1;
            dr[i] = space.makeVector();
            isVisited[i] = false;
        }
    }

    @Override
    public IData getData() {
        if (configStorage.getLastConfigIndex() < 1) return data;
        double[] y = data.getData();
        double[] yErr = errData.getData();
        for (int i = 0; i < percP.length; i++) {
            long M = nSamples[i];
            y[i] = percP[i] / M;
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
    public void newConfigruation() {
        reset(); // reallocates if needed
        long step = configStorage.getSavedSteps()[0];
        Vector[] positions = configStorage.getSavedConfig(0);
        IAtomList atoms = configStorage.getBox().getLeafList();
        for (int i = 0; i <= percP.length; i++) {
            if (step % (1L << i) == 0) {
                Vector[] iPositions = configStorage.getSavedConfig(i);
                atomTest.setConfigIndex(i);
                clusterer.findClusters();
                int[] firstAtom = clusterer.getFirstAtom();
                int[] nextAtom = clusterer.getNextAtom();
                int nClusters = 0;
                for (int j = 0; j < firstAtom.length; j++) {
                    if (firstAtom[j] == -1) break; //
                    nClusters++;
                    clusterSize[j][0] = j;
                    clusterSize[j][1] = 0;
                    for (int ii = firstAtom[j]; ii != -1; ii = nextAtom[ii]) {
                        clusterSize[j][1]++;
                    }
                }
                java.util.Arrays.sort(clusterSize, 0, nClusters, new java.util.Comparator<int[]>() {
                    public int compare(int[] a, int[] b) {
                        return Integer.compare(a[1], b[1]);
                    }
                });

                if(nClusters == 0) continue;

                for (int j = 0; j < nClusters; j++) {
                    if(clusterSize[j][1] > numAtoms/2) {
                        percP[i]++;
                        break; //do not look for further clusters; go to next time interval i.
                    }
                    int c = clusterSize[j][0];
                    int a = firstAtom[c];
                    isVisited[a] = true;
                    clusterStack[0] = a;
                    dr[0] = space.makeVector();
                    boolean test = true;
                    Vector sumDr = space.makeVector();
                    boolean isNewNbrs = true;

                    int[] nbrs0 = clusterer.nbrList[a];

                    int k = 0;
                    outer:
                    while(true){//DBS
                        int[] nbrs = clusterer.nbrList[a];
                        for(int n=0; (n<numAtoms && nbrs[n] !=-1); n++){//loop over nbrs of a
                            int b = nbrs[n];//b=6
                            if(isVisited[b] == false){// 4 is not visited
                                dr[k].Ev1Mv2(positions[a], iPositions[b]);
                                sumDr.PE(dr[k]);
                                int m = 0;
                                while(nbrs0[m] != -1){//loop over nbrs of a0
                                    if(b == nbrs0[m] && Math.abs(sumDr.getX(0))>L.getX(0)/2 && Math.abs(sumDr.getX(1))>L.getX(1)/2 && Math.abs(sumDr.getX(2))>L.getX(2)/2){
                                        percP[i]+=1;
                                        nSamples[i - 1]++;
                                        break outer;
                                    }
                                    m++;
                                }
                                isVisited[b] = true;
                                clusterStack[k] = b;
                                a = b;
                                k++;
                                isNewNbrs = true;
                                break;// go to next k
                            }
                            isNewNbrs = false;
                        }
                        if(isNewNbrs == false){//no new
                            clusterStack[k] = -1;
                            sumDr.ME(dr[k]);
                            k--;
                            a = clusterStack[k];
                            if(k == -1){
                                nSamples[i - 1]++;
                                break outer;//visited all possible paths, but no percolation detected
                            }
                        }
                    }
                }
            }
        }
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

    public interface MSDSink {
        void putMSD(int log2interval, long step, double msd);
    }
}