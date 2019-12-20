/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.data.*;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.math.DoubleRange;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.*;

import java.util.Arrays;

public class DataSourcePercolation implements IDataSource, ConfigurationStorage.ConfigurationStorageListener, DataSourceIndependent {

    protected final ConfigurationStorage configStorage;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction data, immFracData, immFracPercData, chi4Data;
    protected final DataFunction[] immFracDataByType;
    protected final DataFunction.DataInfoFunction[] immFracDataByTypeInfo;
    protected DataFunction.DataInfoFunction dataInfo, immFracDataInfo, immFracPercDataInfo, chi4DataInfo;
    protected double[] percP;
    protected long[] immTotal, imm2Total;
    protected long[][] immFractionByType;
    protected final int[] numAtomsByType;
    protected final DataTag tTag, tag, immFracTag, immFracPercTag, chi4Tag;
    protected final DataTag[] immFracByTypeTag;

    protected long[] nSamples;
    protected final AtomNbrClusterer clusterer;
    protected AtomTestDeviation atomTest;
    protected final int[][] clusterSize;
    protected final int numAtoms;
    protected final int[] clusterStack;
    protected final boolean[] isVisited;
    protected final Space space;
    protected final Vector[] r;
    protected int log2StepMin;
    protected final int numTypes;

    protected final HistogramNotSoSimple histogramImmPerc;

    public DataSourcePercolation(ConfigurationStorage configStorage, AtomTestDeviation atomTest, int log2StepMin) {
        this(configStorage, atomTest, log2StepMin, null);
    }

    public DataSourcePercolation(ConfigurationStorage configStorage, AtomTestDeviation atomTest, int log2StepMin, HistogramNotSoSimple sharedHistogram) {
        this.configStorage = configStorage;
        int nt = 0;
        IAtomList atoms = configStorage.box.getLeafList();
        for (IAtom a : atoms) {
            int t = a.getType().getIndex();
            if (nt <= t) nt = t + 1;
        }
        numTypes = nt;
        numAtomsByType = new int[numTypes];
        for (IAtom a : atoms) {
            int t = a.getType().getIndex();
            numAtomsByType[t]++;
        }
        immFracDataByType = new DataFunction[numTypes];
        immFracDataByTypeInfo = new DataFunction.DataInfoFunction[numTypes];
        immFracByTypeTag = new DataTag[numTypes];
        for (int i = 0; i < numTypes; i++) {
            immFracByTypeTag[i] = new DataTag();
        }

        clusterer = new AtomNbrClusterer(configStorage.getBox(), atomTest, true);
        this.atomTest = atomTest;
        this.log2StepMin = log2StepMin;
        numAtoms = configStorage.getBox().getLeafList().size();
        clusterSize = new int[numAtoms][2];
        clusterStack = new int[numAtoms];
        isVisited = new boolean[numAtoms];
        space = configStorage.box.getSpace();
        percP = new double[0];
        imm2Total = immTotal = new long[0];
        immFractionByType = new long[0][0];
        nSamples = new long[0];
        tag = new DataTag();
        immFracTag = new DataTag();
        chi4Tag = new DataTag();
        tTag = new DataTag();
        r = space.makeVectorArray(numAtoms);
        if (sharedHistogram == null) {
            immFracPercTag = new DataTag();
            histogramImmPerc = new HistogramNotSoSimple(new DoubleRange(0, 1));
            int nbins = histogramImmPerc.getNBins();
            immFracPercData = new DataFunction(new int[]{nbins}, histogramImmPerc.getHistogram());
            immFracPercDataInfo = new DataFunction.DataInfoFunction("percolation fraction", Null.DIMENSION,
                    new DataSourceIndependentSimple(histogramImmPerc.xValues(),
                            new DataDoubleArray.DataInfoDoubleArray("immobile fraction", Null.DIMENSION, new int[]{nbins})));
            immFracPercDataInfo.addTag(immFracPercTag);
        } else {
            histogramImmPerc = sharedHistogram;
            immFracPercData = null;
            immFracPercTag = null;
            immFracPercDataInfo = null;
        }
        reallocate(0);
    }

    public void setLog2StepStart(int newStart) {
        log2StepMin = newStart;
    }

    public int getLog2StepStart() {
        return log2StepMin;
    }

    public void setNbrMax(double nbrMax) {
        clusterer.setNbrMax(nbrMax);
    }

    public double getNbrMax() {
        return clusterer.getNbrMax();
    }

    public void zeroData() {
        Arrays.fill(percP, 0);
        Arrays.fill(immTotal, 0);
        Arrays.fill(imm2Total, 0);
        Arrays.fill(nSamples, 0);
        for (int i = 0; i < immFractionByType.length; i++) {
            Arrays.fill(immFractionByType[i], 0);
        }
        histogramImmPerc.reset();
    }

    public void reallocate(int n) {
        percP = Arrays.copyOf(percP, n);
        immTotal = Arrays.copyOf(immTotal, n);
        imm2Total = Arrays.copyOf(imm2Total, n);
        int oldSize = immFractionByType.length;
        immFractionByType = Arrays.copyOf(immFractionByType, n);
        for (int i = oldSize; i < n; i++) {
            immFractionByType[i] = new long[numTypes];
        }
        nSamples = Arrays.copyOf(nSamples, n);
        data = new DataFunction(new int[]{n});
        immFracData = new DataFunction(new int[]{n});
        chi4Data = new DataFunction(new int[]{n});
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        tDataInfo.addTag(tTag);
        dataInfo = new DataFunction.DataInfoFunction("percProb", new CompoundDimension(new Dimension[]{Length.DIMENSION}, new double[]{2}), this);
        dataInfo.addTag(tag);
        immFracDataInfo = new DataFunction.DataInfoFunction("immFraction", Null.DIMENSION, this);
        immFracDataInfo.addTag(immFracTag);
        chi4DataInfo = new DataFunction.DataInfoFunction("chi4*", Null.DIMENSION, this);
        chi4DataInfo.addTag(chi4Tag);
        for (int i = 0; i < numTypes; i++) {
            immFracDataByType[i] = new DataFunction(new int[]{n});
            immFracDataByTypeInfo[i] = new DataFunction.DataInfoFunction("immFraction", Null.DIMENSION, this);
            immFracDataByTypeInfo[i].addTag(immFracByTypeTag[i]);
        }

        double[] t = tData.getData();
        if (t.length > 0) {
            double[] savedTimes = configStorage.getSavedTimes();
            double dt = savedTimes[0] - savedTimes[1];
            for (int i = 0; i < t.length; i++) {
                t[i] = dt * (1L << i);
            }
        }
    }

    @Override
    public void newConfigruation() {
        long step = configStorage.getSavedSteps()[0];
        Vector[] positions = configStorage.getSavedConfig(0);
        IAtomList atoms = configStorage.getBox().getLeafList();
        for (int i = 0; i < configStorage.getLastConfigIndex(); i++) {
            int x = Math.max(i, log2StepMin);
            if (step % (1L << x) == 0) {
                if (i >= percP.length) reallocate(i + 1); // reallocates if needed
                int immCount = 0;
                for (int j = 0; j < numAtoms; j++) {
                    isVisited[j] = false;
                }

                atomTest.setConfigIndex(i + 1);
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
                        immCount++;
                        int t = atoms.get(ii).getType().getIndex();
                        immFractionByType[i][t]++;
                    }
                }
                immTotal[i] += immCount;
                imm2Total[i] += immCount * immCount;
                java.util.Arrays.sort(clusterSize, 0, nClusters, new java.util.Comparator<int[]>() {
                    public int compare(int[] a, int[] b) {
                        return Integer.compare(b[1], a[1]);
                    }
                });

                boolean percolated = false;
                outer:
                for (int j = 0; j < nClusters; j++) {
                    int c = clusterSize[j][0]; //cluster No.
                    if(clusterSize[j][1] == 1) break;
                    int a = firstAtom[c]; // a is 1st atom in cluster
                    clusterStack[0] = a; //push a
                    isVisited[a] = true; // a is visited
                    r[a].E(positions[a]);
                    Vector tmp = space.makeVector();
                    int k_top = 0;
                    while (k_top > -1) {//BFS
                        a = clusterStack[k_top];
                        k_top--;
                        int[] nbrs = clusterer.nbrList[a];
                        for (int m = 0; m < nbrs.length && nbrs[m] != -1; m++) {
                            int b = nbrs[m];
                            tmp.Ev1Mv2(positions[b], positions[a]);
                            configStorage.getBox().getBoundary().nearestImage(tmp);
                            tmp.PE(r[a]);
                            if(isVisited[b] && r[b].Mv1Squared(tmp) > 1e-8){//percolation
                                percolated = true;
                                break outer;//
                            }
                            if(!isVisited[b]){//New nbrs
                                k_top++;
                                clusterStack[k_top] = b; //push b
                                isVisited[b] = true; // b is visited
                            }
                            r[b].E(tmp);
                        }
                    }
                }//loop over clusters
                histogramImmPerc.addValue(immCount / (double) numAtoms, percolated ? 1 : 0);
                if (percolated) percP[i]++;
                nSamples[i]++;
            }//2^i check
        }// loop over i
    }

    public HistogramNotSoSimple getHistogram() {
        return histogramImmPerc;
    }

    @Override
    public IData getData() {
        if (configStorage.getLastConfigIndex() < 1) return data;
        double[] y = data.getData();
        for (int i = 0; i < percP.length; i++) {
            y[i] = percP[i] / nSamples[i];
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

    public ImmFractionSource makeImmFractionSource() {
        return new ImmFractionSource();
    }

    public ImmFractionByTypeSource makeImmFractionSource(AtomType type) {
        return new ImmFractionByTypeSource(type);
    }

    public PercolationByImmFrac makePerclationByImmFracSource() {
        if (immFracPercData == null) throw new RuntimeException("this meter is using a shared histogram");
        return new PercolationByImmFrac();
    }

    public Chi4Source makeChi4Source() {
        return new Chi4Source();
    }

    public class ImmFractionSource implements IDataSource {

        @Override
        public IData getData() {
            if (configStorage.getLastConfigIndex() < 1) return immFracData;
            double[] yImmFrac = immFracData.getData();
            for (int i = 0; i < immTotal.length; i++) {
                yImmFrac[i] = (double) immTotal[i] / nSamples[i] / numAtoms;

            }
            return immFracData;
        }

        @Override
        public DataTag getTag() {
            return immFracTag;
        }

        @Override
        public IDataInfo getDataInfo() {
            return immFracDataInfo;
        }
    }

    public class Chi4Source implements IDataSource {

        @Override
        public IData getData() {
            if (configStorage.getLastConfigIndex() < 1) return chi4Data;
            double[] chi4 = chi4Data.getData();
            double V = configStorage.getBox().getBoundary().volume();
            for (int i = 0; i < immTotal.length; i++) {
                double imm2Avg = (double) imm2Total[i] / (numAtoms * numAtoms) / nSamples[i];
                double immAvg = (double) immTotal[i] / nSamples[i] / numAtoms;
                chi4[i] = (imm2Avg - immAvg * immAvg) * V;
            }
            return chi4Data;
        }

        @Override
        public DataTag getTag() {
            return chi4Tag;
        }

        @Override
        public IDataInfo getDataInfo() {
            return chi4DataInfo;
        }
    }

    public class ImmFractionByTypeSource implements IDataSource {

        protected final int typeIndex;

        public ImmFractionByTypeSource(AtomType atomType) {
            typeIndex = atomType.getIndex();
        }

        @Override
        public IData getData() {
            DataFunction myData = immFracDataByType[typeIndex];
            if (configStorage.getLastConfigIndex() < 1) return myData;
            int na = numAtomsByType[typeIndex];
            double[] yImmFrac = myData.getData();
            for (int i = 0; i < immFractionByType.length; i++) {
                yImmFrac[i] = (double) immFractionByType[i][typeIndex] / nSamples[i] / na;

            }
            return myData;
        }

        @Override
        public DataTag getTag() {
            return immFracByTypeTag[typeIndex];
        }

        @Override
        public IDataInfo getDataInfo() {
            return immFracDataByTypeInfo[typeIndex];
        }
    }

    public class PercolationByImmFrac implements IDataSource {

        @Override
        public IData getData() {
            histogramImmPerc.getHistogram();
            return immFracPercData;
        }

        @Override
        public DataTag getTag() {
            return immFracPercTag;
        }

        @Override
        public IDataInfo getDataInfo() {
            return immFracPercDataInfo;
        }
    }

}
