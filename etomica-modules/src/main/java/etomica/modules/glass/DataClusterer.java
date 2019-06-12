/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSink;
import etomica.util.random.IRandom;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DataClusterer implements IDataSink {

    protected double[][] centers;
    protected final List<double[]> allData;
    protected boolean first = true;
    protected int[] clusterPop;
    protected int[] clusters;
    protected double nbrDistance;
    protected final IRandom random;
    protected double d2Max = Double.POSITIVE_INFINITY;
    protected int maxIterations = 10;
    protected int nDataP, nSamplesP;
    protected double pSum;
    protected double[] pData;

    public DataClusterer(int nClusters, IRandom random) {
        centers = new double[nClusters][0];
        allData = new ArrayList<>();
        clusterPop = new int[nClusters];
        clusters = new int[0];
        this.random = random;
        pData = new double[0];
    }

    public void setNumClusters(int nClusters) {
        if (nClusters == centers.length) return;
        if (nClusters > centers.length) {
            int oldN = centers.length;
            centers = Arrays.copyOf(centers, nClusters);
            clusterPop = new int[nClusters];
            for (int i = oldN; i < nClusters; i++) {
                centers[i] = new double[centers[0].length];
                int j = random.nextInt(allData.size());
                System.arraycopy(allData.get(j), 0, centers[i], 0, centers[i].length);
            }
        } else {
            centers = Arrays.copyOf(centers, nClusters);
            clusterPop = Arrays.copyOf(clusterPop, nClusters);
        }
    }

    public int getNumClusters() {
        return centers.length;
    }

    @Override
    public void putData(IData data) {
        if (data.getLength() == 0) {
            reset();
            return;
        }
        double[] x = new double[data.getLength()];
        for (int i = 0; i < x.length; i++) {
            x[i] = data.getValue(i);
        }

        if (nSamplesP > 0 && allData.size() != nDataP) {
            allData.clear();
            nDataP = 0;
        }
        allData.add(x);
        if (nSamplesP > 0) {
            if (nDataP == pData.length) {
                pData = Arrays.copyOf(pData, 2 * pData.length + 100);
            }
            pData[nDataP] = pSum / nSamplesP;
            pSum = nSamplesP = 0;
            nDataP++;
        } else {
            if (nDataP == pData.length) {
                pData = Arrays.copyOf(pData, 2 * pData.length + 100);
            }
            pData[nDataP] = Double.NaN;
            pSum = nSamplesP = 0;
            nDataP++;
        }
    }

    public void reset() {
        allData.clear();
        nDataP = nSamplesP = 0;
        pSum = 0;
        if (pData.length > 0) pData = new double[0];
        first = true;
    }

    public void setClusterNeighborDistance(double d) {
        nbrDistance = d / 10;
    }

    public double getClusterNeighborDistance() {
        return nbrDistance * 10;
    }

    public void setMaxClusterDistance(double dMax) {
        d2Max = dMax * dMax;
    }

    public double getMaxClusterDistance() {
        return Math.sqrt(d2Max);
    }

    public int getMaxIterations() {
        return maxIterations;
    }

    public void setMaxIterations(int maxIterations) {
        if (maxIterations < 1 || maxIterations > 100000)
            throw new RuntimeException("max iterations must be positive, no more than 100000");
        this.maxIterations = maxIterations;
    }

    public void findClusters() {
        if (centers.length * 3 > allData.size()) return;
        if (allData.size() > clusters.length) {
            int old = clusters.length;
            clusters = Arrays.copyOf(clusters, allData.size());
            for (int j = old; j < clusters.length; j++) clusters[j] = -1;
        }
        if (first) {
            int jump = allData.size() / centers.length;
            for (int i = 0; i < centers.length; i++) {
                System.arraycopy(allData.get(i * jump), 0, centers[i], 0, centers[i].length);
            }
            for (int j = 0; j < clusters.length; j++) clusters[j] = -1;
        }
        first = false;
        int dim = centers[0].length;
        int[][] iBonds = null;
        int[][] iBondCount = null;
        int totalBondCount = 0;
        double[] clusterP = new double[centers.length];
        for (int outer = 0; outer < maxIterations; outer++) {
            int nChanged = 0;
            for (int i = 0; i < centers.length; i++) clusterPop[i] = 0;
            iBonds = new int[centers.length][0];
            iBondCount = new int[centers.length][0];
            totalBondCount = 0;
            int lastCluster = -1;
            int dataOutside = 0;
            for (int j = 0; j < allData.size(); j++) {
                double[] x = allData.get(j);
                double d2Min = Double.POSITIVE_INFINITY;
                int bestCluster = -1;
                for (int i = 0; i < centers.length; i++) {
                    double d2 = 0;
                    for (int k = 0; k < dim; k++) {
                        double dk = x[k] - centers[i][k];
                        d2 += dk * dk;
                    }
                    if (d2 < d2Min) {
                        d2Min = d2;
                        bestCluster = i;
                    }
                }
                if (clusters[j] != bestCluster) nChanged++;
                clusters[j] = bestCluster;
                if (bestCluster == -1) {
                    dataOutside++;
                    continue;
                }
                clusterPop[bestCluster]++;
                if (nChanged > 0 && outer < maxIterations - 1) continue;
                if (lastCluster > -1 && lastCluster != bestCluster) {
                    int c1 = lastCluster;
                    int c2 = bestCluster;
                    if (c2 < c1) {
                        c1 = bestCluster;
                        c2 = lastCluster;
                    }
                    int check = Arrays.binarySearch(iBonds[c1], c2);
                    if (check < 0) {
                        iBonds[c1] = Arrays.copyOf(iBonds[c1], iBonds[c1].length + 1);
                        iBondCount[c1] = Arrays.copyOf(iBondCount[c1], iBonds[c1].length);
                        int insert = -(check + 1);
                        for (int i = iBonds[c1].length - 1; i > insert; i--) {
                            iBonds[c1][i] = iBonds[c1][i - 1];
                            iBondCount[c1][i] = iBondCount[c1][i - 1];
                        }
                        iBondCount[c1][insert] = 0;
                        iBonds[c1][insert] = c2;
                        check = insert;
                    }
                    iBondCount[c1][check]++;
                    totalBondCount++;
                }
                lastCluster = bestCluster;
            }
            if (dataOutside > 0) {
                System.out.println(dataOutside + " data points outside nbr distance from center");
            }
            for (int i = 0; i < centers.length; i++) {
                for (int k = 0; k < dim; k++) centers[i][k] = 0;
                clusterP[i] = 0;
            }

            for (int j = 0; j < allData.size(); j++) {
                double[] x = allData.get(j);
                int i = clusters[j];
                if (i == -1) continue;
                for (int k = 0; k < dim; k++) centers[i][k] += x[k];
                clusterP[i] += pData[j];
            }
            int nEmpty = 0;
            for (int i = 0; i < centers.length; i++) {
                if (clusterPop[i] > 0) {
                    for (int k = 0; k < dim; k++) centers[i][k] /= clusterPop[i];
                    clusterP[i] /= clusterPop[i];
                } else {
                    clusterP[i] = Double.NaN;
                    nEmpty++;
                    // orphaned cluster.  move it to a center of a particle that was not
                    // the only one in a cluster
                    int j = random.nextInt(allData.size());
                    while (clusters[j] == -1 || clusterPop[clusters[j]] < 2) {
                        j = random.nextInt(allData.size());
                    }
                    System.arraycopy(allData.get(j), 0, centers[i], 0, dim);
                }
            }
            if (nEmpty > 0) System.out.println(nEmpty + " empty");
            System.out.println(nChanged + "/" + (allData.size()) + " changed");
            if (nChanged == 0) break;
        }
        int[] clusterSizes = new int[clusters.length];
        for (int i = 0; i < centers.length; i++) {
            clusterSizes[clusterPop[i]]++;
        }
        int[] totalBonds = new int[centers.length];
        int totalTotalBonds = 0;
//        int[] clusterSet = new int[centers.length];
//        List<Set<Integer>> clustersSets = new ArrayList<>();
//        for (int i=0; i<clusterSet.length; i++) clusterSet[i] = 0;
        try {
            FileWriter fw = new FileWriter("G.dot");
            fw.write("graph G {\n");
            int nearHops = 0;
            for (int i = 0; i < centers.length; i++) {
                int pop = 100 * clusterPop[i] * dim / allData.size();
                double lnPop = Math.log(clusterPop[i] / (allData.size() * Math.log(2)));
                fw.write("c" + i + " [label=\"c" + i + "\",pop=\"" + String.format("%3.1f", lnPop) + "\",p=" + String.format("%3.1f", clusterP[i]) + "]\n");
                for (int jj = 0; jj < iBonds[i].length; jj++) {
                    int j = iBonds[i][jj];
                    double d2 = 0;
                    for (int k = 0; k < dim; k++) {
                        double dk = centers[j][k] - centers[i][k];
                        d2 += dk * dk;
                    }
                    if (d2 / dim < nbrDistance * nbrDistance) {
                        int check = jj;
                        // we made at least 1 transition between i and j
                        int hopFrac = iBondCount[i][check] * 200 / (clusterPop[i] + clusterPop[j]);
                        nearHops += iBondCount[i][check];
                        int lnHop = (int) Math.log(iBondCount[i][check] / (double) allData.size());
                        if (lnHop < -2000) {
                            throw new RuntimeException("oops " + i + " " + check + " " + iBonds[i][check] + " " + iBondCount[i][check]);
                        }
                        fw.write("c" + i + " -- c" + j + " [label=\"" + lnHop + "\"]\n");
                        totalBonds[i]++;
                        totalBonds[j]++;
                        totalTotalBonds++;
                    }
                }
            }
            int forcedBonds = 0;
            for (int i = 0; i < centers.length; i++) {
                if (totalBonds[i] == 0) {
                    double d2min = Double.POSITIVE_INFINITY;
                    int jMin = -1;
                    for (int jj = 0; jj < iBonds[i].length; jj++) {
                        int j = iBonds[i][jj];
                        double d2 = 0;
                        for (int k = 0; k < dim; k++) {
                            double dk = centers[j][k] - centers[i][k];
                            d2 += dk * dk;
                        }
                        if (d2 < d2min) {
                            d2min = d2;
                            jMin = j;
                        }
                    }
                    for (int j = 0; j < i; j++) {
                        int check = Arrays.binarySearch(iBonds[j], i);
                        if (check < 0) continue;
                        double d2 = 0;
                        for (int k = 0; k < dim; k++) {
                            double dk = centers[j][k] - centers[i][k];
                            d2 += dk * dk;
                        }
                        if (d2 < d2min) {
                            d2min = d2;
                            jMin = j;
                        }
                    }
                    // we made at least 1 transition between i and j
                    int ii = i;
                    int jj = jMin;
                    if (ii > jj) {
                        jj = i;
                        ii = jMin;
                    }
                    int check = Arrays.binarySearch(iBonds[ii], jj);
                    int hopFrac = iBondCount[ii][check] * 200 / (clusterPop[ii] + clusterPop[jj]);
                    int lnHop = (int) Math.log(iBondCount[ii][check] / (double) allData.size());
//                    nearHops += iBondCount[ii][check];
                    fw.write("c" + ii + " -- c" + jj + " [label=\"" + lnHop + "\"]\n");
                    totalBonds[ii]++;
                    totalBonds[jj]++;
                    forcedBonds++;
                }
            }
            fw.write("}\n");
            fw.close();
            System.out.println(((double) nearHops) / totalBondCount + " are near");
            System.out.println((totalTotalBonds + forcedBonds) + " bonds, " + (forcedBonds) + " forced");
        } catch (IOException ex) {
            throw new RuntimeException(ex);
        }
        System.out.println("most recent config in cluster " + clusters[allData.size() - 1]);
//        for (int i=0; i<centers.length; i++) {
//            if (totalBonds[i] == 0) {
//                System.out.println("c"+i+" unbonded");
//            }
//        }
//        for (int k=0; k<clusterSizes.length; k++) {
//            if (clusterSizes[k]!=0) System.out.println(k+" "+clusterSizes[k]);
//        }
    }

    @Override
    public void putDataInfo(IDataInfo inputDataInfo) {
        int n = inputDataInfo.getLength();
        if (n != centers[0].length) {
            for (int i = 0; i < centers.length; i++) {
                centers[i] = new double[n];
            }
            reset();
        }
    }

    public PressureSink makePressureSink() {
        return new PressureSink();
    }

    public class PressureSink implements IDataSink {
        public void putDataInfo(IDataInfo inputDataInfo) {
        }

        public void putData(IData data) {
            pSum += data.getValue(0);
            nSamplesP++;
        }
    }
}
