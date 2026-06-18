/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.space.Vector;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Time;
import etomica.util.Statefull;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.Arrays;

public class DataSourceStrings implements IDataSource, ConfigurationStorage.ConfigurationStorageListener, DataSourceIndependent, Statefull {

    protected final ConfigurationStorage configStorage;
    protected DataDoubleArray tData;
    protected DataDoubleArray.DataInfoDoubleArray tDataInfo;
    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected final DataTag tTag, tag;
    protected final AtomTestDeviation atomTestDeviation;
    protected final double[][] dr2;
    protected final int numAtoms;
    protected final double mobFrac, strTol2;
    protected final Box box;
    protected long[] nStrings;
    protected long[] numAtomInString;
    protected final Vector dr;
    protected double nbrMax2 = 1.5 * 1.5;
    protected final int[] strings;
    protected final int[] nextAtom, firstAtoms;
    protected int log2StepMin;

    public DataSourceStrings(ConfigurationStorage configStorage, int log2StepMin) {
        this.configStorage = configStorage;
        this.log2StepMin = log2StepMin;
        box = configStorage.getBox();
        numAtoms = box.getLeafList().size();
        nStrings = new long[0];
        numAtomInString = new long[0];
        dr2 = new double[numAtoms][2];
        tag = new DataTag();
        tTag = new DataTag();
        atomTestDeviation = new AtomTestDeviation(box, configStorage);
        mobFrac = 0.065;
        strTol2 = 0.3*0.3;
        dr = box.getSpace().makeVector();
        int n = box.getLeafList().size();
        strings = new int[n];
        firstAtoms = new int[n];
        nextAtom = new int[n];
        reallocate(0);
    }

    public void reallocate(int n) {
        if (n == numAtomInString.length && data != null) return;
        if (n < 1) n = 0;
        numAtomInString = Arrays.copyOf(numAtomInString, n);
        nStrings = Arrays.copyOf(nStrings, n);
        data = new DataFunction(new int[]{n});
        tData = new DataDoubleArray(new int[]{n});
        tDataInfo = new DataDoubleArray.DataInfoDoubleArray("t", Time.DIMENSION, new int[]{n});
        tDataInfo.addTag(tTag);
        dataInfo = new DataFunction.DataInfoFunction("String", new CompoundDimension(new Dimension[]{Length.DIMENSION}, new double[]{2}), this);
        dataInfo.addTag(tag);
        double[] t = tData.getData();
        if (t.length > 0) {
            double dt = configStorage.getDeltaT();
            for (int i = 0; i < t.length; i++) {
                t[i] = dt * (1L << i);
            }
        }
    }

    @Override
    public IData getData() {
        if (configStorage.getLastConfigIndex() < 1) return data;
        double[] y = data.getData();
        for (int i = 0; i < numAtomInString.length; i++) {
            y[i] = ((double)numAtomInString[i]) / nStrings[i];
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
        long step = configStorage.getSavedSteps()[0];
        if (step % (1L << log2StepMin) != 0) return;
        IAtomList atoms = box.getLeafList();

        for (int i = 0; i < configStorage.getLastConfigIndex(); i++) {
            int x = Math.max(log2StepMin, i);
            if (step % (1L << x) == 0) {
                if (i >= nStrings.length) reallocate(i + 1);
                atomTestDeviation.setConfigIndex(i + 1);
                for (int j = 0; j < numAtoms; j++) {
                    dr2[j][0] = j;
                    dr2[j][1] = atomTestDeviation.getDisplacementSq(atoms.get(j));
                }

                java.util.Arrays.sort(dr2, 0, numAtoms, new java.util.Comparator<double[]>() {
                    public int compare(double[] a, double[] b) {
                        return Double.compare(b[1], a[1]);
                    }
                });

                findStrings(i);

                for (int j = 0; j < firstAtoms.length; j++) {
                    if (firstAtoms[j] == -1) break;
                    nStrings[i]++;
                    for (int ii = firstAtoms[j]; ii != -1; ii = nextAtom[ii]) {
                        numAtomInString[i]++;
                    }
                }
            }
        }
    }

    public void findStrings(int interval) {
        Vector[] positions = configStorage.getSavedConfig(0);
        Vector[] oldPositions = configStorage.getSavedConfig(interval + 1);
        int[] replacedAtom = new int[numAtoms];

        for (int i = 0; i < numAtoms*mobFrac; i++) {
            int ii = (int)dr2[i][0];
            nextAtom[ii] = strings[ii] = replacedAtom[ii] = -1;
        }
        int nClusters = 0;
        firstAtoms[0] = -1;
        for (int ii = 0; ii < numAtoms*mobFrac; ii++) {
            int i = (int)dr2[ii][0];
            Vector ri = positions[i];
            Vector riOld = oldPositions[i];
            for (int jj = 0; jj < numAtoms*mobFrac; jj++) {
                if(ii == jj) continue;
                int j = (int)dr2[jj][0];
                Vector rjOld = oldPositions[j];
                dr.Ev1Mv2(rjOld, riOld);
                box.getBoundary().nearestImage(dr);
                if (dr.squared() > nbrMax2) continue;

                dr.Ev1Mv2(ri, rjOld); //i(t) --> j(0)
                box.getBoundary().nearestImage(dr);
                if(dr.squared() > strTol2) continue;

                replacedAtom[j] = i;
//                System.out.println(interval + " " + i +" replaced " + j);
                if (strings[j] != -1) {
                    nextAtom[i] = j;
//                    System.out.println("nextAtom["+i+"]" + " is now" + j);
                    int iCluster = strings[j];
                    strings[i] = iCluster;
                    firstAtoms[iCluster] = i;
                    if(replacedAtom[i]>-1){
                        int jCluster = strings[replacedAtom[i]];
                        //we have detected a ring
                        if(iCluster == jCluster) continue;
                        nextAtom[replacedAtom[i]] = i;
//                        System.out.println("nextAtom["+replacedAtom[i]+"]" + " is now" + i);
                        for (int m = firstAtoms[iCluster]; m != -1; m = nextAtom[m]) {
                            strings[m] = jCluster;
                        }
                        firstAtoms[iCluster] = -1;

                        int lastCluster = nClusters - 1;
                        if (lastCluster > iCluster) {
                            // fill in our hole
                            firstAtoms[iCluster] = firstAtoms[lastCluster];
                            for (int m = firstAtoms[iCluster]; m != -1; m = nextAtom[m]) {
                                strings[m] = iCluster;
                            }
                        }
                        firstAtoms[lastCluster] = -1;
                        nClusters--;
                    }
                }else if(replacedAtom[i]>-1){
                    int iCluster = strings[replacedAtom[i]];
                    strings[i] = iCluster;
                    nextAtom[replacedAtom[i]] = i;
//                    System.out.println("nextAtom["+replacedAtom[i]+"]" + " is now" + i);
                }else{
                    strings[i] = nClusters;
                    firstAtoms[nClusters] = i;
                    if (nClusters < firstAtoms.length - 1) firstAtoms[nClusters + 1] = -1;
                    nClusters++;
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

    @Override
    public void saveState(Writer fw) throws IOException {
        fw.write(nStrings.length+"\n");
        for (int i=0; i<nStrings.length; i++) {
            fw.write(nStrings[i]+" "+numAtomInString[i]+"\n");
        }
    }

    @Override
    public void restoreState(BufferedReader br) throws IOException {
        int n = Integer.parseInt(br.readLine());
        reallocate(n);
        for (int i=0; i<n; i++) {
            String[] bits = br.readLine().split(" ");
            nStrings[i] = Long.parseLong(bits[0]);
            numAtomInString[i] = Long.parseLong(bits[1]);
        }
    }
}
