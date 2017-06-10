/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import etomica.action.IAction;
import etomica.util.random.IRandom;
import etomica.virial.IntSet.PropertyBin;

/**
 * Measures cluster averages for virial coefficients.  Configurations are
 * classified by a "property" (an array of ints) and assigned a bin.  The
 * average value and standard deviation within that bin is tracked.  If a bin
 * has a small standard deviation, its configurations may be skipped (not
 * measured) in order to spend more time on other bins, or to track how often
 * each bin is visited.
 * 
 * Additionally, this class can be used in a threaded simulation.  Each thread
 * runs a separate simulation, sharing a Map of MyData.  This class ensures
 * that the various instances of the meter do not concurrently use the Map in
 * a way that causes problems.  Many of the operations need only be approximate
 * and the class avoids forcing synchronization for such operations.
 */
public class MeterVirialBinMultiThreaded implements IAction {

    protected final IRandom random;
    protected final ClusterWheatley targetCluster;
    protected BoxCluster box;
    protected final Map<IntSet,MyData> allMyData;
    protected static double tRatio;
    protected long nextReweightStep = 100000L;
    protected final long[] totalCount;
    protected final PropertyBin property;
    protected double nominalWeight = 1;
    protected int iThread;
//    protected final BufferedReader goodBufReader;
//    protected long lastCPairID, lastLastCPairID;
    protected boolean excludeBogusConfigs = false;
    protected static boolean quiet = false;

    /**
     * Constructor for MeterVirial.
     */
    public MeterVirialBinMultiThreaded(ClusterWheatley targetCluster, IRandom random, PropertyBin prop) {
        this(targetCluster, random, prop, new long[1], null, 0, true);
    }

    public MeterVirialBinMultiThreaded(ClusterWheatley targetCluster, IRandom random, PropertyBin prop, long[] totalCount, Map<IntSet,MyData> allMyData, int iThread, boolean doReweight) {
        System.out.println("Bin2");
        this.targetCluster = targetCluster;
        this.random = random;
        this.allMyData = allMyData == null ? new HashMap<IntSet,MyData>() : allMyData;
        property = prop;
        this.totalCount = totalCount;
        this.iThread = iThread;
        if (!doReweight) nextReweightStep = Long.MAX_VALUE;
        /*File f = new File("thread"+iThread+"_raw.dat");
        if (f.exists()) {
            f.delete();
        }*/
//        try {
//            FileReader fileReader = new FileReader("thread"+iThread+"_good.dat");
//            goodBufReader = new BufferedReader(fileReader);
//        }
//        catch (IOException ex) {
//            throw new RuntimeException(ex);
//        }
    }

    public void setDoExcludeBogusConfigs(boolean newDoExclude) {
        excludeBogusConfigs = newDoExclude;
    }
    
    public void setWeight(double newWeight) {
        nominalWeight = newWeight;
        nextReweightStep = Long.MAX_VALUE;
    }

    public static void setTRatio(double newTRatio) {
        tRatio = newTRatio;
    }
    
    public static void setQuiet(boolean newQuiet) {
        quiet = newQuiet;
    }

    public void actionPerformed() {
        
        if (iThread == 0) {
            long tc = 0;
            for (int i=0; i<totalCount.length; i++) {
                tc += totalCount[i];
            }
            if (tc >= nextReweightStep) {
                // this thread will be recomputing weights for all threads
                synchronized (allMyData) {
                    recomputeWeights();
                }
                nextReweightStep = tc*2;
            }
        }
        // this can have thread trouble, but it only matter if we're going to
        // reweight, and (even then) only slightly effects the outcome of
        // reweighting
        totalCount[iThread]++;

        boolean hasNonZeroValue = targetCluster.checkConfig(box);
        if (!hasNonZeroValue) {
//            lastLastCPairID = lastCPairID;
//            lastCPairID = box.getCPairSet().getID();
            return;
        }
        IntSet propValue = property.value();
        MyData amd;
        synchronized (allMyData) {
            // this needs to be synchronized to trying retrieving a value here
            // and while adding that same key/value below.
            amd = allMyData.get(propValue);
            if (amd == null) {
                IntSet pvCopy = propValue.copy();
                amd = new MyData();
                amd.weight = nominalWeight;
                allMyData.put(pvCopy, amd);
//                int foo = allMyData.size();
//                if (foo == 1<<(int)(Math.log(foo)/Math.log(2))) {
//                    System.out.println(foo+" sets, "+totalCount/foo+" steps/set");
//                }
            }

            amd.unscreenedCount++;
        }
        double myWeight = amd.weight;
        if (myWeight < 1 && myWeight < random.nextDouble()) {
//            lastLastCPairID = lastCPairID;
//            lastCPairID = box.getCPairSet().getID();
            return;
        }
        double x = 0;
        double v = targetCluster.calcValue(box);
        if (v != 0) {
//            long cPairIDNow = box.getCPairSet().getID();
            double pi = box.getSampleCluster().value(box);
//            long cPairIDNow2 = box.getCPairSet().getID();
            x = v/pi;
        }

        synchronized (amd) {
            // synchronize to prevent recomputeWeights from reading data now
            amd.addData(x);
        }
//        lastLastCPairID = lastCPairID;
//        lastCPairID = box.getCPairSet().getID();
    }
    
    public void writeData(String filename) {
        long tc = 0;
        for (int i=0; i<totalCount.length; i++) {
            tc += totalCount[i];
        }
        writeData(filename, allMyData, tc);
    }

    public static void writeData(String filename, Map<IntSet,MyData> allMyData, long totalCount) {
        try {
            FileWriter fw = new FileWriter(filename);
            fw.write(""+totalCount+"\n");
            List<IntSet> pvs = new ArrayList<IntSet>();
            pvs.addAll(allMyData.keySet());
            Collections.sort(pvs);
            for (IntSet pv : pvs) {
                MyData amd = allMyData.get(pv);
                fw.write(pv+" "+amd.unscreenedCount+" "+amd.sampleCount+" "+amd.sum+" "+amd.sum2+"\n");
            }
            fw.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public void writeWeights(String filename) {
        writeWeights(filename, allMyData);
    }

    public static void writeWeights(String filename, Map<IntSet,MyData> allMyData) {
        try {
            FileWriter fw = new FileWriter(filename);
            List<IntSet> pvs = new ArrayList<IntSet>();
            pvs.addAll(allMyData.keySet());
            Collections.sort(pvs);
            for (IntSet pv : pvs) {
                fw.write(pv+" "+allMyData.get(pv).weight+"\n");
            }
            fw.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public void mergeData(Map<IntSet,MyData> moreData) {
        Set<IntSet> pvs = moreData.keySet();
        for (IntSet pv : pvs) {
            MyData amd = allMyData.get(pv);
            if (amd == null) {
                amd = new MyData();
                allMyData.put(pv, amd);
            }
            MyData amdMore = moreData.get(pv);
            amd.sum += amdMore.sum;
            amd.sum2 += amdMore.sum2;
            amd.sampleCount += amdMore.sampleCount;
            amd.unscreenedCount += amdMore.unscreenedCount;
        }
    }

    public void readData(String[] filenames) {
        Map<IntSet,Double> sums = new HashMap<IntSet,Double>();
        Map<IntSet,Double> sumSquares = new HashMap<IntSet,Double>();
        Map<IntSet,Double> dsums = new HashMap<IntSet,Double>();
        Map<IntSet,Long> sampleCounts = new HashMap<IntSet,Long>();
        try {
            for (String filename : filenames) {
                File f = new File(filename);
                if (!f.exists()) continue;
                FileReader fr = new FileReader(filename);
                BufferedReader bufReader = new BufferedReader(fr);
                String line = bufReader.readLine();
                totalCount[iThread] += Long.parseLong(line);
                while ((line=bufReader.readLine()) != null) {
                    String pvStr = line.replaceFirst("].*", "").substring(1);
                    String[] pvSplit = pvStr.split("[, ]+");
                    int[] v = new int[pvSplit.length];
                    for (int i=0; i<v.length; i++) {
                        v[i] = Integer.parseInt(pvSplit[i]);
                    }
                    IntSet pv = new IntSet(v);
                    String[] values = line.replaceFirst(".*] ", "").split(" +");
                    long usc = Long.parseLong(values[0]);
                    long sampleCount = Long.parseLong(values[1]);
                    double sum = Double.parseDouble(values[2]);
                    double sumSquare = Double.parseDouble(values[3]);
                    if (allMyData.containsKey(pv)) {
                        MyData amd = allMyData.get(pv);
                        amd.unscreenedCount += usc;
                        sampleCounts.put(pv, sampleCounts.get(pv)+sampleCount);
                        sums.put(pv, sums.get(pv) +sum);
                        sumSquares.put(pv, sumSquares.get(pv) + sumSquare);
                    }
                    else {
                        MyData amd = new MyData();
                        amd.unscreenedCount = usc;
                        amd.weight = nominalWeight;
                        allMyData.put(pv, amd);
                        sampleCounts.put(pv, sampleCount);
                        sums.put(pv, sum);
                        sumSquares.put(pv, sumSquare);
                        dsums.put(pv, Double.parseDouble(values[2]));
                    }
                }
                bufReader.close();
            }
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
        for (IntSet pv : allMyData.keySet()) {
            MyData amd = allMyData.get(pv);
            amd.sampleCount = sampleCounts.get(pv);
            amd.sum = sums.get(pv);
            amd.sum2 = sumSquares.get(pv);
        }
    }

    public void readWeights(String filename) {
        File f = new File(filename);
        if (!f.exists()) return;
        try {
            FileReader fr = new FileReader(filename);
            BufferedReader bufReader = new BufferedReader(fr);
            String line = null;
            while ((line=bufReader.readLine()) != null) {
                String pvStr = line.replaceAll("].*", "").substring(1);
                String[] pvSplit = pvStr.split("[, ]+");
                int[] v = new int[pvSplit.length];
                for (int i=0; i<v.length; i++) {
                    v[i] = Integer.parseInt(pvSplit[i]);
                }
                IntSet pv = new IntSet(v);
                String weightStr = line.replaceAll(".*] ", "");
                MyData amd = new MyData();
                amd.weight = Double.parseDouble(weightStr);
                amd.unscreenedCount = 0;
                allMyData.put(pv, amd);
            }
            bufReader.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
        nextReweightStep = Long.MAX_VALUE;
    }

    public void recomputeWeights() {
        if (nextReweightStep == Long.MAX_VALUE) return;
        long tc = 0;
        for (int i=0; i<totalCount.length; i++) {
            tc += totalCount[i];
        }
        recomputeWeights(allMyData, tc, true);
        nextReweightStep = tc * 2;
    }

    public static void recomputeWeights(Map<IntSet,MyData> allMyData, long totalCount) {
        recomputeWeights(allMyData, totalCount, true);
    }
    
    /**
     * This method will read MyData.unscreenedCount and MyData.accumulator
     * and write to MyData.weight.
     */
    public static void recomputeWeights(Map<IntSet,MyData> allMyData, long totalCount, boolean doPadVar) {
        // tRatio is the ratio of the time needed to compute the biconnected
        // value (and reference value) for one configuration to the time needed
        // to generate a configuration and screen it (and any other overhead).
        long totalSampleCount = 0;
        double totalSqValue = 0;
        for (MyData amd : allMyData.values()) {
            long sc = 0;
            double average = 0;
            double var = 0;
            synchronized (amd) {
                sc = amd.sampleCount;
                average = amd.getAvg();
                var = amd.getVar();
            }
            totalSampleCount += sc;
            totalSqValue += sc * (var + average*average);
        }
        double avgSqValue = totalSqValue / totalSampleCount;
        double t0 = totalCount;
        double t1 = totalSampleCount*tRatio;
        double E0a = 0, E0a2 = 0;
        double E1 = 0;
        // E0 = sum(sci*(steps-sci)/steps * ai^2)
        // E1 = sum(sci*sci*stdev*stdev/sampci)

        Map<IntSet,Double> localWeight = new HashMap<IntSet,Double>();
        
        for (IntSet pv : allMyData.keySet()) {
            MyData amd = allMyData.get(pv);
            long c = amd.unscreenedCount;
            if (c == 0) continue;
            long sampleCount = 0;
            double average = 0;
            double var = 0;
            synchronized (amd) {
                sampleCount = amd.sampleCount;
                average = amd.getAvg();
                var = amd.getVar();
            }
            double lwi = doPadVar ? (avgSqValue/sampleCount) : 0;

            if (average != 0) {
                // E0 = sum(sci*(steps-sci)/steps * ai^2)
                E0a += c*average;
                E0a2 += c*average*average;
            }

            if (sampleCount<2) {
                // we have never seen i bonds, or the configuration was always screened
                // or we just have no statistics
                localWeight.put(pv, Math.sqrt(lwi));
                continue;
            }

            lwi += var;
            localWeight.put(pv, Math.sqrt(lwi));

            // E1 = sum(sci*sci*stdev*stdev/sampci)
            E1 += c*((double)c)/sampleCount * var;
        }
        if (E0a2 == 0) {
            return;
        }
        double E0ave = E0a/totalCount;
        double E0 = E0a2/totalCount - E0ave*E0ave;
        E1 /= totalCount;
        if (E1 == 0 && doPadVar) {
            // no value fluctuations, perhaps B4 or B5?
            E1 = avgSqValue/totalSampleCount;
        }
//        System.out.println("weights");
        double k = Math.sqrt(1/(E0*tRatio));
//        k=1e-8;

        double newT1 = 0;
        long totalUnscreened = 0;
        double newE1 = 0;
        double E1all = 0;
        double allT1 = 0;
        for (IntSet pv : allMyData.keySet()) {
            MyData amd = allMyData.get(pv);
            long c = amd.unscreenedCount;
            if (c == 0) continue;
            double w = localWeight.get(pv)*k;
            if (w > 1 || amd.sampleCount < 2) {
                w = 1;
            }
//            if (i>=targetCluster.pointCount()) System.out.println(String.format("%2d %12d  %6.4f  %6.4f  %6.4f", i, (count[i]-screenedCount[i]), weight[i], w, tRatio * (count[i]-screenedCount[i]) * w / totalCount));
            amd.weight = w;
            newT1 += c * w;
            allT1 += c;
            totalUnscreened += c;
            double s = amd.getVar();
            if (s > 0) {
                newE1 += c*s/w;
                E1all += c*s;
            }
        }
        newT1 *= tRatio/totalCount;
        allT1 *= tRatio/totalCount;
        newE1 /= totalCount;
        E1all /= totalCount;
//        System.out.println(E0+" "+E1+" "+newE1+" "+newT1);
        if (!quiet) {
            System.out.print(String.format("var0 frac %8.5f (opt: %8.5f)  t0 frac %8.5f  k %8.2e  new t0 frac %5.3f   measure frac %7.5f\n", E0/(E0+E1), E0/(E0+newE1), t0/(t0+t1), k, 1/(1+newT1), newT1*(totalCount/tRatio/totalUnscreened)));
            // difficulty:   opt   actual   w=1    w=0
            System.out.print(String.format(" Difficulty: %10.4e %10.4e %10.4e %10.4e\n", Math.sqrt((E0+newE1)*(1+newT1)), Math.sqrt((E0+E1)*(1+t1/t0)), Math.sqrt((E0+E1all)*(1+allT1)), Math.sqrt(E0+E1all)));
        }
        System.out.flush();
    }
    
    public ClusterWheatley getTargetCluster() {
        return targetCluster;
    }
    
    public BoxCluster getBox() {
        return box;
    }
    
    public void setBox(BoxCluster newBox) {
        box = newBox;
    }
    
    public Map<IntSet,MyData> getAllMyData() {
        return allMyData;
    }

    public long getTotalCount() {
        long tc = 0;
        for (int i=0; i<totalCount.length; i++) {
            tc += totalCount[i];
        }
        return tc;
    }

    public static class MyData {
        public long unscreenedCount, sampleCount;
        public double weight;
        public double sum, sum2;
        
        public double getAvg() {
            return sum/sampleCount;
        }
        
        public double getVar() {
            if (sampleCount < 1) return Double.NaN;
            double avg = getAvg();
            double avg2 = avg*avg;
            double var = sum2/sampleCount - avg2;
            if (var < avg2*1e-7) var = 0;
            return var;
        }
        
        public void addData(double value) {
            sum += value;
            sum2 += value*value;
            sampleCount++;
        }
    }}
