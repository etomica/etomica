/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import etomica.action.IAction;
import etomica.api.IRandom;
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
public class MeterVirialEBinMultiThreaded implements IAction {

    protected final int n;
    protected final IRandom random;
    protected final ClusterWheatleyExtendSW targetCluster;
    protected final MathContext mc = new MathContext(40);
    protected BoxCluster box;
    protected final Map<IntSet,MyData> allMyData;
    protected static double tRatio;
    protected long nextReweightStep = 100000L;
    protected final long[] totalCount;
    protected final PropertyBin property;
    protected double nominalWeight = 1;
    protected int iThread;
    protected static boolean quiet = false;
    protected boolean doCov = false;

    /**
     * Constructor for MeterVirial.
     */
    public MeterVirialEBinMultiThreaded(ClusterWheatleyExtendSW targetCluster, IRandom random, PropertyBin prop, int n) {
        this(targetCluster, random, prop, new long[1], null, 0, true, n);
    }

    public MeterVirialEBinMultiThreaded(ClusterWheatleyExtendSW targetCluster, IRandom random, PropertyBin prop, long[] totalCount, Map<IntSet,MyData> allMyData, int iThread, boolean doReweight, int n) {
        this.n = n;
        this.targetCluster = targetCluster;
        this.random = random;
        this.allMyData = allMyData == null ? new HashMap<IntSet,MyData>() : allMyData;
        property = prop;
        this.totalCount = totalCount;
        this.iThread = iThread;
        if (!doReweight) nextReweightStep = Long.MAX_VALUE;
    }

    public boolean getDoCov() {
        return doCov;
    }

    public void setDoCov(boolean newDoCov) {
        this.doCov = newDoCov;
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
                amd = makeData(1+n*(n-1)/2);
                amd.weight = nominalWeight;
                allMyData.put(pvCopy, amd);
            }

            amd.unscreenedCount++;
        }
        double myWeight = amd.weight;
        if (myWeight < 1 && myWeight < random.nextDouble()) {
            return;
        }
        double[] v = targetCluster.calcValue(box);
//        System.out.println("value "+Arrays.toString(v));
        double pi = -1;
        for (int i=0; i<v.length; i++) {
            if (pi==-1) {
            	pi = box.getSampleCluster().value(box);
//            	System.out.println("pi "+pi);
            }
            v[i] /= pi;
        }

//        System.out.println(propValue+" "+Arrays.toString(v));
        synchronized (amd) {
            // synchronize to prevent recomputeWeights from reading data now
            amd.addData(v);
        }
    }
    
    public void writeData(String filename) {
        long tc = 0;
        for (int i=0; i<totalCount.length; i++) {
            tc += totalCount[i];
        }
        writeData(filename, allMyData, tc, n);
    }

    public static void writeData(String filename, Map<IntSet,MyData> allMyData, long totalCount, int n) {
        try {
            FileWriter fw = new FileWriter(filename);
            fw.write(""+totalCount+"\n");
            List<IntSet> pvs = new ArrayList<IntSet>();
            pvs.addAll(allMyData.keySet());
            Collections.sort(pvs);
            int nn = 1+n*(n-1)/2;
            double[] zsum = new double[nn];
            int nnn = nn*(nn-1)/2;
            double[] zPairSum = new double[nnn];
            for (IntSet pv : pvs) {
                MyData amd = allMyData.get(pv);
                if (amd.unscreenedCount==0) continue;
                fw.write(pv+" "+amd.unscreenedCount+" "+amd.sampleCount);
                double[] s = amd.sum != null ? amd.sum : zsum;
                double[] s2 = amd.sum2 != null ? amd.sum2 : zsum;
                for (int i=0; i<nn; i++) {
                	fw.write(" "+s[i]+" "+s2[i]);
                }
                if (amd instanceof MyDataCov) {
                    double[] pairSum = amd.sum != null ? ((MyDataCov)amd).pairSum : zPairSum;
                    for (int i=0; i<nnn; i++) {
                        fw.write(" "+pairSum[i]);
                    }
                }
                fw.write("\n");
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

    public void mergeData(Map<IntSet,MeterVirialEBinMultiThreaded.MyData> moreData) {
        Set<IntSet> pvs = moreData.keySet();
        for (IntSet pv : pvs) {
            MyData amd = allMyData.get(pv);
            if (amd == null) {
                amd = makeData(1+n*(n-1)/2);
                allMyData.put(pv, amd);
            }
            MyData amdMore = moreData.get(pv);
            // null sum just means the bin was never visited
            // it was just there to hold the weight.
            if (amdMore.sum == null) continue;
            if (amd.sum != null) {
                for (int i=0; i<1+n*(n-1)/2; i++) {
    	            amd.sum[i] += amdMore.sum[i];
    	            amd.sum2[i] += amdMore.sum2[i];
                }
                if (amd instanceof MyDataCov) {
                    double[] pairSum = ((MyDataCov)amd).pairSum;
                    double[] pairSumMore = ((MyDataCov)amdMore).pairSum;
                    for (int i=0; i<pairSum.length; i++) {
                        pairSum[i] += pairSumMore[i];
                    }
                }
                amd.sampleCount += amdMore.sampleCount;
                amd.unscreenedCount += amdMore.unscreenedCount;
            }
            else {
                allMyData.put(pv, amdMore);
            }
        }
    }

    public void readData(String[] filenames, int n) {
        Map<IntSet,double[]> sums = new HashMap<IntSet,double[]>();
        Map<IntSet,double[]> sumSquares = new HashMap<IntSet,double[]>();
        Map<IntSet,double[]> pairSums = new HashMap<IntSet,double[]>();
        Map<IntSet,Long> sampleCounts = new HashMap<IntSet,Long>();
        try {
            boolean first = true;
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
                    if (first) {
                        int nn = 1+n*(n-1)/2;
                        int nnn = nn*(nn-1)/2;
                        if (values.length == 2+2*nn) doCov = false;
                        else if (values.length == 2+2*nn+nnn) doCov = true;
                        else {
                            bufReader.close();
                            throw new RuntimeException("I expect to see "+nn+" values for !doCov and "+(2+2*nn+nnn)+" values for doCov, but I actually found "+values.length+" values");
                        }
                        first = false;
                    }
                    long usc = Long.parseLong(values[0]);
                    long sampleCount = Long.parseLong(values[1]);
                    if (allMyData.containsKey(pv)) {
                        MyData amd = allMyData.get(pv);
                        amd.unscreenedCount += usc;
                        sampleCounts.put(pv, sampleCounts.get(pv)+sampleCount);
                        double[] x = sums.get(pv);
                        double[] x2 = sumSquares.get(pv);
                        for (int i=0; i<1+n*(n-1)/2; i++) {
                        	x[i] += Double.parseDouble(values[2+2*i]);
                        	x2[i] += Double.parseDouble(values[2+2*i+1]);
                        }
                        if (amd instanceof MyDataCov) {
                            int nn = (1+n*(n-1)/2);
                            int nnOffset = 2+2*nn;
                            double[] pairSum = pairSums.get(pv);
                            for (int i=0; i<pairSum.length; i++) {
                                pairSum[i] += Double.parseDouble(values[nnOffset+i]);
                            }
                        }
                    }
                    else {
                        MyData amd = makeData(1+n*(n-1)/2);
                        amd.unscreenedCount = usc;
                        amd.weight = nominalWeight;
                        allMyData.put(pv, amd);
                        sampleCounts.put(pv, sampleCount);
                        double[] x = new double[1+n*(n-1)/2];
                        double[] x2 = new double[1+n*(n-1)/2];
                        for (int i=0; i<1+n*(n-1)/2; i++) {
                        	x[i] = Double.parseDouble(values[2+2*i]);
                        	x2[i] = Double.parseDouble(values[2+2*i+1]);
                        }
                        if (amd instanceof MyDataCov) {
                            int nn = (1+n*(n-1)/2);
                            int nnOffset = 2+2*nn;
                            int nnn = nn*(nn-1)/2;
                            double[] pairSum = new double[nnn];
                            for (int i=0; i<pairSum.length; i++) {
                                pairSum[i] += Double.parseDouble(values[nnOffset+i]);
                            }
                            pairSums.put(pv, pairSum);
                        }
                        sums.put(pv, x);
                        sumSquares.put(pv, x2);
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
            if (amd instanceof MyDataCov) {
                ((MyDataCov)amd).pairSum = pairSums.get(pv);
            }
        }
    }

    public void readWeights(String filename, int n) {
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
                MyData amd = makeData(1+n*(n-1)/2);
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
        recomputeWeights(allMyData, tc, true, n);
        nextReweightStep = tc * 2;
    }

    public static void recomputeWeights(Map<IntSet,MyData> allMyData, long totalCount, int n) {
        recomputeWeights(allMyData, totalCount, true, n);
    }
    
    /**
     * This method will read MyData.unscreenedCount and MyData.accumulator
     * and write to MyData.weight.
     */
    public static void recomputeWeights(Map<IntSet,MyData> allMyData, long totalCount, boolean doPadVar, int n){
        //MeterVirialBDBinMultiThreadedOld.recomputeWeights(allMyData, totalCount, doPadVar);
        // tRatio is the ratio of the time needed to compute the biconnected
        // value (and reference value) for one configuration to the time needed
        // to generate a configuration and screen it (and any other overhead).

        for (MyData amd : allMyData.values()) {
        	amd.weight = 0;
        }
        double E0 = 0;
        long totalSampleCount = 0;
        double E1 = 0;
        double totTotalSqValue = 0;
        // weight each bin based on its contribution to the final virial coefficient at Y=1
        for (int i=0; i<1+n*(n-1)/2; i++) {
            double totalSqValue = 0;
            for (MyData amd : allMyData.values()) {
                long sc = 0;
                double average = 0;
                double var = 0;
                synchronized (amd) {
                    sc = amd.sampleCount;
                    average = amd.getAvg(i);
                    var = amd.getVar(i);
                }
                if (i==0) totalSampleCount += sc;
                totalSqValue += sc * (var + average*average);
            }
            totTotalSqValue += totalSqValue;
            double avgSqValue = totalSqValue / totalSampleCount;
            double E0a = 0, E0a2 = 0;
            // E0 = sum(sci*(steps-sci)/steps * ai^2)
            // E1 = sum(sci*sci*stdev*stdev/sampci)
        
            for (IntSet pv : allMyData.keySet()) {
                MyData amd = allMyData.get(pv);
                long c = amd.unscreenedCount;
                if (c == 0) continue;
                long sampleCount = 0;
                double average = 0;
                double var = 0;
                synchronized (amd) {
                    sampleCount = amd.sampleCount;
                    average = amd.getAvg(i);
                    var = amd.getVar(i);
                }
                if (average != 0) {
                    // E0 = sum(sci*(steps-sci)/steps * ai^2)
                    E0a += c*average;
                    E0a2 += c*average*average;
                }

                if (doPadVar) {
                    amd.weight += avgSqValue/sampleCount;
                }
                if (sampleCount<2) {
                    // we have never seen i bonds, or the configuration was always screened
                    // or we just have no statistics
                    continue;
                }

                amd.weight += var;

                // E1 = sum(sci*sci*stdev*stdev/sampci)
                E1 += c*((double)c)/sampleCount * var;
            }

            double E0ave = E0a/totalCount;
            // subtract 1 here to force E0 to be finite, even if sample is perfect so far
            double iE0 = E0a2/(totalCount-1) - E0ave*E0ave;
            E0 += iE0;
    	}
    	
    	// we've considered contributions to all Y-expansion coefficients, now
    	// we can compue weights for individual bins

        E1 /= totalCount;
//      sum1 /= totalCount;
        if (E1 == 0 && doPadVar) {
          // no value fluctuations, perhaps B4 or B5?
          E1 = (totTotalSqValue/totalSampleCount)/totalSampleCount;
        }

        double k = Math.sqrt(1/(E0*tRatio));

        double newT1 = 0;
        long totalUnscreened = 0;
        double newE1 = 0;
        double E1all = 0;
        double allT1 = 0;
        double t0 = totalCount;
        double t1 = totalSampleCount*tRatio;
        for (IntSet pv : allMyData.keySet()) {
            MyData amd = allMyData.get(pv);
            long c = amd.unscreenedCount;
            if (c == 0) {
                amd.weight = 1;
                continue;
            }
            double var = amd.weight;
            double w = Math.sqrt(var)*k;

            if (w > 1 || amd.sampleCount < 2) {
                w = 1;
            }
//            if (i>=targetCluster.pointCount()) System.out.println(String.format("%2d %12d  %6.4f  %6.4f  %6.4f", i, (count[i]-screenedCount[i]), weight[i], w, tRatio * (count[i]-screenedCount[i]) * w / totalCount));
            amd.weight = w;
            newT1 += c * w;
            allT1 += c;
            totalUnscreened += c;
            if (var > 0) {
                newE1 += c*var/w;
                E1all += c*var;
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
            System.out.print(String.format("   Difficulty: %10.4e %10.4e %10.4e %10.4e\n", Math.sqrt((E0+newE1)*(1+newT1)), Math.sqrt((E0+E1)*(1+t1/t0)), Math.sqrt((E0+E1all)*(1+allT1)), Math.sqrt(E0+E1all)));
            // D = printed * sqrt((ts*totalCount)/totalCount)
            //   = printed * sqrt(ts)
        }


        System.out.flush();
    }
    
    public ClusterWheatleyExtendSW getTargetCluster() {
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
    
    public MathContext getMathContext() {
        return mc;
    }
    
    public long getTotalCount() {
        long tc = 0;
        for (int i=0; i<totalCount.length; i++) {
            tc += totalCount[i];
        }
        return tc;
    }
    
    protected MyData makeData(int n) {
        return doCov ? new MyDataCov(n) : new MyData(n);
    }

    public static class MyData {
        public long unscreenedCount, sampleCount;
        public double weight;
        public byte n;
        public double[] sum, sum2;

        public MyData(int n) {
            this.n = (byte)n;
        }

        public double getAvg(int i) {
            if (sampleCount < 1) return 0;
            return sum[i]/sampleCount;
        }

        public double getVar(int i) {
            if (sampleCount < 1) return Double.NaN;
            double avg = getAvg(i);
            double avg2 = avg*avg;
            double var = sum2[i]/sampleCount - avg2;
            if (var < avg2*1e-7) var = 0;
            return var;
        }

        public void addData(double[] v) {
            if (sum == null) {
                sum = new double[n];
                sum2 = new double[n];
            }
        	for (int i=0; i<v.length; i++) {
        		sum[i] += v[i];
        		sum2[i] += v[i]*v[i];
        	}
            sampleCount++;
        }
    }

    public static class MyDataCov extends MyData {
        public double[] pairSum;

        public MyDataCov(int n) {
            super(n);
        }

        public void addData(double[] v) {
            if (pairSum==null) {
                pairSum = new double[n*(n-1)/2];
            }
            super.addData(v);
            int k = 0;
            for (int i=0; i<n-1; i++) {
                for (int j=i+1; j<n; j++) {
                    pairSum[k] += v[i]*v[j];
                    k++;
                }
            }
        }

        public double getCov(int i, int j) {
            if (i==j) return getVar(i);
            if (i>j) return getCov(j,i);
            int k = (2*n-i-1)*i/2 + (j-i-1);
            double avgi = sum[i]/sampleCount;
            double avgj = sum[j]/sampleCount;
            double cov = pairSum[k]/sampleCount - avgi*avgj;
            if (Math.abs(cov/(avgi*avgj)) < 1e-7) cov = 0;
            return cov;
        }
    }
}
