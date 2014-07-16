package etomica.virial;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import etomica.action.IAction;
import etomica.api.IRandom;
import etomica.virial.MeterVirialBDBinMulti.IntSet;
import etomica.virial.MeterVirialBDBinMulti.PropertyBin;

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
public class MeterVirialBDBinMultiThreaded implements IAction {

    protected final IRandom random;
    protected final ClusterWheatley targetCluster;
    protected final MathContext mc = new MathContext(40);
    protected BoxCluster box;
    protected final Map<IntSet,MyData> allMyData;
    protected static double tRatio;
    protected static final BigDecimal BDZERO = new BigDecimal(0);
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
    public MeterVirialBDBinMultiThreaded(ClusterWheatley targetCluster, IRandom random, PropertyBin prop) {
        this(targetCluster, random, prop, new long[1], new HashMap<IntSet,MyData>(), 0, true);
    }

    public MeterVirialBDBinMultiThreaded(ClusterWheatley targetCluster, IRandom random, PropertyBin prop, long[] totalCount, Map<IntSet,MyData> allMyData, int iThread, boolean doReweight) {
        this.targetCluster = targetCluster;
        this.random = random;
        this.allMyData = allMyData;
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
        
        if (excludeBogusConfigs) {
            if (MCMoveClusterAtomHSRing.fullWouldReject) return;
        }
        
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
        BigDecimal x = BDZERO;
        double v = targetCluster.calcValue(box);
        if (v != 0) {
//            long cPairIDNow = box.getCPairSet().getID();
            double pi = box.getSampleCluster().value(box);
//            long cPairIDNow2 = box.getCPairSet().getID();
            x = new BigDecimal(v, mc).divide(new BigDecimal(pi, mc), mc);
//            try {
//              String line = goodBufReader.readLine();
//              if (!x.toString().equals(line)) {
//                System.out.println("hi I'm thread "+iThread);
//                System.out.println("we're in bin "+propValue);
//                System.out.println("this is unscreened config "+amd.unscreenedCount);
//                System.out.println("and sample # "+(amd.sampleCount+1));
//                System.out.println("we expected "+line);
//                System.out.println("but we calculated "+x.toString());
//                long sig = ((ClusterWheatleyPartitionScreening)targetCluster).calcSignature(box);
//                System.out.println("configuration signature is "+sig);
//                System.out.println("prev prev cPairID is "+lastLastCPairID);
//                System.out.println("prev cPairID is "+lastCPairID);
//                System.out.println("cPairID is "+box.getCPairSet().getID());
//                System.out.println("cPairID now and now2 "+cPairIDNow+" "+cPairIDNow2);
//                System.out.println("target value "+v+" ("+((ClusterWheatleyPartitionScreening)targetCluster).getCPairID()+")");
//                System.out.println("previous target value "+((ClusterWheatleyPartitionScreening)targetCluster).getLastValue()+" ("+((ClusterWheatleyPartitionScreening)targetCluster).getLastCPairID()+")");
//                System.out.println("ref value "+pi+" ("+((ClusterChainHS)((ClusterWeightUmbrella)box.getSampleCluster()).getClusters()[0]).getCPairID()+","+((ClusterSinglyConnected)((ClusterWeightUmbrella)box.getSampleCluster()).getClusters()[1]).getCPairID()+")");
//                System.out.println("previous ClusterChainHS value "+((ClusterChainHS)((ClusterWeightUmbrella)box.getSampleCluster()).getClusters()[0]).getLastValue()+" ("+((ClusterChainHS)((ClusterWeightUmbrella)box.getSampleCluster()).getClusters()[0]).getLastCPairID()+")");
//                System.out.println("previous ClusterSC value "+((ClusterSinglyConnected)((ClusterWeightUmbrella)box.getSampleCluster()).getClusters()[1]).getLastValue()+" ("+((ClusterSinglyConnected)((ClusterWeightUmbrella)box.getSampleCluster()).getClusters()[1]).getLastCPairID()+")");
//                System.out.println("individual ref values");
//                ((ClusterWeightUmbrella)box.getSampleCluster()).DEBUGME = true;
//                double newpi = box.getSampleCluster().value(box);
//                System.out.println("recalc ref value is "+newpi);
//                // force recalculation
//                box.trialNotify();
//                ((ClusterChainHS)((ClusterWeightUmbrella)box.getSampleCluster()).getClusters()[0]).DEBUGME = true;
//                System.out.println("forced recalculation");
//                double v2 = targetCluster.calcValue(box);
//                System.out.println("target value "+v2);
//                System.out.println("individual ref values");
//                newpi = box.getSampleCluster().value(box);
//                System.out.println("recalc ref value is "+newpi);
//                System.out.println("new ratio is "+v2/newpi);
//                box.acceptNotify();
//                throw new RuntimeException("oops");
//              }
//            }
//            catch (IOException ex) {
//              throw new RuntimeException(ex);
//            }
        }
//        else {
//            try {
//              String line = goodBufReader.readLine();
//              if (!x.toString().equals("0")) {
//                System.out.println("hi I'm thread "+iThread);
//                System.out.println("we're in bin "+propValue);
//                System.out.println("this is unscreened config "+amd.unscreenedCount);
//                System.out.println("and sample # "+(amd.sampleCount+1));
//                System.out.println("we expected "+line);
//                System.out.println("but we calculated 0");
//                long sig = ((ClusterWheatleyPartitionScreening)targetCluster).calcSignature(box);
//                System.out.println("configuration signature is "+sig);
//                box.trialNotify();
//                System.out.println("forced recalculation");
//                double v2 = targetCluster.calcValue(box);
//                System.out.println("target value "+v2);
//                box.acceptNotify();
//                throw new RuntimeException("oops");
//              }
//            }
//            catch (IOException ex) {
//              throw new RuntimeException(ex);
//            }
//        }
        
//        if (propValue.v[0] == 33 && propValue.v[1] == 197) {
//            System.out.println("here we go!");
//            for (int i=0; i<1462563; i++) {
//                amd.addData(x, mc);
//                double var = amd.getAvgDouble();
//                if (var > 0) {
//                    throw new RuntimeException(i+" oops");
//                }
//            }
//            System.out.println("success!");
//            System.exit(1);
//        }

        synchronized (amd) {
            // synchronize to prevent recomputeWeights from reading data now
            amd.addData(x, mc, propValue, iThread);
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
                fw.write(pv+" "+amd.unscreenedCount+" "+amd.sampleCount+" "+amd.sum.toString()+" "+amd.sum2.toString()+" "+amd.dsum+"\n");
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

    public void readData(String[] filenames) {
        Map<IntSet,BigDecimal> sums = new HashMap<IntSet,BigDecimal>();
        Map<IntSet,BigDecimal> sumSquares = new HashMap<IntSet,BigDecimal>();
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
                    BigDecimal sum = new BigDecimal(values[2], mc);
                    BigDecimal sumSquare = new BigDecimal(values[3], mc);
                    if (allMyData.containsKey(pv)) {
                        MyData amd = allMyData.get(pv);
                        amd.unscreenedCount += usc;
                        sampleCounts.put(pv, sampleCounts.get(pv)+sampleCount);
                        sums.put(pv, sums.get(pv).add(sum, mc));
                        sumSquares.put(pv, sumSquares.get(pv).add(sumSquare, mc));
                        dsums.put(pv, dsums.get(pv)+Double.parseDouble(values[4]));
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
            amd.dsum = dsums.get(pv);
        }
    }
    
    public static void mergeData(MathContext mc, Map<IntSet,MeterVirialBDBinMultiThreaded.MyData> masterData, Map<IntSet,MeterVirialBDBinMultiThreaded.MyData> singleData) {
        Set<IntSet> pvs = singleData.keySet();
        for (IntSet pv : pvs) {
            MyData amdMaster = masterData.get(pv);
            if (amdMaster == null) {
                amdMaster = new MyData();
                masterData.put(pv, amdMaster);
            }
            MyData amdSingle = singleData.get(pv);
            amdMaster.dsum += amdSingle.dsum;
            amdMaster.dsum2 += amdSingle.dsum2;
            amdMaster.sampleCount += amdSingle.sampleCount;
            amdMaster.sum = amdMaster.sum.add(amdSingle.sum, mc);
            amdMaster.sum2 = amdMaster.sum2.add(amdSingle.sum2, mc);
            amdMaster.unscreenedCount += amdSingle.unscreenedCount;
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
        recomputeWeights(allMyData, tc);
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
                average = amd.getAvgDouble();
                var = amd.getVarDouble();
            }
            totalSampleCount += sc;
            totalSqValue += sc * (var + average*average);
        }
        double avgSqValue = totalSqValue / totalSampleCount;
        double t0 = totalCount;
        double t1 = totalSampleCount*tRatio;
        double E0 = 0;
        double E0a = 0, E0a2 = 0;
        double E1 = 0;
        // E0 = sum(sci*(steps-sci)/steps * ai^2)
        // E1 = sum(sci*sci*stdev*stdev/sampci)

        double maxw = 0;
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
                average = amd.getAvgDouble();
                var = amd.getVarDouble();
            }
            double lwi = doPadVar ? avgSqValue/sampleCount : 0;

            if (average != 0) {
                // E0 = sum(sci*(steps-sci)/steps * ai^2)
                E0 += c*((double)(totalCount-c))/totalCount * average*average;
                E0a += c*average;
                E0a2 += c*average*average;
            }

            if (sampleCount<2) {
                // we have never seen i bonds, or the configuration was always screened
                // or we just have no statistics
                localWeight.put(pv, lwi);
                continue;
            }

            if (sampleCount > 1) {
                lwi += var;
                if (sampleCount > 10 && lwi > maxw) maxw = lwi;
            }
            localWeight.put(pv, lwi);

            // E1 = sum(sci*sci*stdev*stdev/sampci)
            E1 += c*((double)c)/sampleCount * var;
        }
        double E0ave = E0a/totalCount;
        double E00 = E0a2/totalCount - E0ave*E0ave;
        if (E0 == 0) {
            return;
        }
        if (maxw==0) {
            // nonsense, but whatever
            maxw = avgSqValue/totalSampleCount;
        }
        if (E1 == 0 && doPadVar) {
            // no value fluctuations, perhaps B4?
            E1 = avgSqValue/totalSampleCount;
        }
        E0 /= totalCount;
//        System.out.println("E0 (old,new): "+E0+" "+E00);
        E1 /= totalCount;

//        System.out.println("weights");
        double newT1 = 0;
        for (IntSet pv : allMyData.keySet()) {
            MyData amd = allMyData.get(pv);
            long c = amd.unscreenedCount;
            if (c==0) continue;
            double w = localWeight.get(pv)/maxw;
            amd.weight = w;
            newT1 += tRatio * c * w / totalCount;
        }
        double y1 = Math.sqrt(E1*t1/(E0*t0));
        // newT1 is the new fraction of time we would spend calculating cluster values
        // y is the optimal fraction of time.
        // if y<newT, then we scale everything down (happy)
        // if y>newT1, then we scale up, but some weights>1.
        //    we'll still calculate all of them, but can't actually cause it to visit more than 100%
        double x1 = y1==0 ? 0 : y1/newT1;
        newT1 = 0;
        long totalUnscreened = 0;
        double newE1 = 0;
        for (IntSet pv : allMyData.keySet()) {
            MyData amd = allMyData.get(pv);
            double w = amd.weight*x1;
            if (w > 1) w = 1;
            else if (amd.sampleCount < 2) w = 1;
//            if (i>=targetCluster.pointCount()) System.out.println(String.format("%2d %12d  %6.4f  %6.4f  %6.4f", i, (count[i]-screenedCount[i]), weight[i], w, tRatio * (count[i]-screenedCount[i]) * w / totalCount));
            amd.weight = w;
            long c = amd.unscreenedCount;
            newT1 += c * w;
            totalUnscreened += c;
            double s = amd.getVarDouble();
            if (s>0) newE1 += c*s/w;
        }
        newE1 /= totalCount;
        newT1 *= tRatio/totalCount;
//        System.out.println(E00+" "+newE1+" "+newT1);
        if (!quiet) System.out.print(String.format("var0 frac %8.5f   t0 frac %8.5f  ideal t0 frac %5.3f  new t0 frac %5.3f   measure frac %5.3f  difficulty %10.4e\n", E0/(E0+E1), t0/(t0+t1), 1/(1+y1), 1/(1+newT1), newT1*(totalCount/tRatio/totalUnscreened), Math.sqrt((E00+newE1)*(1+newT1))));
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

    public static class MyData {
        public long unscreenedCount, sampleCount;
        public double weight;
        public BigDecimal sum;
        public double dsum, dsum2;
        public BigDecimal sum2;
//        public boolean locked;
        public MyData() {
            sum = BDZERO;
            sum2 = BDZERO;
        }
        
        public double getAvgDouble() {
            return sum.doubleValue()/sampleCount;
        }
        
        public double getDoubleAvg() {
            return dsum/sampleCount;
        }
        
        public double getDoubleVar() {
            if (sampleCount < 1) return Double.NaN;
            double avg = getDoubleAvg();
            double avg2 = avg*avg;
            double var = dsum2/sampleCount - avg2;
            if (var < avg2*1e-7) var = 0;
            return var;
        }
        
        public double getVarDouble() {
            if (sampleCount < 1) return Double.NaN;
            double avg = getAvgDouble();
            double avg2 = avg*avg;
            double var = sum2.doubleValue()/sampleCount - avg2;
            if (var < avg2*1e-13) var = 0;
            return var;
        }

        public double getVar(MathContext mc) {
            if (sampleCount < 1) return Double.NaN;
            BigDecimal avg = sum.divide(new BigDecimal(sampleCount, mc), mc);
            BigDecimal avg2 = avg.multiply(avg, mc);
            BigDecimal varBD = sum2.divide(new BigDecimal(sampleCount, mc), mc).subtract(avg2, mc);
            double var = varBD.doubleValue();
            if (var < avg2.doubleValue()*1e-25) var = 0;
            return var;
        }

        
        public void addData(BigDecimal value, MathContext mc, IntSet pv, int iThread) {
/*            if (Math.abs((sum.doubleValue()-dsum)/(2*dsum)) > 1e-6) {
                System.out.println(sum.toString()+" "+dsum);
                throw new RuntimeException("wtf");
            } */
/*            try {
                FileWriter fw = new FileWriter("thread"+iThread+"_raw.dat", true);
                fw.write(Arrays.toString(pv.v).replaceAll(" ","")+" "+value.toString()+"\n");
                fw.close();
            }
            catch (IOException ex) {
                throw new RuntimeException(ex);
            } */
//            BigDecimal oldsum = sum;
//            double olddsum = dsum;
            sum = sum.add(value, mc);
            dsum += value.doubleValue();
/*            if (Math.abs((sum.doubleValue()-dsum)/(2*dsum)) > 1e-6) {
                System.out.println(oldsum.toString()+" + "+value.toString()+" = "+sum.toString());
                System.out.println(olddsum+" + "+value.doubleValue()+" = "+dsum);
                throw new RuntimeException("oopsie");
            } */
            dsum2 += value.doubleValue() * value.doubleValue();
            value = value.multiply(value, mc);
            sum2 = sum2.add(value, mc);
            sampleCount++;
//            locked = false;
        }
    }
}
