/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.MathContext;

import etomica.action.IAction;
import etomica.util.random.IRandom;
import etomica.data.AccumulatorAverageBlockless;
import etomica.data.AccumulatorAverageBlockless.AllData;
import etomica.data.DataTag;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDoubleBDArray;
import etomica.data.types.DataDoubleBDArray.DataInfoDoubleBDArray;
import etomica.data.types.DataGroup;
import etomica.units.dimensions.Null;

/**
 * Measures value of clusters in a box and returns the values
 * divided by the sampling bias from the sampling cluster.
 */
public class MeterVirialBDBin implements IAction {

    protected final IRandom random;
    protected final ClusterWheatley targetCluster;
    private final DataDoubleBDArray data;
    private final MathContext mc;
    private final IEtomicaDataInfo dataInfo;
    private final DataTag tag;
    private BoxCluster box;
    protected final AccumulatorAverageBlockless[] accumulators;
    protected final long[] unscreenedCount;
    protected final double[] weight;
    protected double maxSampleValue;
    protected double tRatio;
    protected static final BigDecimal BDZERO = new BigDecimal(0);
    protected static final BigDecimal BDONE = new BigDecimal(1);
    protected long nextReweightStep = 100000L;
    protected long totalCount = 0;

    /**
	 * Constructor for MeterVirial.
	 */
	public MeterVirialBDBin(ClusterWheatley targetCluster, IRandom random) {
		this.targetCluster = targetCluster;
		this.random = random;
		mc = new MathContext(40);
        data = new DataDoubleBDArray(1, 40);
        dataInfo = new DataInfoDoubleBDArray("Cluster Value",Null.DIMENSION, new int[]{1}, 40);
        tag = new DataTag();
        dataInfo.addTag(tag);
        DataDoubleBDArray.addBDOne(BDONE);
        DataDoubleBDArray.addBDZero(BDZERO);
        File vpi = new File("vpi.dat");
        if (vpi.exists()) vpi.delete();
        int nPairs = targetCluster.pointCount()*(targetCluster.pointCount()-1)/2;
        unscreenedCount = new long[nPairs+1];
        weight = new double[nPairs+1];
        accumulators = new AccumulatorAverageBlockless[nPairs+1];
        for (int i = 0; i<accumulators.length; i++) {
            accumulators[i] = new AccumulatorAverageBlockless();
            accumulators[i].putDataInfo(dataInfo);
            accumulators[i].setPushInterval(100000000);
            weight[i] = 1;
        }
	}

	public void setWeight(double newWeight) {
	    for (int i=0; i<weight.length; i++) {
	        weight[i] = newWeight;
	    }
	    nextReweightStep = Long.MAX_VALUE;
	}

	public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public void setTRatio(double newTRatio) {
        tRatio = newTRatio;
    }

    public void actionPerformed() {
        
        totalCount++;
        boolean hasNonZeroValue = targetCluster.checkConfig(box);
        if (!hasNonZeroValue) {
            return;
        }
        int edgeCount = targetCluster.getEdgeCount();
        unscreenedCount[edgeCount]++;
        if (weight[edgeCount] < 1 && weight[edgeCount] < random.nextDouble()) {
            return;
        }
        BigDecimal[] x = data.getData();
        double v = targetCluster.calcValue(box);
        if (v != 0) {
            double pi = box.getSampleCluster().value(box);
            x[0] = new BigDecimal(v, mc).divide(new BigDecimal(pi, mc), mc);

//          try {
//             FileWriter fw = new FileWriter("vpi.dat", true);
//              fw.write(pi+" "+Math.abs(v)+" "+edgeCount+"\n");
//              fw.close();
//          }
//          catch (IOException e) {throw new RuntimeException(e);}

        }
        else {
            x[0] = BDZERO;
        }
        
        accumulators[edgeCount].addData(data);
        
        if (totalCount >= nextReweightStep) {
            recomputeWeights();
            nextReweightStep = totalCount*2;
        }
    }
    
    public void recomputeWeights() {
        // tRatio is the ratio of the time needed to compute the biconnected
        // value (and reference value) for one configuration to the time needed
        // to generate a configuration and screen it (and any other overhead).
        long totalSampleCount = 0;
        double totalSqValue = 0;
        for (int i=0; i<accumulators.length; i++) {
            AccumulatorAverageBlockless ai = accumulators[i];
            long sc = ai.getSampleCount();
            totalSampleCount += sc;

            DataGroup allYourBase = (DataGroup)ai.getData();
            double average = ((DataDoubleBDArray)allYourBase.getData(ai.AVERAGE.index)).getData()[0].doubleValue();
            double stdev = ((DataDoubleBDArray)allYourBase.getData(ai.STANDARD_DEVIATION.index)).getData()[0].doubleValue();
            totalSqValue += sc * (stdev*stdev + average*average);
        }
        double avgSqValue = totalSqValue / totalSampleCount;
        double t0 = totalCount;
        double t1 = totalSampleCount*tRatio;
        double E0 = 0;
        double E1 = 0;
        // E0 = sum(sci*(steps-sci)/steps * ai^2)
        // E1 = sum(sci*sci*stdev*stdev/sampci)

        double maxw = 0;
        double[] localWeight = new double[unscreenedCount.length];
        for (int i=0; i<unscreenedCount.length; i++) {
            long c = unscreenedCount[i];
            long sampleCount = accumulators[i].getSampleCount();
            if (sampleCount == 0) {
                // we have never seen i bonds, or the configuration was always screened
                localWeight[i] = avgSqValue;
                continue;
            }
            localWeight[i] = avgSqValue/sampleCount;
            DataGroup allYourBase = (DataGroup)accumulators[i].getData();
            DataDoubleBDArray averageData = (DataDoubleBDArray)allYourBase.getData(accumulators[i].AVERAGE.index);
            DataDoubleBDArray stdevData = (DataDoubleBDArray)allYourBase.getData(accumulators[i].STANDARD_DEVIATION.index);

            BigDecimal avg = averageData.getData()[0];
            if (avg.doubleValue() == 0) {
                // we have seen i bonds, but the cluster value is zero
                continue;
            }

            // E0 = sum(sci*(steps-sci)/steps * ai^2)
            E0 += c*((double)(totalCount-c))/totalCount * avg.doubleValue()*avg.doubleValue();

            if (sampleCount<2) {
                continue;
            }
            
            BigDecimal stdev = stdevData.getData()[0];

            if (sampleCount > 1) {
                BigDecimal lw = stdev.multiply(stdev, mc);
                localWeight[i] += lw.doubleValue();
                if (sampleCount > 10 && localWeight[i] > maxw) maxw = localWeight[i];
            }

            // E1 = sum(sci*sci*stdev*stdev/sampci)
            E1 += c*((double)c)/accumulators[i].getSampleCount() * stdev.doubleValue()*stdev.doubleValue();
        }
        if (E0 == 0) {
            nextReweightStep *= 2;
            return;
        }
        if (E1 == 0) {
            // no value fluctuations, perhaps B4?
            E1 = avgSqValue/totalSampleCount;
        }
        E0 /= (totalCount*totalCount);
        E1 /= (totalCount*totalCount);
//        System.out.println("weights");
        double newT1 = 0;
        for (int i=0; i<unscreenedCount.length; i++) {
            weight[i] = localWeight[i]/maxw;
            newT1 += tRatio * unscreenedCount[i] * weight[i] / totalCount;
        }
        double y1 = Math.sqrt(E1*t1/(E0*t0));
        // newT1 is the new fraction of time we would spend calculating cluster values
        // y is the optimal fraction of time.
        // if y<newT, then we scale everything down (happy)
        // if y>newT1, then we scale up, but some weights>1.
        //    we'll still calculate all of them, but can't actually cause it to visit more than 100%
        double x1 = y1 / newT1;
        newT1 = 0;
        long totalUnscreened = 0;
        for (int i=0; i<unscreenedCount.length; i++) {
            double w = weight[i]*x1;
            if (w > 1) w = 1;
            else if (accumulators[i].getSampleCount() < 2) w = 1;
            long c = unscreenedCount[i];
//            if (i>=targetCluster.pointCount()) System.out.println(String.format("%2d %12d  %6.4f  %6.4f  %6.4f", i, (count[i]-screenedCount[i]), weight[i], w, tRatio * (count[i]-screenedCount[i]) * w / totalCount));
            weight[i] = w;
            newT1 += tRatio * c * w / totalCount;
            totalUnscreened += c;
        }
        System.out.print(String.format("var0 frac %8.5f   t0 frac %8.5f  ideal t0 frac %5.3f  new t0 frac %5.3f   measure frac %5.3f\n", E0/(E0+E1), t0/(t0+t1), 1/(1+y1), 1/(1+newT1), newT1*(totalCount/tRatio/totalUnscreened)));
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
    
    public long[] getUnscreenedCount() {
        return unscreenedCount;
    }
    
    public AccumulatorAverageBlockless[] getAccumulators() {
        return accumulators;
    }
    
    public MathContext getMathContext() {
        return mc;
    }
    
    public void writeData(String filename) {
        try {
            FileWriter fw = new FileWriter(filename);
            fw.write(""+totalCount+"\n");
            for (int i=0; i<unscreenedCount.length; i++) {
                long sampleCount = accumulators[i].getSampleCount();
                if (sampleCount == 0) continue;
                long usc = unscreenedCount[i];
                AllData ad = accumulators[i].getRawData();
                fw.write(i+" "+usc+" "+ad.count+" "+((DataDoubleBDArray)ad.sum).getData()[0].toString()+" "+((DataDoubleBDArray)ad.sumSquare).getData()[0].toString()+"\n");
            }
            fw.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public void writeWeights(String filename) {
        try {
            FileWriter fw = new FileWriter(filename);
            for (int i=0; i<unscreenedCount.length; i++) {
                long sampleCount = accumulators[i].getSampleCount();
                if (sampleCount == 0) continue;
                fw.write(i+" "+weight[i]+"\n");
            }
            fw.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public void readData(String[] filenames) {
        int nPairs = targetCluster.pointCount()*(targetCluster.pointCount()-1)/2;
        BigDecimal[] sums = new BigDecimal[nPairs+1];
        BigDecimal[] sumSquares = new BigDecimal[nPairs+1];
        long[] sampleCounts = new long[nPairs+1];
        try {
            for (String filename : filenames) {
                File f = new File(filename);
                if (!f.exists()) continue;
                FileReader fr = new FileReader(filename);
                BufferedReader bufReader = new BufferedReader(fr);
                String line = bufReader.readLine();
                totalCount += Long.parseLong(line);
                while ((line=bufReader.readLine()) != null) {
                    String iStr = line.replaceAll(" .*", "");
                    int i = Integer.parseInt(iStr);
                    String[] values = line.replaceAll("^[^ ]* ", "").split(" +");
                    long usc = Long.parseLong(values[0]);
                    long sampleCount = Long.parseLong(values[1]);
                    BigDecimal sum = new BigDecimal(values[2], mc);
                    BigDecimal sumSquare = new BigDecimal(values[3], mc);
                    if (unscreenedCount[i] > 0) {
                        unscreenedCount[i] += usc;
                        sampleCounts[i] += sampleCount;
                        sums[i] = sums[i].add(sum, mc);
                        sumSquares[i] = sumSquares[i].add(sumSquare, mc);
                    }
                    else {
                        unscreenedCount[i] = usc;
                        sampleCounts[i] = sampleCount;
                        weight[i] = 1.0;
                        sums[i] = sum;
                        sumSquares[i] = sumSquare;
                    }
                }
                bufReader.close();
            }
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
        for (int i=0; i<unscreenedCount.length; i++) {
            long usc = unscreenedCount[i];
            if (usc == 0) continue;
            DataDoubleBDArray sumData = new DataDoubleBDArray(1, 40);
            DataDoubleBDArray sumSquareData = new DataDoubleBDArray(1, 40);
            sumData.E(new BigDecimal[]{sums[i]});
            sumSquareData.E(new BigDecimal[]{sumSquares[i]});
            AccumulatorAverageBlockless acc = new AccumulatorAverageBlockless();
            acc.putDataInfo(dataInfo);
            acc.setPushInterval(100000000);
            accumulators[i] = acc;
            acc.setRawData(new AllData(sampleCounts[i], sumData, sumSquareData));
        }
        recomputeWeights();
    }

    public void readWeights(String filename) {
        File f = new File(filename);
        if (!f.exists()) return;
        try {
            FileReader fr = new FileReader(filename);
            BufferedReader bufReader = new BufferedReader(fr);
            String line = null;
            while ((line=bufReader.readLine()) != null) {
                String iStr = line.replaceAll(" .*", "");
                int i = Integer.parseInt(iStr);
                String weightStr = line.replaceAll("^[^ ]* ", "");
                weight[i] = Double.parseDouble(weightStr);
                unscreenedCount[i] = 0;
                AccumulatorAverageBlockless acc = new AccumulatorAverageBlockless();
                acc.putDataInfo(dataInfo);
                acc.setPushInterval(100000000);
                accumulators[i] = acc;
            }
            bufReader.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
        nextReweightStep = Long.MAX_VALUE;
    }

}
