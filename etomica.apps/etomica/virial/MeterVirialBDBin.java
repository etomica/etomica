/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import java.io.File;
import java.math.BigDecimal;
import java.math.MathContext;

import etomica.action.IAction;
import etomica.api.IRandom;
import etomica.data.AccumulatorAverageBlockless;
import etomica.data.DataTag;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDoubleBDArray;
import etomica.data.types.DataDoubleBDArray.DataInfoDoubleBDArray;
import etomica.data.types.DataGroup;
import etomica.units.Null;

/**
 * Measures value of clusters in a box and returns the values
 * divided by the sampling bias from the sampling cluster.
 */
public class MeterVirialBDBin implements IAction {

    protected static final BigDecimal BDZERO = new BigDecimal(0);
    protected static final BigDecimal BDONE = new BigDecimal(1);
    protected long nextReweightStep = 10000;
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
        count = new long[nPairs+1];
        screenedCount = new long[nPairs+1];
        weight = new double[nPairs+1];
        accumulators = new AccumulatorAverageBlockless[nPairs+1];
        for (int i = 0; i<accumulators.length; i++) {
            accumulators[i] = new AccumulatorAverageBlockless();
            accumulators[i].putDataInfo(dataInfo);
            accumulators[i].setPushInterval(100000000);
            weight[i] = 1;
        }
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
        BigDecimal[] x = data.getData();
        boolean hasNonZeroValue = targetCluster.checkConfig(box);
        int edgeCount = targetCluster.getEdgeCount();
        count[edgeCount]++;
        if (!hasNonZeroValue) {
            screenedCount[edgeCount]++;
            return;
        }
        if (weight[edgeCount] < 1 && weight[edgeCount] < random.nextDouble()) {
            return;
        }
        double v = targetCluster.calcValue(box);
        if (v != 0) {
            double pi = box.getSampleCluster().value(box);
            if (Math.abs(v/pi) > maxSampleValue) maxSampleValue = Math.abs(v/pi);
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
        
        if (totalCount >= nextReweightStep && maxSampleValue > 0) {
            // tRatio is the ratio of the time needed to compute the biconnected
            // value (and reference value) for one configuration to the time needed
            // to generate a configuration and screen it (and any other overhead).
            long totalSampleCount = 0;
            for (int i=0; i<count.length; i++) {
                totalSampleCount += accumulators[i].getSampleCount();
            }
            double t0 = totalCount;
            double t1 = totalSampleCount*tRatio;
            double E0 = 0;
            double E1 = 0;
            // E0 = sum(sci*(steps-sci)/steps * ai^2)
            // E1 = sum(sci*sci*stdev*stdev/sampci)

            double maxw = 0;
            double[] localWeight = new double[count.length];
            double maxAvg = 0;
            for (int i=0; i<count.length; i++) {
                long c = count[i] - screenedCount[i];
                long sampleCount = accumulators[i].getSampleCount();
                if (sampleCount == 0) sampleCount = 1;
                localWeight[i] = maxSampleValue*maxSampleValue/sampleCount;
                if (sampleCount==1) {
                    // we have never seen i bonds, or the configuration was always screened
                    // or we just have no statistics
                    continue;
                }
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

                if (Math.abs(avg.doubleValue()) > maxAvg) maxAvg = Math.abs(avg.doubleValue());
                BigDecimal stdev = stdevData.getData()[0];

                if (sampleCount > 1) {
                    BigDecimal lw = stdev.multiply(stdev, mc);
                    localWeight[i] += lw.doubleValue();
                    if (sampleCount > 10 && localWeight[i] > maxw) maxw = localWeight[i];
                }

                // E1 = sum(sci*sci*stdev*stdev/sampci)
                E1 += c*((double)c)/accumulators[i].getSampleCount() * stdev.doubleValue()*stdev.doubleValue();
            }
            if (E0 == 0 || E1 == 0) {
                nextReweightStep *= 2;
                return;
            }
            E0 /= (totalCount*totalCount);
            E1 /= (totalCount*totalCount);
            if (E1 == 0) {
                // no value fluctuations, perhaps B4?
                E1 = 0.1*E0;
            }
//            System.out.println("weights");
            double newT1 = 0;
            for (int i=0; i<count.length; i++) {
                if (localWeight[i] > maxw) {
                    weight[i] = 1;
                }
                else {
                    weight[i] = localWeight[i] / maxw;
                }
                newT1 += tRatio * (count[i]-screenedCount[i]) * weight[i] / totalCount;
            }
            double y = 1.0/(1.0 + Math.sqrt(E0*t0/(E1*t1)));
            System.out.print(String.format("var0 frac %8.5f  t0 frac %8.5f   y %5.3f   new T0 frac %5.3f\n", E0/(E0+E1), t0/(t0+t1), 1-y, 1/(1+newT1)));
            // newT1 is the new fraction of time we would spend calculating cluster values
            // y is the optimal fraction of time.
            // if y<newT, then we scale everything down (happy)
            // if y>newT1, then we scale up, but some weights>1.
            //    we'll still calculate all of them, but can't actually cause it to visit more than 100%
            y /= newT1/(newT1+1);
            for (int i=0; i<count.length; i++) {
                double w = weight[i]*y;
                if (w > 1) w = 1;
                else if (accumulators[i].getSampleCount() < 10) w = 1;
                if (i>=targetCluster.pointCount()) System.out.println(String.format("%2d %12d  %6.4f  %6.4f  %6.4f", i, (count[i]-screenedCount[i]), weight[i], w, tRatio * (count[i]-screenedCount[i]) * w / totalCount));
                weight[i] = w;
            }
            nextReweightStep = totalCount*2;
        }
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
    
    public long[] getCount() {
        return count;
    }
    
    public long[] getScreenedCount() {
        return screenedCount;
    }
    
    public AccumulatorAverageBlockless[] getAccumulators() {
        return accumulators;
    }
    
    public MathContext getMathContext() {
        return mc;
    }

    protected final IRandom random;
    protected final ClusterWheatley targetCluster;
	private final DataDoubleBDArray data;
	private final MathContext mc;
	private final IEtomicaDataInfo dataInfo;
    private final DataTag tag;
    private BoxCluster box;
    protected final AccumulatorAverageBlockless[] accumulators;
    protected final long[] count, screenedCount;
    protected final double[] weight;
    protected double maxSampleValue;
    protected double tRatio;
}
