/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;

/**
 * Class that does the work of AccumulatorAverageFixed and also calculates the
 * covariance of the incoming data values (the covariance of each value with
 * each other value, as well as itself (which is actually its standard
 * deviation)) as well as the covariance of the block averages of the different
 * data quantities.
 */
public class AccumulatorAverageCovariance extends AccumulatorAverageFixed {

    public AccumulatorAverageCovariance() {
        this(true);
    }

    /**
     * @param fullCovariance enables calculation of covariance for each pair of
     * data values.  With fullCovariance off, covariance is computed between
     * the first data value and all others.
     */
    public AccumulatorAverageCovariance(boolean fullCovariance) {
        this.fullCovariance = fullCovariance;
    }

    public AccumulatorAverageCovariance(long blockSize) {
        this(blockSize, true);
    }

    public AccumulatorAverageCovariance(long blockSize, boolean fullCovariance) {
        super(blockSize);
        this.fullCovariance = fullCovariance;
    }

    public boolean addData(IData data) {
        if (data.isNaN())
            return false;
        super.addData(data);

        double[] x = covSum.getData();
        int n = data.getLength();
        if (fullCovariance) {
            for (int i=0; i<n; i++) {
                double ix = data.getValue(i);
                x[i*n+i] += ix*ix;
                for (int j=i+1; j<n; j++) {
                    double ijx = ix*data.getValue(j);
                    x[i*n+j] += ijx;
                    x[j*n+i] += ijx;
                }
            }
        }
        else {
            double x0 = data.getValue(0);
            for (int j=0; j<n; j++) {
                double j0x = x0*data.getValue(j);
                x[j] += j0x;
            }
        }
        return true;
    }
    
    protected void doBlockSum() {
        // need to do this first since blockSum gets zero'd by super.doBlockSum()
        double[] x = blockCovSum.getData();
        int n = currentBlockAvg.getLength();

        if (fullCovariance) {
            for (int i=0; i<n; i++) {
                double ix = currentBlockAvg.getValue(i);
                x[i * n + i] += ix * ix;
                for (int j=i+1; j<n; j++) {
                    double ijx = ix * currentBlockAvg.getValue(j);
                    x[i*n+j] += ijx;
                    x[j*n+i] += ijx;
                }
            }
        }
        else {
            double x0 = currentBlockAvg.getValue(0);
            for (int j=0; j<n; j++) {
                double j0x = x0 * currentBlockAvg.getValue(j);
                x[j] += j0x;
            }
        }

        super.doBlockSum();
    }

    public IData getData() {
        if (average == null)
            return null;
        super.getData();

        int n = average.getLength();
        covariance.E(covSum);
        long nTotalData = count*blockSize + (blockCountDown-blockSize);
        covariance.TE(1.0/nTotalData);
        double[] x = covariance.getData();

        if (fullCovariance) {
            for (int i=0; i<n; i++) {
                double ix = average.getValue(i);
                x[i*n+i] -= ix*ix;
                for (int j=i+1; j<n; j++) {
                    double ijx = ix*average.getValue(j);
                    x[i*n+j] -= ijx;
                    x[j*n+i] -= ijx;
                }
            }
        }
        else {
            double x0 = average.getValue(0);
            for (int j=0; j<n; j++) {
                double j0x = x0*average.getValue(j);
                x[j] += j0x;
            }
        }

        if (count > 1) {
            blockCovariance.E(blockCovSum);
            blockCovariance.TE(1.0/count);
            x = blockCovariance.getData();
            if (fullCovariance) {
                for (int i=0; i<n; i++) {
                    double ix = average.getValue(i);
                    x[i * n + i] -= ix * ix;
                    for (int j=i+1; j<n; j++) {
                        double ijx = ix * average.getValue(j);
                        x[i*n+j] -= ijx;
                        x[j*n+i] -= ijx;
                    }
                }
            }
            else {
                double x0 = average.getValue(0);
                for (int j=0; j<n; j++) {
                    double j0x = x0 * average.getValue(j);
                    x[j] += j0x;
                }
            }
        }
        else {
            blockCovariance.E(Double.NaN);
        }
        return dataGroup;
    }
    
    public void reset() {
        super.reset();
        if (covSum == null) {
            return;
        }
        covSum.E(0);
        blockCovSum.E(0);
    }

    public IDataInfo processDataInfo(IDataInfo incomingDataInfo) {
        int n = incomingDataInfo.getLength();
        if (fullCovariance) {
            covSum = new DataDoubleArray(new int[]{n,n});
            blockCovSum = new DataDoubleArray(new int[]{n,n});
            covariance = new DataDoubleArray(new int[]{n,n});
            blockCovariance = new DataDoubleArray(new int[]{n,n});
        }
        else {
            covSum = new DataDoubleArray(new int[]{n});
            blockCovSum = new DataDoubleArray(new int[]{n});
            covariance = new DataDoubleArray(new int[]{n});
            blockCovariance = new DataDoubleArray(new int[]{n});
        }
        super.processDataInfo(incomingDataInfo);

        int nSuper = dataGroup.getNData();
        IData[] myData = new IData[nSuper+2];
        for (int i=0; i<nSuper; i++) {
            myData[i] = dataGroup.getData(i);
        }
        myData[nSuper] = covariance;
        myData[nSuper+1] = blockCovariance;
        dataGroup = new DataGroup(myData);

        String incomingLabel = incomingDataInfo.getLabel();
        IDataInfo[] myInfo = new IDataInfo[nSuper+2];
        for (int i=0; i<nSuper; i++) {
            myInfo[i] = ((DataInfoGroup)dataInfo).getSubDataInfo(i);
        }
        if (fullCovariance) {
            myInfo[nSuper] = new DataDoubleArray.DataInfoDoubleArray(
                incomingLabel+" covariance", new CompoundDimension(
                        new Dimension[]{incomingDataInfo.getDimension()},new double[]{2}), new int[]{n,n});
            myInfo[nSuper+1] = new DataDoubleArray.DataInfoDoubleArray(
                incomingLabel+" blk covariance", new CompoundDimension(
                        new Dimension[]{incomingDataInfo.getDimension()},new double[]{2}), new int[]{n,n});
        }
        else {
            myInfo[nSuper] = new DataDoubleArray.DataInfoDoubleArray(
                    incomingLabel+" covariance", new CompoundDimension(
                            new Dimension[]{incomingDataInfo.getDimension()},new double[]{2}), new int[]{n});
                myInfo[nSuper+1] = new DataDoubleArray.DataInfoDoubleArray(
                    incomingLabel+" blk covariance", new CompoundDimension(
                            new Dimension[]{incomingDataInfo.getDimension()},new double[]{2}), new int[]{n});
        }
        dataInfo = new DataInfoGroup(incomingLabel, incomingDataInfo.getDimension(), myInfo);
        dataInfo.addTag(tag);
        return dataInfo;
    }

    public static StatType[] statChoices() {
        return new StatType[] {MOST_RECENT,AVERAGE,ERROR,STANDARD_DEVIATION,BLOCK_CORRELATION,
                                                  COVARIANCE,BLOCK_COVARIANCE};
    }
    
    public static final StatType COVARIANCE = new StatType("Covariance", 5);
    public static final StatType BLOCK_COVARIANCE = new StatType("Block covariance", 6);

    protected DataDoubleArray covSum, blockCovSum;
    protected DataDoubleArray covariance, blockCovariance;
    protected final boolean fullCovariance;
}
