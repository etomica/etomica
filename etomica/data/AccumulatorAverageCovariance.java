package etomica.data;

import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;

/**
 * Class that does the work of AccumulatorAverageFixed and also calculates the
 * covariance of the incoming data values (the covariance of each value with
 * each other value, as well as itself (which is actually its standard
 * deviation)) as well as the covariance of the block averages of the different
 * data quantities.
 */
public class AccumulatorAverageCovariance extends AccumulatorAverageFixed {

    public AccumulatorAverageCovariance() {
    }

    public AccumulatorAverageCovariance(long blockSize) {
        super(blockSize);
    }

    public void addData(IData data) {
        if (data.isNaN())
            return;
        super.addData(data);

        double[] x = covSum.getData();
        int n = data.getLength();
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                x[i*n+j] += data.getValue(i)*data.getValue(j);
            }
        }
    }
    
    protected void doBlockSum() {
        // need to do this first since blockSum gets zero'd by super.doBlockSum()
        double[] x = blockCovSum.getData();
        int n = currentBlockSum.getLength();
        long blockSizeSq = blockSize * blockSize;
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                x[i*n+j] += currentBlockSum.getValue(i)*currentBlockSum.getValue(j) / blockSizeSq;
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
        
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                x[i*n+j] -= average.getValue(i)*average.getValue(j);
            }
        }

        if (count > 1) {
            blockCovariance.E(blockCovSum);
            blockCovariance.TE(1.0/count);
            x = blockCovariance.getData();
            long totalCount = count*blockSize;
            double countSq = totalCount*totalCount;
            for (int i=0; i<n; i++) {
                for (int j=0; j<n; j++) {
                    x[i*n+j] -= sum.getValue(i)*sum.getValue(j) / countSq;
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

    public IEtomicaDataInfo processDataInfo(IEtomicaDataInfo incomingDataInfo) {
        int n = incomingDataInfo.getLength();
        covSum = new DataDoubleArray(new int[]{n,n});
        blockCovSum = new DataDoubleArray(new int[]{n,n});
        covariance = new DataDoubleArray(new int[]{n,n});
        blockCovariance = new DataDoubleArray(new int[]{n,n});
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
        DataDoubleArray.DataInfoDoubleArray covarianceInfo = new DataDoubleArray.DataInfoDoubleArray(
                incomingLabel+" covariance", new CompoundDimension(
                        new Dimension[]{incomingDataInfo.getDimension()},new double[]{2}), new int[]{n,n});
        DataDoubleArray.DataInfoDoubleArray blockCovarianceInfo = new DataDoubleArray.DataInfoDoubleArray(
                incomingLabel+" blk covariance", new CompoundDimension(
                        new Dimension[]{incomingDataInfo.getDimension()},new double[]{2}), new int[]{n,n});
        IEtomicaDataInfo[] myInfo = new IEtomicaDataInfo[nSuper+2];
        for (int i=0; i<nSuper; i++) {
            myInfo[i] = ((DataInfoGroup)dataInfo).getSubDataInfo(i);
        }
        myInfo[nSuper] = covarianceInfo;
        myInfo[nSuper+1] = blockCovarianceInfo;
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

    private static final long serialVersionUID = 1L;
    protected DataDoubleArray covSum, blockCovSum;
    protected DataDoubleArray covariance, blockCovariance;
}
