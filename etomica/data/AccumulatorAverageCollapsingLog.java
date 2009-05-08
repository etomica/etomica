package etomica.data;

import java.util.ArrayList;

import etomica.api.IData;
import etomica.api.IRandom;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.Null;
import etomica.units.Quantity;
import etomica.util.Arrays;
import etomica.util.RandomNumberGenerator;


/**
 * AccumulatorAverage that keeps track of block averages and standard
 * deviations corresponding to all block sizes from 1 to the n where n is the
 * largest power of 2 less than the current number of samples.  In addition to
 * block averages and standard deviations of the sampled quantity, the natural
 * log of the block averages and standard deviation of the natural log is also
 * calculated (this is what makes this class actually useful).  getData returns
 * the average natural log values for all block sizes as a DataFunction.
 * <p>
 * This class also uses a method to improve the block statistics for large.
 * All block averages are kept for a particular block size (n1 blocks of size
 * b1).  Statistics for that block are computed from that data.  Block
 * averages for blocks of the next larger size (2*n1) are computed by randomly
 * selecting n1/2 pairs of block values from the b1 set (without including any
 * block value more than once).  After using all block values, n1/2 additional
 * pairs are taken from the full b1 set, such that n1 blocks of size b1*2 are
 * obtained, instead of n1/2 blocks that would normally exist.  Statistics for
 * block size b1*2 are computed from these n1 blocks.  This process is
 * repeated for each block size such that the largest block size (which would
 * normally have 1 block) has n1 blocks.
 * <p>
 * This accumulator can only operate on Data with a single value.
 */
public class AccumulatorAverageCollapsingLog extends DataAccumulator implements DataSourceIndependent {

    public AccumulatorAverageCollapsingLog(IRandom random) {
        this.random = random;
        setPushInterval(100);
        setNumRawDataDoubles(14);
    }

    /**
     * Sets the number of times the random block generation method is used.
     * The number of blocks used will be 2^numRawDataDoubles.
     */
    public void setNumRawDataDoubles(int newNumRawDataDoubles) {
        nRawDataDoubles = newNumRawDataDoubles;
        intArrayList = new ArrayList<Integer>(1<<nRawDataDoubles);
        reset();
    }

    /**
     * Checks that incoming Data implements Data, and returns null if
     * this is so. Otherwise throws a ClassCastException, as there is no data
     * caster to Data.
     */
    public DataPipe getDataCaster(IEtomicaDataInfo incomingDataInfo) {
        if (incomingDataInfo.getLength() > 1) {
            throw new RuntimeException("AccumulatorAverageCollapsingLog can only handle single data");
        }
        return null;
    }

    public IData processData(IData inputData) {
        if(!active) return null;
        int oldSize = sums.length;
        addData(inputData);
        if (sums.length > oldSize) {
            return getData();
        }
        return null;
    }

    /**
     * Add the given values to the sums and block sums. If any of the given data
     * values is NaN, method returns with no effect on accumulation sums.
     */
    public void addData(IData data) {
        if (data.isNaN())
            return;
        double value = data.getValue(0);
        if (rawDataBlockSize == 0) {
            // special case early part where we keep every sample as raw data
            rawData[nRawData] = value;
            nRawData++;
            count++;
            int nDoubles = 63-Long.numberOfLeadingZeros(count);
            if (nDoubles >= averages.getLength()) {
                resizeData(nDoubles+1);
            }
            if (nRawData == 1<<nRawDataDoubles) {
                // we filled our raw data.  do summation for blockSize=0,
                // prepare for blockSize=1
                for (int i=0; i<nRawData; i+=2) {
                    doSums(sums[0], rawData[i]);
                    doSums(lSums[0], Math.log(rawData[i]));
                    doSums(sums[0], rawData[i+1]);
                    doSums(lSums[0], Math.log(rawData[i+1]));
                    rawData[i/2] = 0.5*(rawData[i]+rawData[i+1]);
                }
                for (int i=nRawData/2; i<nRawData; i++) {
                    rawData[i] = 0;
                }
                rawDataBlockSize++;
                nRawData /= 2;
            }
            return;
        }
        // process sample as blockSize=0 data
        doSums(sums[0], value);
        doSums(lSums[0], Math.log(value));
        blockSums[1] += value;
        count++;
        int nDoubles = Long.numberOfTrailingZeros(count);
        for (int i=1; i<nDoubles+1; i++) {
            blockSums[i] *= 0.5;
            if (rawDataBlockSize == i) {
                // we consider this block average to be raw data
                rawData[nRawData] = blockSums[i];
                blockSums[i] = 0;
                nRawData++;
                if (nRawData == 1<<nRawDataDoubles) {
                    // we filled our raw data
                    resizeData(nDoubles+1);
                    for (int j=0; j<nRawData; j+=2) {
                        doSums(sums[i], rawData[j]);
                        doSums(lSums[i], Math.log(rawData[j]));
                        doSums(sums[i], rawData[j+1]);
                        doSums(lSums[i], Math.log(rawData[j+1]));
                        rawData[j/2] = 0.5*(rawData[j]+rawData[j+1]);
                    }
                    for (int j=nRawData/2; j<nRawData; j++) {
                        rawData[j] = 0;
                    }
                    rawDataBlockSize++;
                    nRawData /= 2;
                }
                return;
            }
            // process block average normally
            doSums(sums[i], blockSums[i]);
            doSums(lSums[i], Math.log(blockSums[i]));
            blockSums[i+1] += blockSums[i];
            blockSums[i] = 0;
        }
    }
    
    protected void resizeData(int newSize) {
        blockSums = Arrays.resizeArray(blockSums, newSize+1);
        sums = (double[][])Arrays.resizeArray(sums, newSize);
        sums[newSize-1] = new double[3];
        lSums = (double[][])Arrays.resizeArray(lSums, newSize);
        lSums[newSize-1] = new double[3];
        averages = new DataFunction(new int[]{newSize});
        stdev = new DataFunction(new int[]{newSize});
        avgDataInfo = new DataInfoFunction(avgDataInfo.getLabel(), avgDataInfo.getDimension(), this);
        nDataInfo = new DataInfoDoubleArray("block size", Quantity.DIMENSION, new int[]{newSize});
        nData = new DataDoubleArray(newSize);
        double[] n = nData.getData();
        for (int i=0; i<n.length; i++) {
            n[i] = 1L<<i;
        }
        
        if (dataSink != null) {
            dataSink.putDataInfo(avgDataInfo);
        }
    }

    /*
     * Returns the average logs
     */
    public IData getData() {
        return getAverageLogs();
    }

    /**
     * Returns the averages for all block sizes
     */
    public IData getAverages() {
        double[] av = averages.getData();
        for (int i=0; i<rawDataBlockSize; i++) {
            av[i] = sums[i][1] / (count/(1L<<i));
        }
        int n = averages.getLength();
        double iSum = 0;
        for (int j=0; j<nRawData; j++) {
            iSum += rawData[j];
        }
        for (int i=rawDataBlockSize; i<n; i++) {
            av[i] = iSum / nRawData;
        }
        return averages;
    }

    /**
     * Returns the average natural logs of the block averages, for all block
     * sizes.
     */
    public IData getAverageLogs() {
        double[] av = averages.getData();
        for (int i=0; i<rawDataBlockSize; i++) {
            av[i] = lSums[i][1] / (count/(1L<<i));
        }
        double lSum = 0;
        for (int j=0; j<nRawData; j++) {
            lSum += Math.log(rawData[j]);
            rawData2[j] = rawData[j];
        }
        av[rawDataBlockSize] = lSum/nRawData;
        for (int i=rawDataBlockSize+1; i<rawDataBlockSize+(32-Integer.numberOfLeadingZeros(nRawData)); i++) {
            lSum = 0;
            reblockData();
            for (int j=0; j<nRawData; j++) {
                lSum += Math.log(rawData2[j]);
            }
            av[i] = lSum/nRawData;
        }
        return averages;
    }

    /**
     * Reblock data from rawData2 into rawData3, then copy it back to rawData2 
     */
    protected void reblockData() {
        intArrayList.clear();
        for (int j=0; j<nRawData; j++) {
            ensureHappyIntArrayList();
            int idxidx1 = random.nextInt(intArrayList.size());
            int idx1 = intArrayList.get(idxidx1);
            intArrayList.remove(idxidx1);
            ensureHappyIntArrayList();
            int idxidx2 = random.nextInt(intArrayList.size());
            int idx2 = intArrayList.get(idxidx2);
            intArrayList.remove(idxidx2);
            double sum12 = (rawData2[idx1] + rawData2[idx2])/2;
            rawData3[j] = sum12;
        }
        for (int j=0; j<nRawData; j++) {
            rawData2[j] = rawData3[j];
        }
    }

    /**
     * Refills inArrayList if it's empty.
     */
    protected void ensureHappyIntArrayList() {
        if (intArrayList.size() == 0) {
            for (int j=0; j<nRawData; j++) {
                intArrayList.add(j);
            }
        }
    }
    
    /**
     * Returns the standard deviation of the block averages, for all block
     * sizes.
     */
    public IData getStdev() {
        double[] std = stdev.getData();
        for (int i=0; i<rawDataBlockSize; i++) {
            double av = sums[i][1] / (count/(1L<<i));
            std[i] = Math.sqrt(sums[i][2] / (count/(1L<<i)) - av*av);
        }
        double iSum = 0, iSum2 = 0;
        for (int j=0; j<nRawData; j++) {
            iSum += rawData[j];
            iSum2 += rawData[j]*rawData[j];
            rawData2[j] = rawData[j];
        }
        std[rawDataBlockSize] = Math.sqrt(iSum2/nRawData - (iSum*iSum)/(nRawData*nRawData));
        for (int i=rawDataBlockSize+1; i<rawDataBlockSize+(32-Integer.numberOfLeadingZeros(nRawData)); i++) {
            iSum = 0;
            iSum2 = 0;
            reblockData();
            for (int j=0; j<nRawData; j++) {
                iSum += rawData2[j];
                iSum2 += rawData2[j]*rawData2[j];
            }
            std[i] = Math.sqrt(iSum2/nRawData - (iSum*iSum)/(nRawData*nRawData));
        }
        return stdev;
    }
    
    /**
     * Returns the standard deviation of the log of the block averages, for all
     * block sizes.
     */
    public IData getStdevLog() {
        double[] std = stdev.getData();
        for (int i=0; i<rawDataBlockSize; i++) {
            double av = lSums[i][1] / (count/(1L<<i));
            std[i] = Math.sqrt(lSums[i][2] / (count/(1L<<i)) - av*av);
        }
        double iSum = 0, iSum2 = 0;
        for (int j=0; j<nRawData; j++) {
            double l = Math.log(rawData[j]);
            iSum += l;
            iSum2 += l*l;
            rawData2[j] = rawData[j];
        }
        std[rawDataBlockSize] = Math.sqrt(iSum2/nRawData - (iSum*iSum)/(nRawData*nRawData));
        for (int i=rawDataBlockSize+1; i<rawDataBlockSize+(32-Integer.numberOfLeadingZeros(nRawData)); i++) {
            iSum = 0;
            iSum2 = 0;
            reblockData();
            for (int j=0; j<nRawData; j++) {
                double l = Math.log(rawData2[j]);
                iSum += l;
                iSum2 += l*l;
            }
            std[i] = Math.sqrt(iSum2/nRawData - (iSum*iSum)/(nRawData*nRawData));
        }
        return stdev;
    }
    /**
     * Performs the block sum after
     */
    protected void doSums(double[] subSums, double v) {
        double vpj = 1;
        for (int i=1; i<3; i++) {
            vpj *= v;
            subSums[i] += vpj;
        }
    }
    /**
     * Resets all sums to zero. All statistics are cleared.
     */
    public void reset() {
        count = 0;
        blockSums = new double[1];
        sums = new double[0][3];
        lSums = new double[0][3];
        averages = new DataFunction(new int[]{0});
        stdev = new DataFunction(new int[]{0});
        if (avgDataInfo != null) {
            avgDataInfo = new DataInfoFunction(avgDataInfo.getLabel(), avgDataInfo.getDimension(), this);
        }
        nDataInfo = new DataInfoDoubleArray("block size", Quantity.DIMENSION, new int[]{0});
        rawData = new double[1<<nRawDataDoubles];
        rawData2 = new double[1<<nRawDataDoubles];
        rawData3 = new double[1<<nRawDataDoubles];
        rawDataBlockSize = 0;
        nRawData = 0;
    }

    public IEtomicaDataInfo processDataInfo(IEtomicaDataInfo incomingDataInfo) {
        avgDataInfo = new DataInfoFunction(incomingDataInfo.getLabel(), incomingDataInfo.getDimension(), this);
        reset();
        return avgDataInfo;
    }
    
    public int getIndependentArrayDimension() {
        return 1;
    }

    public DataDoubleArray getIndependentData(int i) {
        return nData;
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return nDataInfo;
    }

    private static final long serialVersionUID = 1L;
    protected long count;
    protected double[] blockSums;
    protected double[][] sums, lSums;
    protected DataFunction averages, stdev;
    protected DataInfoFunction avgDataInfo;
    protected DataInfoDoubleArray nDataInfo;
    protected DataDoubleArray nData;
    protected double[] rawData, rawData2, rawData3;
    protected int nRawData, rawDataBlockSize;
    protected int nRawDataDoubles;
    protected ArrayList<Integer> intArrayList;
    protected final IRandom random;
    
    public static void main(String[] args) {
        AccumulatorAverageCollapsingLog ac = new AccumulatorAverageCollapsingLog(new RandomNumberGenerator());
        ac.putDataInfo(new DataInfoDouble("foo", Null.DIMENSION));
        DataDouble d = new DataDouble();
        RandomNumberGenerator rng = new RandomNumberGenerator();
        for (long i=0; i<1L<<25; i++) {
            d.x = Math.exp(rng.nextGaussian());
            ac.putData(d);
        }
        System.out.println(ac.getAverages().getValue(0)+" "+ac.getAverages().getValue(1)+" "+ac.getAverages().getValue(2)+" "+ac.getAverages().getValue(3));
        System.out.println(ac.getStdev().getValue(0)+" "+ac.getStdev().getValue(1)+" "+ac.getStdev().getValue(2)+" "+ac.getStdev().getValue(3));
        System.out.println(ac.getAverageLogs().getValue(0)+" "+ac.getAverageLogs().getValue(1)+" "+ac.getAverageLogs().getValue(2)+" "+ac.getAverageLogs().getValue(3));
        System.out.println(ac.getStdevLog().getValue(0)+" "+ac.getStdevLog().getValue(1)+" "+ac.getStdevLog().getValue(2)+" "+ac.getStdevLog().getValue(3));
    }
}
