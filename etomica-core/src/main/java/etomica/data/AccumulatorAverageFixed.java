/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.math.function.Function;
import etomica.math.function.IFunction;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

/**
 * AccumulatorAverage that maintains a fixed block size.
 * This is standard accumulator for collecting averages in production runs.
 * <p>
 * This accumulator accepts any type of Data
 * <p>
 * The average value is based on block averages -- it does not include samples
 * that are not part of blocks.  Average and error are computed using a
 * formulation that reduces numeric roundoff.
 * <p>
 * The standard deviation does include samples not included in a complete block
 * and uses a simpler formulation that is more vulnerable to roundoff.
 */
public class AccumulatorAverageFixed extends AccumulatorAverage {

    protected final IFunction negativeChop, sanityCheckBC;
    protected IData blockVarSum;//sum((blockAvg-totalAvg)^2)
    protected IData currentBlockAvg;
    protected IData sumSquare;//sum(value^2)
    protected IData mostRecentBlock, correlationSum, firstBlock;
    protected IData work, work2;
    protected IDataSink blockDataSink;
    protected String blockFilename;
    protected double[][] savedBlockData;
    protected int numSavedBlockData;

    /**
     * Default constructor sets block size to 1000 and sets the
     * interval for pushing the output data (pushInterval) to 100.
     */
    public AccumulatorAverageFixed() {
        this(1000);
    }

    /**
     * Constructs with interval for pushing the output data set to 100.
     *
     * @param blockSize size of the block.
     */
    public AccumulatorAverageFixed(long blockSize) {
        super(blockSize);
        // this chops negative values to be 0, used in cases where value
        // is non-zero due to roundoff and happens to be negative
        negativeChop = new IFunction() {
            public double f(double x) {
                return !(x <= 0) ? x : 0;
            }
        };

        // this provides a sanity check for the block correlation, which must
        // be between -1 and 1.  It can take strange (nonsense) values,
        // especially when the value is consistent (stdev=0).
        sanityCheckBC = new IFunction() {
            public double f(double x) {
                return (Double.isNaN(x) || x <= -1 || x >= 1) ? 0 : x;
            }
        };
    }

    /**
     * @return value of DoStrictBlockData
     */
    public boolean getDoStrictBlockData() {
        return true;
    }

    /**
     * @return the data sink where block data are passed.
     */
    public IDataSink getBlockDataSink() {
        return blockDataSink;
    }

    /**
     * Configures the accumulator to write block data to the given file.
     * The data can be read back in later so long as the underlying data
     * is DataDoubleArray.
     */
    public void setWriteBlocks(String filename) {
        blockFilename = filename;
        numSavedBlockData = 0;
        savedBlockData = new double[0][0];
    }

    /**
     * Reads block data from the given file.  The accumulator's state (for
     * block data) will match that when the file was last written to.
     */
    public void readBlockData(String filename) {
        if (!(currentBlockAvg instanceof DataDoubleArray)) {
            throw new RuntimeException("Can only read into DataDoubleArray");
        }
        try {
            FileReader fr = new FileReader(filename);
            BufferedReader bufReader = new BufferedReader(fr);
            String line;
            while ((line = bufReader.readLine()) != null) {
                String[] bits = line.split(" ");
                if (bits.length != mostRecentBlock.getLength()) {
                    throw new RuntimeException("Excepted " + mostRecentBlock.getLength() + " values, but found " + bits.length + " in file " + filename);
                }
                double[] x = ((DataDoubleArray) currentBlockAvg).getData();
                for (int i = 0; i < x.length; i++) {
                    x[i] = Double.parseDouble(bits[i]);
                }
                doBlockSum();
            }
            fr.close();
        } catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    public void writeBlockData() {
        if (blockFilename == null) throw new RuntimeException("no block file specified");
        try {
            FileWriter fw = new FileWriter(blockFilename, true);
            for (int i = 0; i < numSavedBlockData; i++) {
                fw.write("" + savedBlockData[i][0]);
                for (int j = 1; j < savedBlockData[i].length; j++) {
                    fw.write(" " + savedBlockData[i][j]);
                }
                fw.write("\n");
            }
            numSavedBlockData = 0;
            fw.close();
        } catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    /**
     * Allows individual block data to be passed to a data sink. Each block is
     * pushed to the sink after it is complete (as determined by block size).
     *
     * @param blockDataSink the sink where the blocks are to be passed. For
     *                      multiple streams, this could be a DataFork.
     */
    public void setBlockDataSink(IDataSink blockDataSink) {
        this.blockDataSink = blockDataSink;
        if (dataInfo != null) {
            blockDataSink.putDataInfo(((DataInfoGroup) dataInfo).getSubDataInfo(MOST_RECENT.index));
        }
    }

    /**
     * Add the given values to the sums and block sums. If any of the given data
     * values is NaN, method returns with no effect on accumulation sums.
     *
     * @param data data to be added.
     * @return true if data are not NaN.
     */
    public boolean addData(IData data) {
        if (data.isNaN())
            return false;

        mostRecent.E(data);
        work.E(data);
        work.ME(currentBlockAvg);
        work.TE(1.0 / (blockSize - blockCountDown + 1));
        currentBlockAvg.PE(work);
        work.E(data);
        work.TE(data);
        sumSquare.PE(work);
        if (--blockCountDown == 0) {//count down to zero to determine
            // completion of block
            doBlockSum();
            if (blockFilename != null) {
                if (numSavedBlockData == savedBlockData.length) {
                    savedBlockData = Arrays.copyOf(savedBlockData, numSavedBlockData + 1);
                    savedBlockData[numSavedBlockData] = new double[mostRecentBlock.getLength()];
                }
                for (int i = 0; i < mostRecentBlock.getLength(); i++) {
                    savedBlockData[numSavedBlockData][i] = mostRecentBlock.getValue(i);
                }
                numSavedBlockData++;
            }
            if (blockDataSink != null) {
                blockDataSink.putData(mostRecentBlock);
            }
        }
        return true;
    }
    
    public void putDataInfo(IDataInfo inputDataInfo) {
        super.putDataInfo(inputDataInfo);
        if (blockDataSink != null) {
            blockDataSink.putDataInfo(inputDataInfo);
        }
    }

    /**
     * Performs the block sum after <tt>blockSize</tt> calls to addData.
     */
    protected void doBlockSum() {
        count++;
        // we may have set average to NaN in getData()
        if (count == 1) average.E(0);
        work.E(currentBlockAvg);
        work.ME(average);
        work.TE(1.0 / count);
        work2.E(currentBlockAvg);
        work2.ME(average);
        average.PE(work);
        work.E(currentBlockAvg);
        work.E(currentBlockAvg);
        work.ME(average);
        work2.TE(work);
        blockVarSum.PE(work2);
        blockCountDown = blockSize;
        if (count > 1) {
            work.E(currentBlockAvg);
            work.TE(mostRecentBlock);
            correlationSum.PE(work);
        } else {
            firstBlock.E(currentBlockAvg);
        }
        //reset blocks
        mostRecentBlock.E(currentBlockAvg);
        currentBlockAvg.E(0.0);
    }

    public IData getData() {
        if (average == null)
            return null;
        if (count == 0) {
            average.E(Double.NaN);
        }
        if (count > 1) {
            // calculate other block properties (these require 2 or more blocks)

            error.E(blockVarSum);
            error.DE(count);
            error.map(negativeChop);

            // error's intermediate value is useful for calculating block correlation
            work.E(average);
            work.TE(average);
            blockCorrelation.E(average);
            blockCorrelation.TE(-2 * count);
            blockCorrelation.PE(firstBlock);
            blockCorrelation.PE(mostRecentBlock);
            blockCorrelation.TE(average);
            blockCorrelation.PE(correlationSum);
            blockCorrelation.DE(count - 1);
            blockCorrelation.PE(work);
            blockCorrelation.DE(error);
            blockCorrelation.map(sanityCheckBC);

            // ok, now finish up with error
            error.DE(count - 1);
            error.map(Function.Sqrt.INSTANCE);

            if (doIncludeACInError && count > 3) {
                work.E(1); // 1+c
                work.PE(blockCorrelation);
                work2.E(1);
                work2.ME(blockCorrelation); // 1-c
                work.DE(work2); // (1+c)/(1-c)
                work.map(Function.Sqrt.INSTANCE);
                error.TE(work);
            }
        } else {
            error.E(Double.NaN);
            blockCorrelation.E(Double.NaN);
        }

        long nTotalData = count * blockSize + (blockSize - blockCountDown);
        if (nTotalData > 0) {
            // now use *all* of the data
            work.E(currentBlockAvg);
            work.ME(average);
            work.TE((double) (blockSize - blockCountDown) / (blockSize * (1 + count) - blockCountDown));
            work.PE(average);
            work.TE(work);
            standardDeviation.E(sumSquare);
            standardDeviation.DE(nTotalData);
            standardDeviation.ME(work);
            standardDeviation.map(negativeChop);
            standardDeviation.map(Function.Sqrt.INSTANCE);
        } else {
            standardDeviation.E(Double.NaN);
        }
        return dataGroup;
    }

    public void reset() {
        super.reset();
        if (average == null) {
            return;
        }
        average.E(0);
        blockVarSum.E(0);
        currentBlockAvg.E(0);
        sumSquare.E(0);
        correlationSum.E(0);
        firstBlock.E(Double.NaN);
        mostRecentBlock.E(Double.NaN);
    }

    public IDataInfo processDataInfo(IDataInfo incomingDataInfo) {
        blockVarSum = incomingDataInfo.makeData();
        currentBlockAvg = incomingDataInfo.makeData();
        sumSquare = incomingDataInfo.makeData();
        firstBlock = incomingDataInfo.makeData();
        correlationSum = incomingDataInfo.makeData();
        mostRecentBlock = incomingDataInfo.makeData();
        work = incomingDataInfo.makeData();
        work2 = incomingDataInfo.makeData();
        return super.processDataInfo(incomingDataInfo);
    }
}
