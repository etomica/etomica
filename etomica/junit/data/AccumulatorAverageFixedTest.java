package etomica.junit.data;

import etomica.api.IData;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.types.DataDouble;
import etomica.units.Null;
import etomica.util.RandomNumberGenerator;

public class AccumulatorAverageFixedTest extends AccumulatorAverageTestBase {

    public AccumulatorAverageFixedTest() {
        super(new AccumulatorAverageFixed());
    }

    public void testSimple() {
        super.testSimple();
        
        // should cause a reset
        accumulator.setBlockSize(500);

        assertTrue(accumulator.getData().isNaN());
        assertEquals(accumulator.getBlockCount(), 0);

        // back to the default
        accumulator.setBlockSize(1000);

        simpleTest();
    }

    public void testCorrelated() {
        accumulator.putDataInfo(new DataDouble.DataInfoDouble("test", Null.DIMENSION));
        DataDouble rawData = new DataDouble();
        RandomNumberGenerator rng = new RandomNumberGenerator();
        accumulator.setBlockSize(10);

        for (int i=0; i<1000000; i++) {
            rawData.x = 0.5 + (rawData.x-0.5)*0.95 + (rng.nextDouble() - 0.5)*0.05;
            accumulator.putData(rawData);
        }
        
        IData accData = accumulator.getData();
        double avg = accData.getValue(AccumulatorAverage.StatType.AVERAGE.index);
        assertTrue("average "+avg, Math.abs(avg-0.5) < 0.01);
        
        double blockCorrelation = accData.getValue(AccumulatorAverage.StatType.BLOCK_CORRELATION.index);
        // block correlation should be ~0, but actual value will depend on # of blocks 
        assertTrue("block correlation "+blockCorrelation, Math.abs(blockCorrelation-0.7196) < 0.005);
        
        double stdev = accData.getValue(AccumulatorAverage.StatType.STANDARD_DEVIATION.index);
        assertTrue("stdev "+stdev, Math.abs(stdev-0.046345) < 5.e-4);

        double error = accData.getValue(AccumulatorAverage.StatType.ERROR.index);
        assertTrue("error "+error, error/1.35e-4 + 1.35e-4/error - 2 < 0.02);
    }
}
