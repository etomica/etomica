/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import junit.framework.TestCase;
import etomica.data.types.DataDouble;
import etomica.units.dimensions.Null;
import etomica.util.random.RandomNumberGenerator;

public abstract class AccumulatorAverageTestBase extends TestCase {

    public AccumulatorAverageTestBase(AccumulatorAverage accumulator) {
        this.accumulator = accumulator;
    }

    public void testSimple() {
        accumulator.putDataInfo(new DataDouble.DataInfoDouble("test", Null.DIMENSION));

        assertTrue(accumulator.getData().isNaN());
        assertEquals(accumulator.getBlockCount(), 0);

        simpleTest();

        accumulator.reset();
        assertTrue(accumulator.getData().isNaN());
        assertEquals(accumulator.getBlockCount(), 0);

        simpleTest();

        accumulator.putDataInfo(new DataDouble.DataInfoDouble("test", Null.DIMENSION));
        assertTrue(accumulator.getData().isNaN());
        assertEquals(accumulator.getBlockCount(), 0);

        simpleTest();
    }

    // test accumulator on simple uncorrelated data
    public void simpleTest() {
        DataDouble rawData = new DataDouble();
        RandomNumberGenerator rng = new RandomNumberGenerator();

        for (int i=0; i<1000000; i++) {
            rawData.x = rng.nextDouble();
            accumulator.putData(rawData);
        }

        IData accData = accumulator.getData();
        double avg = accData.getValue(accumulator.AVERAGE.index);
        assertTrue("average "+avg, Math.abs(avg-0.5) < 0.01);

        double blockCorrelation = accData.getValue(accumulator.BLOCK_CORRELATION.index);
        // block correlation should be ~0, but actual value will depend on # of blocks
        assertTrue("block correlation "+blockCorrelation, Math.abs(blockCorrelation) < 4.0/Math.sqrt(accumulator.getBlockCount()));

        double stdev = accData.getValue(accumulator.STANDARD_DEVIATION.index);
        assertTrue("standard devation "+stdev, Math.abs(stdev-Math.sqrt(1.0/12.0)) < 5.e-4);

        double error = accData.getValue(accumulator.ERROR.index);
        assertTrue("error "+error, error/2.9e-4 + 2.9e-4/error - 2 < 0.02);
    }

    public void testSingleValue() {
        accumulator.putDataInfo(new DataDouble.DataInfoDouble("test", Null.DIMENSION));

        DataDouble rawData = new DataDouble();
        rawData.x = 5.6;
        for (int i=0; i<12345; i++) {
            accumulator.putData(rawData);
        }
        
        IData accData = accumulator.getData();
        double avg = accData.getValue(accumulator.AVERAGE.index);
        assertTrue("average "+avg, Math.abs(avg-rawData.x) < 1.e-10);

        //double blockCorrelation = accData.getValue(AccumulatorAverage.StatType.BLOCK_CORRELATION.index);
        // block correlation should be 0/0, actual value might be 0, some number, Infinity or NaN 
        
        double stdev = accData.getValue(accumulator.STANDARD_DEVIATION.index);
        assertTrue("standard deviation ", stdev < 1.e-7);

        double error = accData.getValue(accumulator.ERROR.index);
        assertTrue("error ", error < 1.e-6);
    }

    protected final AccumulatorAverage accumulator;
}
