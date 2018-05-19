/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataDouble;
import etomica.units.dimensions.Null;
import etomica.util.random.RandomMersenneTwister;

import static org.junit.jupiter.api.Assertions.*;

public abstract class AccumulatorAverageTestBase {

    public static void testSimple(AccumulatorAverage accumulator) {
        accumulator.putDataInfo(new DataDouble.DataInfoDouble("test", Null.DIMENSION));

        assertTrue(accumulator.getData().isNaN());
        assertEquals(accumulator.getBlockCount(), 0);

        simpleTest(accumulator);

        accumulator.reset();
        assertTrue(accumulator.getData().isNaN());
        assertEquals(accumulator.getBlockCount(), 0);

        simpleTest(accumulator);

        accumulator.putDataInfo(new DataDouble.DataInfoDouble("test", Null.DIMENSION));
        assertTrue(accumulator.getData().isNaN());
        assertEquals(accumulator.getBlockCount(), 0);

        simpleTest(accumulator);
    }

    // test accumulator on simple uncorrelated data
    public static void simpleTest(AccumulatorAverage accumulator) {
        DataDouble rawData = new DataDouble();
        RandomMersenneTwister rng = new RandomMersenneTwister(4);

        for (int i=0; i<1000000; i++) {
            rawData.x = rng.nextDouble();
            accumulator.putData(rawData);
        }

        IData accData = accumulator.getData();
        double avg = accData.getValue(accumulator.AVERAGE.index);
        assertTrue(Math.abs(avg-0.5) < 0.01, "average "+avg);

        double blockCorrelation = accData.getValue(accumulator.BLOCK_CORRELATION.index);
        // block correlation should be ~0, but actual value will depend on # of blocks
        assertTrue(Math.abs(blockCorrelation) < 4.0/Math.sqrt(accumulator.getBlockCount()), "block correlation "+blockCorrelation);

        double stdev = accData.getValue(accumulator.STANDARD_DEVIATION.index);
        assertTrue(Math.abs(stdev-Math.sqrt(1.0/12.0)) < 5.e-4, "standard devation "+stdev);

        double error = accData.getValue(accumulator.ERROR.index);
        assertTrue(error/2.9e-4 + 2.9e-4/error - 2 < 0.02, "error "+error);
    }

    public static void testSingleValue(AccumulatorAverage accumulator) {
        accumulator.putDataInfo(new DataDouble.DataInfoDouble("test", Null.DIMENSION));

        DataDouble rawData = new DataDouble();
        rawData.x = 5.6;
        for (int i=0; i<12345; i++) {
            accumulator.putData(rawData);
        }
        
        IData accData = accumulator.getData();
        double avg = accData.getValue(accumulator.AVERAGE.index);
        assertTrue(Math.abs(avg-rawData.x) < 1.e-10, "average "+avg);

        //double blockCorrelation = accData.getValue(AccumulatorAverage.StatType.BLOCK_CORRELATION.index);
        // block correlation should be 0/0, actual value might be 0, some number, Infinity or NaN 
        
        double stdev = accData.getValue(accumulator.STANDARD_DEVIATION.index);
        assertTrue(stdev < 1.e-7, "standard deviation ");

        double error = accData.getValue(accumulator.ERROR.index);
        assertTrue(error < 1.e-6, "error ");
    }
}
