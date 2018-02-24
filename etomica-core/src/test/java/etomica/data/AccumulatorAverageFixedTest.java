/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataDouble;
import etomica.units.dimensions.Null;
import etomica.util.random.RandomMersenneTwister;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static etomica.data.AccumulatorAverageTestBase.simpleTest;

class AccumulatorAverageFixedTest {
    private AccumulatorAverage accumulator;

    @BeforeEach
    public void setUp() {
        this.accumulator = new AccumulatorAverageFixed();
    }

    @Test
    public void testSimple() {
        AccumulatorAverageTestBase.testSimple(accumulator);
        
        // should cause a reset
        accumulator.setBlockSize(500);

        Assertions.assertTrue(accumulator.getData().isNaN());
        Assertions.assertEquals(accumulator.getBlockCount(), 0);

        // back to the default
        accumulator.setBlockSize(1000);

        simpleTest(accumulator);
    }

    @Test
    public void testCorrelated() {
        accumulator.putDataInfo(new DataDouble.DataInfoDouble("test", Null.DIMENSION));
        DataDouble rawData = new DataDouble();
        RandomMersenneTwister rng = new RandomMersenneTwister(4);
        accumulator.setBlockSize(10);

        for (int i=0; i<1000000; i++) {
            rawData.x = 0.5 + (rawData.x-0.5)*0.95 + (rng.nextDouble() - 0.5)*0.05;
            accumulator.putData(rawData);
        }
        
        IData accData = accumulator.getData();
        double avg = accData.getValue(accumulator.AVERAGE.index);
        Assertions.assertTrue(Math.abs(avg - 0.5) < 0.01, "average " + avg);
        
        double blockCorrelation = accData.getValue(accumulator.BLOCK_CORRELATION.index);
        Assertions.assertTrue(Math.abs(blockCorrelation - 0.7195) < 0.006, "block correlation " + blockCorrelation);
        
        double stdev = accData.getValue(accumulator.STANDARD_DEVIATION.index);
        Assertions.assertTrue(Math.abs(stdev - 0.046345) < 5.e-4, "stdev " + stdev);

        double error = accData.getValue(accumulator.ERROR.index);
        Assertions.assertTrue(error / 1.35e-4 + 1.35e-4 / error - 2 < 0.02, "error " + error);
    }
}
