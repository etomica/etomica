/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataDouble;
import etomica.units.dimensions.Null;
import etomica.util.random.RandomNumberGenerator;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static etomica.data.AccumulatorAverageTestBase.simpleTest;

class AccumulatorAverageCollapsingTest {
    private AccumulatorAverage accumulator;

    @BeforeEach
    public void setUp() {
        accumulator = new AccumulatorAverageCollapsing(400);
    }

    @Test
    public void testSimple() {
        AccumulatorAverageTestBase.testSimple(accumulator);
        
        accumulator.setBlockSize(500);

        Assertions.assertTrue(accumulator.getData().isNaN());
        Assertions.assertEquals(accumulator.getBlockCount(), 0);

        // the default
        accumulator.setBlockSize(1);

        Assertions.assertTrue(accumulator.getData().isNaN());
        Assertions.assertEquals(accumulator.getBlockCount(), 0);
        
        simpleTest(accumulator);
    }
    
    @Test
    public void testCorrelated() {
        accumulator.putDataInfo(new DataDouble.DataInfoDouble("test", Null.DIMENSION));
        DataDouble rawData = new DataDouble();
        RandomNumberGenerator rng = new RandomNumberGenerator();
        accumulator.setBlockSize(10);

        for (int i=0; i<10000; i++) {
            rawData.x = 0.5 + (rawData.x-0.5)*0.95 + (rng.nextDouble() - 0.5)*0.05;
            accumulator.putData(rawData);
        }
        
        // might change... we'd also need to update expectation values below
        Assertions.assertEquals(accumulator.getBlockCount(), 250, "block count " + accumulator.getBlockCount());
        
        IData accData = accumulator.getData();
        double avg = accData.getValue(accumulator.AVERAGE.index);
        Assertions.assertTrue(Math.abs(avg - 0.5) < 0.015, "average " + avg);
        
        double blockCorrelation = accData.getValue(accumulator.BLOCK_CORRELATION.index);
        Assertions.assertTrue(Math.abs(blockCorrelation / 0.27 + 0.27 / blockCorrelation - 2) < 0.4, "block correlation " + blockCorrelation);
        
        double stdev = accData.getValue(accumulator.STANDARD_DEVIATION.index);
        Assertions.assertTrue(Math.abs(stdev / 0.046345 + 0.046345 / stdev - 2) < 0.02, "stdev " + stdev);

        double error = accData.getValue(accumulator.ERROR.index);
        Assertions.assertTrue(error / 0.0023 + 0.0023 / error - 2 < 0.2, "error " + error);
    }
}
