/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataDouble;
import etomica.units.dimensions.Null;
import etomica.util.random.RandomNumberGenerator;

public class AccumulatorAverageCollapsingTest extends AccumulatorAverageTestBase {

    public AccumulatorAverageCollapsingTest() {
        super(new AccumulatorAverageCollapsing(400));
    }
    
    public void testSimple() {
        super.testSimple();
        
        accumulator.setBlockSize(500);

        assertTrue(accumulator.getData().isNaN());
        assertEquals(accumulator.getBlockCount(), 0);

        // the default
        accumulator.setBlockSize(1);

        assertTrue(accumulator.getData().isNaN());
        assertEquals(accumulator.getBlockCount(), 0);
        
        simpleTest();
    }
    
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
        assertEquals("block count "+accumulator.getBlockCount(), accumulator.getBlockCount(), 250);
        
        IData accData = accumulator.getData();
        double avg = accData.getValue(accumulator.AVERAGE.index);
        assertTrue("average "+avg, Math.abs(avg-0.5) < 0.015);
        
        double blockCorrelation = accData.getValue(accumulator.BLOCK_CORRELATION.index);
        assertTrue("block correlation "+blockCorrelation, Math.abs(blockCorrelation/0.27 + 0.27/blockCorrelation - 2) < 0.4);
        
        double stdev = accData.getValue(accumulator.STANDARD_DEVIATION.index);
        assertTrue("stdev "+stdev, Math.abs(stdev/0.046345 + 0.046345/stdev - 2) < 0.02);

        double error = accData.getValue(accumulator.ERROR.index);
        assertTrue("error "+error, error/0.0023 + 0.0023/error - 2 < 0.2);
    }
}
