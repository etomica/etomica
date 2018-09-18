/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.types;

import etomica.data.IData;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

class DataGroupTest {

    protected final DataGroup data;

    public DataGroupTest() {
        DataDoubleArray data0 = new DataDoubleArray(new int[]{2}, new double[]{0,1});
        DataDoubleArray data1 = new DataDoubleArray(new int[]{3}, new double[]{2,3,4});
        DataDoubleArray data2 = new DataDoubleArray(new int[]{4}, new double[]{5,6,7,8});
        data = new DataGroup(new IData[]{data0, data1, data2});
    }

    @Test
    public void testDataGroup() {
        testGetValue(data, new double[]{0,1,2,3,4,5,6,7,8});
    }

    public void testGetValue(IData someData, double[] expectedValues) {
        Assertions.assertEquals(someData.getLength(), expectedValues.length);
        for (int i=0; i<expectedValues.length; i++) {
            Assertions.assertEquals(expectedValues[i], someData.getValue(i));
        }
    }
}
