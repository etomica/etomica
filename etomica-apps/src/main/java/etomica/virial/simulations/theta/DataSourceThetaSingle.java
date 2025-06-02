/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial.simulations.theta;

import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.potential.compute.PotentialCompute;
import etomica.units.dimensions.Null;

public class DataSourceThetaSingle implements IDataSource {

    protected final Box box;
    protected final PotentialCompute pcLJ, pcdudk;
    protected final DataDoubleArray data;
    protected final DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected final double temperature, kBend;
    protected final DataTag tag;

    public DataSourceThetaSingle(Box box, double temperature, double kBend, PotentialCompute pcLJ, PotentialCompute pcdudk) {
        this.box = box;
        this.pcLJ = pcLJ;
        this.pcdudk = pcdudk;
        this.temperature = temperature;
        this.kBend = kBend;
        data = new DataDoubleArray(4);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("Stuff", Null.DIMENSION, new int[]{4});
        tag = new DataTag();
    }

    @Override
    public IData getData() {
        double[] y = data.getData();
        double ulj = pcLJ.computeAll(false);
        double dudk = pcdudk == null ? 0 : pcdudk.computeAll(false);
        double u1 = ulj + kBend*dudk;
        y[0] = -u1;
        y[1] = -dudk/temperature;
        y[2] = (u1/temperature - 1) * dudk;
        y[3] = u1*u1;
        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
}
