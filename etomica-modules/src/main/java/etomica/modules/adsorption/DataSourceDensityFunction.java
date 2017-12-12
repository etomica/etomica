/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.adsorption;

import etomica.data.*;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Quantity;

public class DataSourceDensityFunction implements IDataSource {

    protected final EOSSW eos;
    protected final DataFunction data;
    protected final DataInfoFunction dataInfo;
    protected final DataSourceIndependentSimple xDataSource;
    protected final DataTag tag;
    protected final boolean isLiquid;
    protected double pressure;
    
    public DataSourceDensityFunction(EOSSW eos, double xMin, double xMax) {
        this(eos, xMin, xMax, false);
    }
    
    public DataSourceDensityFunction(EOSSW eos, double xMin, double xMax, boolean isLiquid) {
        this.eos = eos;
        this.isLiquid = isLiquid;
        data = new DataFunction(new int[]{2});
        double[] xData = new double[]{xMin,xMax};
        DataInfoDoubleArray xDataInfo = new DataInfoDoubleArray("x", Length.DIMENSION, new int[]{2});
        xDataSource = new DataSourceIndependentSimple(xData, xDataInfo);
        dataInfo = new DataInfoFunction("homogenous density", Quantity.DIMENSION, xDataSource);
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public void setPressure(double p) {
        pressure = p;
    }
    
    public IData getData() {
        double[] yData = data.getData();
        double p = pressure;
        if (isLiquid) {
            p = eos.pSat()*1.000001;
            if (p < pressure) p = pressure;
        }
        double rho = eos.rhoForPressure(p);
        yData[0] = rho;
        yData[1] = rho;
        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

}
