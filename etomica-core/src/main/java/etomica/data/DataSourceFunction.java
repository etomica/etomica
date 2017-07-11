/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import java.io.Serializable;

import etomica.math.function.IFunction;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;
import etomica.math.function.Function;

/**
 * Datasource formed as a wrapper around a function.  Uses a DataSourceUniform
 * to generate x values, and takes a Function that computes
 * the y values from them.  This DataSource returns a DataGroup with
 * two Data components; the first (0) is x, and the second (1) is y.
 * Useful for displaying a fixed function on a plot.
 */
public class DataSourceFunction implements IEtomicaDataSource, DataSourceIndependent, Serializable {
    
    public DataSourceFunction() {
        this(new Function.Constant(0.0));
    }
    public DataSourceFunction(IFunction function) {
        this("y(x)", Null.DIMENSION, function, 100);
    }
    
    public DataSourceFunction(String label, Dimension dimension, IFunction function, int nValues) {
        this(label, dimension, function, nValues, "x", Null.DIMENSION);
    }
    
    public DataSourceFunction(String label, Dimension dimension, IFunction function, int nValues,
            String xLabel, Dimension xDimension) {
        xSource = new DataSourceUniform(xLabel, xDimension,nValues,0,1);
        this.function = function;
        setupData(label, dimension);
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    /**
     * Returns the DataSourceUniform instance that generates the x values.
     * The range and spacing of the x values can be adjusted through the
     * methods of the returned instance.
     */
    public DataSourceUniform getXSource() {
        return xSource;
    }

    /**
     * Returns the DataFunction made by this source.
     */
    public IData getData() { 
        return data;
    }
    
    /**
     * @return Returns the function.
     */
    public IFunction getFunction() {
        return function;
    }
    /**
     * @param function The function to set.
     */
    public void setFunction(IFunction function) {
        this.function = function;
        updateF();
    }
    /**
     * Recalculates the y values from the current x values.  This must be
     * invoked if the methods of the DataSourceUniform given by getXSource
     * are used to change the x values.
     *
     */
    protected void setupData(String label, Dimension dimension) {
        boolean needUpdate = false;
        if (xData != xSource.getData()) {
            xData = (DataDoubleArray)xSource.getData();
            needUpdate = true;
        }
        double[] x = xData.getData();
        if (data == null || data.getLength() != x.length) {
            needUpdate = true;
        }
        if (needUpdate) {
            data = new DataFunction(new int[]{xData.getArrayShape(0)});
            dataInfo = new DataInfoFunction(label, dimension, this);
        }
        updateF();
    }
    
    public DataDoubleArray getIndependentData(int i) {
        return xData;
    }
    
    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataInfoDoubleArray)xSource.getDataInfo();
    }
    
    public int getIndependentArrayDimension() {
        return 1;
    }
    
    public DataTag getIndependentTag() {
        return xSource.getTag();
    }
    
    /**
     * Updates the wrapped Data and DataInfo for change to the xDataSource
     */
    public void update() {
        setupData(dataInfo.getLabel(), dataInfo.getDimension());
    }
    
    public void updateF() {
        double[] x = xData.getData();
        double[] y = data.getData();
        for(int i=0; i<x.length; i++) {
            y[i] = function.f(x[i]);
        }
    }
    
    private static final long serialVersionUID = 1L;
    private DataFunction data;
    private IEtomicaDataInfo dataInfo;
    private final DataSourceUniform xSource;
    private DataDoubleArray xData;
    private IFunction function;
    protected final DataTag tag;
}
