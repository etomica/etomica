/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.data.*;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Length;
import etomica.util.random.IRandom;

/**
 * Meter that takes a (scalar) Meter and records its property as a
 * 1-dimensional function of position in the simulation volume. The measured
 * property must be a quantity that can be associated with a position in the
 * box. The position coordinate lies along one dimension (x,y,z).
 * 
 * @author Rob Riggleman
 * @author Andrew Schultz
 */
public class MeterProfile implements IDataSource, DataSourceIndependent, java.io.Serializable {
    
    /**
     * Default constructor sets profile along the x-axis, with 100 points in
     * the profile.
     */
    public MeterProfile(Space space, IRandom random) {
        xDataSource = new DataSourceUniform("x", Length.DIMENSION);
        tag = new DataTag();
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);
        positionSource = new RandomPositionSourceRectangular(space, random);
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }

    /**
     * The meter that defines the profiled quantity
     */
    public DataSourcePositioned getDataSource() {return meter;}
    
    /**
     * Accessor method for the meter that defines the profiled quantity.
     */
    public void setDataSource(DataSourcePositioned m) {
        if (!(m.getPositionDataInfo() instanceof DataInfoDouble)) {
            throw new IllegalArgumentException("data source must return a DataDouble");
        }
        meter = m;
        if (box != null) {
            meter.setBox(box);
        }
        reset();
    }
    
    /**
     * Accessor method for vector describing the direction along which the profile is measured.
     * Each atom position is dotted along this vector to obtain its profile abscissa value.
     */
    public int getProfileDim() {return profileDim;}
    
    /**
     * Accessor method for vector describing the direction along which the profile is measured.
     * Each atom position is dotted along this vector to obtain its profile abscissa value.
     * The given vector is converted to a unit vector, if not already.
     */
    public void setProfileDim(int dim) {
        profileDim = dim;
        reset();
    }
    
    /**
     * Returns the profile for the current configuration.
     */
    public IData getData() {
        data.E(0);
        double[] y = data.getData();
        IData xData = xDataSource.getData();
        
        for (int i=0; i<y.length; i++) {
            double x = xData.getValue(i);
            Vector pos = positionSource.randomPosition();
            pos.setX(profileDim, x);
            y[i] = meter.getData(pos).getValue(0);
        }
        return data;
    }

    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray)xDataSource.getData();
    }
    
    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataInfoDoubleArray)xDataSource.getDataInfo();
    }
    
    public int getIndependentArrayDimension() {
        return 1;
    }

    public DataTag getIndependentTag() {
        return xDataSource.getTag();
    }

    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
    public void setBox(Box box) {
        this.box = box;
        if (meter != null) {
            meter.setBox(box);
        }
        positionSource.setBox(box);
    }
    
    public void reset() {
        if (box == null) return;
        
        double halfBox = 0.5*box.getBoundary().getBoxSize().getX(profileDim);
        xDataSource.setXMin(-halfBox);
        xDataSource.setXMax(halfBox);
        
        if (meter != null) {
            data = new DataFunction(new int[] {xDataSource.getNValues()});
            dataInfo = new DataInfoFunction(meter.getPositionDataInfo().getLabel()+" Profile", meter.getPositionDataInfo().getDimension(), this);
            dataInfo.addTag(meter.getTag());
            dataInfo.addTag(tag);
        }
    }

    /**
     * Sets a new RandomPositionSource for this meter to use.  By default, a
     * position source is used which assumes rectangular boundaries.
     */
    public void setPositionSource(RandomPositionSource newPositionSource) {
        positionSource = newPositionSource;
        positionSource.setBox(box);
    }

    /**
     * Returns the RandomPositionSource used by this meter.
     */
    public RandomPositionSource getPositionSource() {
        return positionSource;
    }

    private static final long serialVersionUID = 1L;
    private Box box;
    private DataSourceUniform xDataSource;
    private DataFunction data;
    private IEtomicaDataInfo dataInfo;
    protected RandomPositionSource positionSource;
    /**
     * Vector describing the orientation of the profile.
     * For example, (1,0) is along the x-axis.
     */
    protected int profileDim;
    /**
     * Meter that defines the property being profiled.
     */
    protected DataSourcePositioned meter;
    protected final DataTag tag;
}
