/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.units.dimensions.Dimension;

/**
 * Particular data source for which the data is a simple scalar of type double.
 */

public abstract class DataSourceScalar implements IDataSource {
    
    public DataSourceScalar(String label, Dimension dimension) {
        data = new DataDouble();
        dataInfo = new DataInfoDouble(label, dimension);
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    /**
     * Returns a single scalar value as the measurement for the given box.
     * Subclasses define this method to specify the measurement they make.
     */
	public abstract double getDataAsScalar();
	
    
    /**
     * Causes the single getDataAsScalar(Box) value to be computed and
     * returned for the given box. In response to a getData() call,
     * MeterAbstract superclass will loop over all boxs previously specified
     * via setBox and collect these values into a vector and return them in
     * response to a getData() call.
     */
	public final DataDouble getDataDouble() {
		data.x = getDataAsScalar();
		return data;
	}
    
    public final IData getData() {
        return getDataDouble();
    }
	
	protected final DataDouble data;
    protected final IDataInfo dataInfo;
    protected final DataTag tag;
}
