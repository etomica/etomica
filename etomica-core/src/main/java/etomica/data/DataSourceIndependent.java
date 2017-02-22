/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;

/**
 * Interface for an object that can provide X Data objects on request.  This
 * provides the independent data part of a DataFunction, which downstream data
 * elements can retrieve when they need it.
 */ 
public interface DataSourceIndependent {

    /**
     * Returns the X data for the given dimension
     */
	public DataDoubleArray getIndependentData(int i);

    /**
     * Returns the DataInfo for the given dimension
     */
    public DataInfoDoubleArray getIndependentDataInfo(int i);
    
    /**
     * Returns the number of independent data dimensions
     */
    public int getIndependentArrayDimension();
    
    /**
     * Returns the tag associated with this DataSource.
     */
    public DataTag getIndependentTag();
}