/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;



/**
 * Interface for a class that performs some action in response to
 * a change in a DataSet instance.  Classes implementing this interface can
 * register with a DataSet via the addDataListener method.  Typically
 * this interface is used to cause a display element to update when the
 * data is updated.
 * 
 * @author Andrew Schultz
 */
public interface DataSetListener {

    /**
     * Method called when one or more pieces of data have changed.
     */
    public void dataChanged(DataSet dataSet);
    
    /**
     * Method called when a Data object is added.  Firing is
     * performed after the object is added to the DataSet.
     * @param newData the Data object that has been added
     */
    public void dataCountChanged(DataSet dataSet);
}
