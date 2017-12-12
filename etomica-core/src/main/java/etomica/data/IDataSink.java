/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;



/**
 * A recipient of Data.  Data goes in and might (or might not) come out.
 */
public interface IDataSink {

    /**
     * Gives data to DataSink for processing, display, or whatever it does.
     */
    void putData(IData data);

    /**
     * Informs the DataSink of the type of Data it should expect to receive.
     */
    void putDataInfo(IDataInfo dataInfo);
}
