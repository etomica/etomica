/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;


/**
 * Interface for an object that can provide Data objects
 * on request.  Normally a DataSource heads a stream that processes
 * and/or records the Data as it passes through different segments.  
 * Data is pulled from the DataSource by a DataPump and pushed down the
 * data stream.
 */
 
public interface IEtomicaDataSource extends IDataSource {

    public DataTag getTag();

    public IEtomicaDataInfo getDataInfo();
}