/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

/**
 * Interface for an object that can provide IData objects
 * on request.  Normally an IDataSource heads a stream that processes
 * and/or records the IData as it passes through different segments.
 * IData is pulled from the IDataSource by a DataPump and pushed down the
 * data stream.
 */
public interface IDataSource {

    /**
     * @return the data given by this source
     */
    IData getData();

    /**
     * @return a tag that uniquely identifies the IDataInfo
     */
    DataTag getTag();

    /**
     * @returns the IEtomicaDataInfo instance associated with this source.
     * This information is useful for setting up the data stream and for
     * providing annotation when displaying or writing the IData.
     */
    IEtomicaDataInfo getDataInfo();
}
