/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

public interface IDataSource {

    /**
     * @return the data given by this source
     */
    IData getData();

    /**
     * Returns the DataInfo instance that will be held by Data
     * given by this source.  This information is useful for
     * setting up the data stream and for providing annotation when
     * displaying or writing the Data.
     */
//    public IDataInfo getDataInfo();

}
