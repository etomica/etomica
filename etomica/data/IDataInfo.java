/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

public interface IDataInfo {

    /**
     * @return Returns the descriptive label of the data, as given in
     *         constructor or at last call to setLabel.
     */
    public String getLabel();

    /**
     * Returns "label (dimension)", where label and dimension are the values
     * held by this instance.
     */
    public String toString();

    /**
     * Returns the number of numerical values held by the IData object.
     */
    public int getLength();

    /**
     * Returns a Data object appropriate for this DataInfo instance.
     */
    public IData makeData();

}