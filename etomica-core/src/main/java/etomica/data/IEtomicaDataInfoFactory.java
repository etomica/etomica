/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.units.dimensions.Dimension;

public interface IEtomicaDataInfoFactory {

    /**
     * Creates a new DataInfo object using the information held by this factory.
     */
    public IEtomicaDataInfo makeDataInfo();

    /**
     * Sets the label
     */
    public void setLabel(String newLabel);

    /**
     * Returns the label
     */
    public String getLabel();

    /**
     * Sets the dimension
     */
    public void setDimension(Dimension newDimension);

    /**
     * Returns the dimension
     */
    public Dimension getDimension();

    /**
     * Returns the factory's tags as an array.  The array and its elements
     * should not be modified.  If modifications are needed, setTags should
     * be used.
     */
    public DataTag[] getTags();

    /**
     * Sets the factory's tags to those in the given array.
     */
    public void setTags(DataTag[] newTags);

}
