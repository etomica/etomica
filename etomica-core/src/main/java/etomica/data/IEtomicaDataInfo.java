/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import com.fasterxml.jackson.annotation.JsonIgnore;
import etomica.units.dimensions.Dimension;

public interface IEtomicaDataInfo extends IDataInfo {

    /**
     * @return Returns the dimension given at construction.
     */
    public Dimension getDimension();
    
    /**
     * Adds the tag to this IDataInfo object.  You should only call this method
     * if you created this IDataInfo (You should not call this for IDataInfo
     * you receive.  Make a new instance with the DataInfoFactory and then add
     * your tag to that).
     */
    public void addTag(DataTag newTag);

    /**
     * Convenience method to add all of the given tags.  You should only call
     * this method if you created this IDataInfo (You should not call this for
     * IDataInfo you receive.  Make a new instance with the DataInfoFactory and
     * then add your tag to that).
     */
    public void addTags(DataTag[] newTags);

    /**
     * Removes all tags.  You should only call this method if you created this
     * IDataInfo.
     */
    public void clearTags();

    /**
     * Returns the array of tags held by this DataInfo instance.  The same
     * instance of the array can be returned each time the method is called,
     * so the array elements should not be modified.
     */
    @JsonIgnore
    public DataTag[] getTags();

    /**
     * Returns a mutable factory that can make copies of this instance of
     * DataInfo.
     */
    @JsonIgnore
    public IEtomicaDataInfoFactory getFactory();

}
