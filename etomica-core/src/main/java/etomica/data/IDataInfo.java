/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

//import com.fasterxml.jackson.annotation.JsonIgnore;
import etomica.units.dimensions.Dimension;

public interface IDataInfo {

    /**
     * @return Returns the descriptive label of the data, as given in
     *         constructor or at last call to setLabel.
     */
    String getLabel();

    /**
     * Returns "label (dimension)", where label and dimension are the values
     * held by this instance.
     */
    String toString();

    /**
     * Returns the number of numerical values held by the IData object.
     */
    int getLength();

    /**
     * Returns a Data object appropriate for this DataInfo instance.
     */
    IData makeData();

    /**
     * @return Returns the dimension given at construction.
     */
    Dimension getDimension();

    /**
     * Adds the tag to this IDataInfo object.  You should only call this method
     * if you created this IDataInfo (You should not call this for IDataInfo
     * you receive.  Make a new instance with the DataInfoFactory and then add
     * your tag to that).
     */
    void addTag(DataTag newTag);

    /**
     * Convenience method to add all of the given tags.  You should only call
     * this method if you created this IDataInfo (You should not call this for
     * IDataInfo you receive.  Make a new instance with the DataInfoFactory and
     * then add your tag to that).
     */
    void addTags(DataTag[] newTags);

    /**
     * Removes all tags.  You should only call this method if you created this
     * IDataInfo.
     */
    void clearTags();

    /**
     * Returns the array of tags held by this DataInfo instance.  The same
     * instance of the array can be returned each time the method is called,
     * so the array elements should not be modified.
     */
//    @JsonIgnore
    DataTag[] getTags();

    /**
     * Returns a mutable factory that can make copies of this instance of
     * DataInfo.
     */
//    @JsonIgnore
    IDataInfoFactory getFactory();
}
