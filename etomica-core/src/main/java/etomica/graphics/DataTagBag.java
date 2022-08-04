/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * 
 */
package etomica.graphics;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import etomica.data.DataTag;

/**
 * Simple class that associates a set of tags with an object (unit, label, etc) 
 * @author andrew
 */
public class DataTagBag {

    public DataTagBag(DataTag[] tags, Object obj) {
        this.tags = tags.clone();
        object = obj;
    }
    public final DataTag[] tags;
    public final Object object;

    /**
     * Returns the first DataTagBag from the given list that has no non-matching
     * tags with the given desiredTags
     */
    public static DataTagBag getDataTagBag(List<DataTagBag> dataTagBags, DataTag[] desiredTags) {
        // match the first set of tags that has no non-matching tags
        DataTagBag bestBag = null;
        int bestNumMatches = 0;
        for (DataTagBag tagUnit : dataTagBags) {
            int thisNumMatches = 0;
            DataTag[] tags = tagUnit.tags;
            boolean found = true;
            for (DataTag tag : tags) {
                found = false;
                for (DataTag desiredTag : desiredTags) {
                    if (desiredTag == tag) {
                        found = true;
                        thisNumMatches++;
                        break;
                    }
                }
                if (!found) {
                    break;
                }
            }
            if (found) {
                if (thisNumMatches > bestNumMatches) {
                    bestBag = tagUnit;
                    bestNumMatches = thisNumMatches;
                }
            }
        }
        
        return bestBag;
    }

    /**
     * Return the DataTagBag with the given set of tags (exactly) if it exists,
     * or null of no such DataTagBag exists.
     */
    public static DataTagBag getDataTagBagExact(LinkedList<DataTagBag> dataTagBags, DataTag[] exactTags) {
        Iterator<DataTagBag> iterator = dataTagBags.iterator();
        // match the first set of tags that has no non-matching tags
        while (iterator.hasNext()) {
            DataTagBag tagUnit = iterator.next();
            DataTag[] tags = tagUnit.tags;
            if (tags.length != exactTags.length) {
                continue;
            }
            boolean found = true;
            for (int j=0; j<tags.length; j++) {
                found = false;
                for (int k=0; k<exactTags.length; k++) {
                    if (exactTags[k] == tags[j]) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    break;
                }
            }
            if (found) {
                return tagUnit;
            }
        }
        
        return null;
    }
 
}