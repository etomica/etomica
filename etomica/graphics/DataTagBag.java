/**
 * 
 */
package etomica.graphics;

import java.util.Iterator;
import java.util.LinkedList;

import etomica.data.DataTag;

/**
 * Simple class that associates a set of tags with an object (unit, label, etc) 
 * @author andrew
 */
public class DataTagBag {

    public DataTagBag(DataTag[] tags, Object obj) {
        this.tags = (DataTag[])tags.clone();
        object = obj;
    }
    public final DataTag[] tags;
    public final Object object;

    /**
     * Returns the first DataTagBag from the given list that has no non-matching
     * tags with the given desiredTags
     */
    public static DataTagBag getDataTagBag(LinkedList dataTagBags, DataTag[] desiredTags) {
        Iterator iterator = dataTagBags.iterator();
        // match the first set of tags that has no non-matching tags
        DataTagBag bestBag = null;
        int bestNumMatches = 0;
        while (iterator.hasNext()) {
            int thisNumMatches = 0;
            DataTagBag tagUnit = (DataTagBag)iterator.next();
            DataTag[] tags = tagUnit.tags;
            boolean found = true;
            for (int j=0; j<tags.length; j++) {
                found = false;
                for (int k=0; k<desiredTags.length; k++) {
                    if (desiredTags[k] == tags[j]) {
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
    
}