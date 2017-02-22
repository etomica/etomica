/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.iterators.filters;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;
import etomica.graph.model.Metadata;
import etomica.graph.model.Node;

/**
 * Filters by number of field nodes.  Useful for truncating a series, perhaps
 * after performing multiplication.
 * 
 * @author Andrew Schultz
 */
public class FieldNodeCount extends LocalFilter {

    public FieldNodeCount(GraphIterator iterator, int maxCount) {
        super(iterator);
        this.maxCount = maxCount;
    }

    protected boolean accept(Graph g) {
        int fieldCount = 0;
        for (Node node : g.nodes()) {
            if (node.getType() == Metadata.TYPE_NODE_FIELD) {
                fieldCount++;
            }
        }
        return fieldCount <= maxCount;
    }

    protected final int maxCount;
}
