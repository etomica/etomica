/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.model.comparators;

import static etomica.graph.model.Metadata.TYPE_NODE_ROOT;

import java.util.Comparator;

import etomica.graph.model.Graph;
import etomica.graph.model.Node;

public class ComparatorNumRootNodes implements Comparator<Graph> {

    public int compare(Graph g1, Graph g2) {
        int rootCount1 = 0;
        int rootCount2 = 0;
        for (Node node : g1.nodes()) {
            if (node.getType() == TYPE_NODE_ROOT) {
                rootCount1++;
            }
        }
        for (Node node : g2.nodes()) {
            if (node.getType() == TYPE_NODE_ROOT) {
                rootCount2++;
            }
        }
        return rootCount1 - rootCount2;
    }
}