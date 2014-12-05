/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.model.comparators;

import java.util.Comparator;

import etomica.graph.model.Graph;

public class ComparatorNumEdges implements Comparator<Graph> {

    public int compare(Graph g1, Graph g2) {
        return g2.edgeCount() - g1.edgeCount();
    }
}