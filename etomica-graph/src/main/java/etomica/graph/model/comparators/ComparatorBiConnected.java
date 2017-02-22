/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.model.comparators;

import java.util.Comparator;

import etomica.graph.model.Graph;
import etomica.graph.property.IsBiconnected;

public class ComparatorBiConnected implements Comparator<Graph> {

    public int compare(Graph g1, Graph g2) {
        boolean isBiCon1 = isBi.check(g1);
        boolean isBiCon2 = isBi.check(g2);
        if (isBiCon1 == isBiCon2) return 0;
        return isBiCon1 ? -1 : 1;
    }
    
    protected final IsBiconnected isBi = new IsBiconnected();
}