/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.iterators.DefaultIterator;
import etomica.graph.iterators.filters.PropertyFilter;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;
import etomica.graph.model.GraphList;
import etomica.graph.model.comparators.ComparatorBiConnected;
import etomica.graph.model.comparators.ComparatorChain;
import etomica.graph.model.comparators.ComparatorNumEdges;
import etomica.graph.model.comparators.ComparatorNumFieldNodes;
import etomica.graph.model.comparators.ComparatorNumNodes;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.IsoFree;
import etomica.graph.property.IsBiconnected;
import etomica.graph.viewer.ClusterViewer;

/**
 * Main method that simply generates and displays the graphs in the direct
 * correlation function.
 * 
 * @author andrew
 */
public class DCF {

    public static void main(String[] args) {
        IsBiconnected isBi = new IsBiconnected();
        ComparatorChain comp = new ComparatorChain();
        comp.addComparator(new ComparatorNumFieldNodes());
        comp.addComparator(new ComparatorBiConnected());
        comp.addComparator(new ComparatorNumEdges());
        comp.addComparator(new ComparatorNumNodes());
        GraphList graphList = new GraphList(comp);
        for (int i=2; i<6; i++) {
            GraphIterator iter = new PropertyFilter(new DefaultIterator((byte)i, (byte)2), isBi);
            while (iter.hasNext()) {
                Graph g = iter.next();
                graphList.add(g);
            }
        }
        MetadataImpl.rootPointsSpecial = true;
        Set<Graph> condensed = new HashSet<Graph>();
        IsoFree isoFree = new IsoFree();
        condensed.addAll(isoFree.apply(graphList, null));
        graphList.clear();
        graphList.addAll(condensed);
        for (Graph g : graphList) {
            System.out.println(g);
        }
        ClusterViewer.createView("DCF", graphList);
    }
}
