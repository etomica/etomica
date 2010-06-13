package etomica.graph.model.comparators;

import java.util.Comparator;

import etomica.graph.model.Graph;

public class ComparatorNumEdges implements Comparator<Graph> {

    public int compare(Graph g1, Graph g2) {
        return g2.edgeCount() - g1.edgeCount();
    }
}