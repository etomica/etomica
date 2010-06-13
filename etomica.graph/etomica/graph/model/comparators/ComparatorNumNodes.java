package etomica.graph.model.comparators;

import java.util.Comparator;

import etomica.graph.model.Graph;

public class ComparatorNumNodes implements Comparator<Graph> {

    public int compare(Graph g1, Graph g2) {
        return g1.nodeCount() - g2.nodeCount();
    }
}