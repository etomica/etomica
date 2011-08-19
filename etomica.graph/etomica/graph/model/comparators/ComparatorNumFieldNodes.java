package etomica.graph.model.comparators;

import java.util.Comparator;

import etomica.graph.model.Graph;
import etomica.graph.property.NumFieldNodes;

public class ComparatorNumFieldNodes implements Comparator<Graph> {

    public int compare(Graph g1, Graph g2) {
        int fieldCount1 = NumFieldNodes.value(g1);
        int fieldCount2 = NumFieldNodes.value(g2);
        return fieldCount1 - fieldCount2;
    }
}