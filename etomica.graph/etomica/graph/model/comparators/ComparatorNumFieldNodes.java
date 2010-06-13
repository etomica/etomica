package etomica.graph.model.comparators;

import static etomica.graph.model.Metadata.TYPE_NODE_FIELD;

import java.util.Comparator;

import etomica.graph.model.Graph;
import etomica.graph.model.Node;

public class ComparatorNumFieldNodes implements Comparator<Graph> {

    public int compare(Graph g1, Graph g2) {
        int fieldCount1 = 0;
        int fieldCount2 = 0;
        for (Node node : g1.nodes()) {
            if (node.getType() == TYPE_NODE_FIELD) {
                fieldCount1++;
            }
        }
        for (Node node : g2.nodes()) {
            if (node.getType() == TYPE_NODE_FIELD) {
                fieldCount2++;
            }
        }
        return fieldCount1 - fieldCount2;
    }
}