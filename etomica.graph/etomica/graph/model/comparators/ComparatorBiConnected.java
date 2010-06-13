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