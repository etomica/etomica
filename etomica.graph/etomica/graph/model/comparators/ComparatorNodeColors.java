/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.model.comparators;

import java.util.Comparator;

import etomica.graph.model.Graph;
import etomica.graph.model.Node;

/**
 * Compares based on graph factors (when available) or node colors.  If factors
 * exist, they are used.
 * (2,0) < (1,0) < (0,1)
 * 
 * If factors do not exist, then the number of nodes of each color are counted.
 * A graph with more nodes of the first color will come first.
 * 
 * @author Andrew Schultz
 */
public class ComparatorNodeColors implements Comparator<Graph> {

    protected final char[] colors;
    protected final int[][] colorCounts;
    
    public ComparatorNodeColors() {
        colors = null;
        colorCounts = null;
    }
    
    public ComparatorNodeColors(char[] colors) {
        this.colors = colors;
        colorCounts = new int[2][colors.length];
    }
    
    public int compare(Graph g1, Graph g2) {
        int[] factors1 = g1.factors();
        int[] factors2 = g2.factors();
        if (factors1.length > 0) {
            for (int i=0; i<factors1.length; i++) {
                int c = factors2[i] - factors1[i];
                if (c != 0) return c;
            }
            return 0;
        }
        if (colors == null) {
            return 0;
        }
        for (int i=0; i<colors.length; i++) {
            colorCounts[0][i] = 0;
            colorCounts[1][i] = 0;
        }
        for (Node node : g1.nodes()) {
            char nodeColor = node.getColor();
            for (int i=0; i<colors.length; i++) {
                if (nodeColor == colors[i]) {
                    colorCounts[0][i]++;
                }
            }
        }
        for (Node node : g2.nodes()) {
            char nodeColor = node.getColor();
            for (int i=0; i<colors.length; i++) {
                if (nodeColor == colors[i]) {
                    colorCounts[1][i]++;
                }
            }
        }
        for (int i=0; i<colors.length; i++) {
            int c = colorCounts[1][i] - colorCounts[0][i];
            if (c != 0) return c;
        }
        return 0;
    }
}