/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.Metadata;
import etomica.graph.model.Node;
import etomica.graph.operations.Parameters;
import etomica.graph.operations.Unary;
import etomica.graph.traversal.CVisitor;

/**
 * This operation replaces groups of points connected (internally) by
 * e-exchange bonds with a single point.  Bonds with other points (or groups of
 * points) are assumed to be f-exchange bonds (meaning that the group of points
 * has a single f-bond with the other point (or groups of points) and these
 * f-exchange bonds are replaced by a single f-bond.
 *
 * Exchange group points will have a color that corresponds to the number of
 * exchange points, color1 + numPoints.  color1 is assumed to be upper case;
 * Where possible, root points will be moved out of exchange groups.  Where not
 * possible, the resulting color will be a lower-case version of
 * color1 + numPoints.  color1 is assumed to be upper-case.
 * 
 * @author Andrew Schultz
 */
public class CondenseExchange implements Unary {

    public Set<Graph> apply(Set<Graph> argument, Parameters params) {
        Set<Graph> result = new HashSet<Graph>();
        for (Graph g : argument) {
            result.add(apply(g, (CondenseExchangeParameters)params));
        }
        return result;
    }

    public Graph apply(Graph g, CondenseExchangeParameters params) {
        Graph gCopy = g.copy();
        byte n = g.nodeCount();
        for (byte i=0; i<n*(n-1)/2; i++) {
            if (gCopy.hasEdge(i) && gCopy.getEdge(i).getColor() != params.eExchange) gCopy.deleteEdge(i);
        }
        List<List<Byte>> comps = CVisitor.getComponents(gCopy);
        // comps is our list of "molecules"
        // visit each component and find the bonds coming from that component.
        // we only need to look at one node from each component (they must
        // all have the same bonds).
        
        if (comps.size() == n) return g.copy();

        byte extraRoots = 0;
        for (byte i=0; i<comps.size(); i++) {
            List<Byte> iComp = comps.get(i);
            if (iComp.size() == 1 || iComp.size() == n) continue;
            boolean rooty = false;
            for (byte j=0; j<iComp.size(); j++) {
                if (g.getNode(iComp.get(j)).getType() == Metadata.TYPE_NODE_ROOT) {
                    rooty = true;
                    break;
                }
            }
            if (rooty) {
                extraRoots++;
            }
        }
        
        char color1 = params.color1;
        char color1Root = (char)(params.color1 + ('a' - 'A'));
        
        Graph gNew = GraphFactory.createGraph((byte)comps.size());
        for (byte i=0; i<comps.size(); i++) {
            List<Byte> iComp = comps.get(i);
            if (iComp.size() == 1 && g.getNode(iComp.get(0)).getType() == Metadata.TYPE_NODE_ROOT) {
                // it was rooty in the original graph
                gNew.getNode(i).setType(Metadata.TYPE_NODE_ROOT);
            }
            char nodeColor = g.getNode(iComp.get(0)).getColor();
            if (iComp.size() > 1) {
                // we want to use the original color unless this is an exchanged component
                nodeColor = (char)(color1 + (iComp.size()-1));
                for (byte j=0; j<iComp.size(); j++) {
                    if (g.getNode(iComp.get(j)).getType() == Metadata.TYPE_NODE_ROOT) {
                        nodeColor = (char)(color1Root + (iComp.size()-1));
                        break;
                    }
                }
            }
            gNew.getNode(i).setColor(nodeColor);

            byte iNode = iComp.get(0);
            for (byte j=(byte)(i+1); j<comps.size(); j++) {
                List<Byte> jComp = comps.get(j);
                byte jNode = jComp.get(0);
                if (!g.hasEdge(iNode, jNode)) continue;
                gNew.putEdge(i,j);
                char oldColor = g.getEdge(iNode,jNode).getColor();
                if (oldColor == params.fExchange) {
                    gNew.getEdge(i,j).setColor(params.f);
                }
                else {
                    // these are two single points
                    if (iComp.size() != 1 || jComp.size() != 1) {
                        throw new RuntimeException("oops");
                    }
                    gNew.getEdge(i,j).setColor(oldColor);
                }
            }                
        }
        byte rootsAdded = 0;
        for (byte i=0; i<comps.size(); i++) {
            if (extraRoots == rootsAdded) break;
            Node iNode = gNew.getNode(i);
            if (iNode.getType() == Metadata.TYPE_NODE_ROOT || iNode.getColor() != color1) continue;
            iNode.setType(Metadata.TYPE_NODE_ROOT);
            rootsAdded++;
        }
        for (byte i=0; i<comps.size(); i++) {
            if (rootsAdded == 0) break;
            Node iNode = gNew.getNode(i);
            char iColor = iNode.getColor();
            if (iColor >= 'a' && iColor <= 'z') {
                // found a root point.  convert back to upper-case (not root)
                iNode.setColor((char)(iColor + ('A' - 'a')));
                rootsAdded--;
            }
        }
        gNew.coefficient().multiply(g.coefficient());
        return gNew;
    }

    public static class CondenseExchangeParameters implements Parameters {
        public final char fExchange, eExchange, f, color1;
        public CondenseExchangeParameters(char fExchange, char eExchange, char f, char color1) {
            this.fExchange = fExchange;
            this.eExchange = eExchange;
            this.f = f;
            this.color1 = color1;
        }
    }
}
