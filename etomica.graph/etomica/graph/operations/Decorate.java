/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.Node;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;

/**
 * Performs decoration of diagrams by points.
 *
 * @author Andrew Schultz
 */
public class Decorate implements Binary {

  public Set<Graph> apply(Set<Graph> argument, Set<Graph> argument2, Parameters params) {

    assert (params instanceof DecorateParameters);
    char color = ((DecorateParameters)params).color;
    int factor = ((DecorateParameters)params).factor;
    MulFlexibleParameters mfp = ((DecorateParameters)params).mfp;
    Set<Graph> result = new HashSet<Graph>();
    List<Set<Graph>> allSet1 = new ArrayList<Set<Graph>>();
    for (Graph g : argument) {
      int colorCount = 0;
      if (factor == -1) {
        for (byte i=0; i<g.nodeCount(); i++) {
          if (g.getNode(i).getColor() == color) {
            colorCount++;
          }
        }
      }
      else {
        colorCount = g.factors()[factor];
      }
      if (colorCount == 0) {
        result.add(g);
        continue;
      }
      while (colorCount > allSet1.size()-1) {
        allSet1.add(new HashSet<Graph>());
      }
      allSet1.get(colorCount).add(g);
    }
    // set2Pow holds (set2)^i
    Set<Graph> set2Pow = new HashSet<Graph>();
    set2Pow.addAll(argument2);
    MulFlexible mulFlex = new MulFlexible();
    IsoFree isoFree = new IsoFree();
    for (int i=1; i<allSet1.size(); i++) {
      result.addAll(mulFlex.apply(allSet1.get(i), set2Pow, mfp));
      result = isoFree.apply(result, null);

      if (i+1<allSet1.size()) {
        // we're going to make another pass.  calculate (set2)^(i+1)
        set2Pow = isoFree.apply(mulFlex.apply(set2Pow, argument2, mfp), null);
      }
    }
    if (factor != -1) {
      for (Graph g : result) {
        g.factors()[factor] = 0;
      }
    }
    return result;
  }

  public static class DecorateParameters implements Parameters {
    public final char color;
    public final int factor;
    public final MulFlexibleParameters mfp;
    public DecorateParameters(char color, MulFlexibleParameters mfp) {
      this.color = color;
      this.mfp = mfp;
      factor = -1;
    }
    public DecorateParameters(int factor, MulFlexibleParameters mfp) {
      this.factor = factor;
      this.color = '-';
      this.mfp = mfp;
    }
  }
}
