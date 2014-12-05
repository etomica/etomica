/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph;

import java.awt.Color;

public class ColorMap {

  public static Color COLOR_0 = Color.BLACK;
  public static Color COLOR_1 = Color.RED;
  public static Color COLOR_2 = Color.GREEN;
  public static Color COLOR_3 = Color.BLUE;
  public static Color COLOR_4 = Color.MAGENTA;
  public static Color COLOR_5 = Color.CYAN;
  public static Color COLOR_6 = Color.YELLOW;
  public static Color COLOR_7 = Color.ORANGE;
  public static Color COLOR_8 = Color.PINK;
  public static Color COLOR_9 = Color.LIGHT_GRAY;

  public static Color getNodeColor(char ch) {

    switch (ch) {
      case Nodes.NODE_COLOR_0: return COLOR_0;
      case Nodes.NODE_COLOR_1: return COLOR_1;
      case Nodes.NODE_COLOR_2: return COLOR_2;
      case Nodes.NODE_COLOR_3: return COLOR_3;
      case Nodes.NODE_COLOR_4: return COLOR_4;
      case Nodes.NODE_COLOR_5: return COLOR_5;
      case Nodes.NODE_COLOR_6: return COLOR_6;
      case Nodes.NODE_COLOR_7: return COLOR_7;
      case Nodes.NODE_COLOR_8: return COLOR_8;
      case Nodes.NODE_COLOR_9: return COLOR_9;
      default:
        return COLOR_0;
    }
  }

  public static Color getEdgeColor(char ch) {

    switch (ch) {
      case Edges.EDGE_COLOR_0: return COLOR_0;
      case Edges.EDGE_COLOR_1: return COLOR_1;
      case Edges.EDGE_COLOR_2: return COLOR_2;
      case Edges.EDGE_COLOR_3: return COLOR_3;
      case Edges.EDGE_COLOR_4: return COLOR_4;
      case Edges.EDGE_COLOR_5: return COLOR_5;
      case Edges.EDGE_COLOR_6: return COLOR_6;
      case Edges.EDGE_COLOR_7: return COLOR_7;
      case Edges.EDGE_COLOR_8: return COLOR_8;
      case Edges.EDGE_COLOR_9: return COLOR_9;
      default:
        return COLOR_0;
    }
  }
}
