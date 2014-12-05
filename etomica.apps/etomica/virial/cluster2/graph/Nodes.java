/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph;

import java.util.List;
import java.util.Set;

public interface Nodes {

  public static final byte NODE_CLASS_FIELD = 0;
  public static final byte NODE_CLASS_ROOT = 1;
  public static final byte NODE_CLASS_DEFAULT = NODE_CLASS_FIELD;
  
  public static final char NODE_COLOR_0 = 'A';
  public static final char NODE_COLOR_1 = 'B';
  public static final char NODE_COLOR_2 = 'C';
  public static final char NODE_COLOR_3 = 'D';
  public static final char NODE_COLOR_4 = 'E';
  public static final char NODE_COLOR_5 = 'F';
  public static final char NODE_COLOR_6 = 'G';
  public static final char NODE_COLOR_7 = 'H';
  public static final char NODE_COLOR_8 = 'I';
  public static final char NODE_COLOR_9 = 'J';
  public static final char NODE_COLOR_DEFAULT = NODE_COLOR_0;
  
  // how many nodes are there?
  public byte count();

  // get the attributes of nodeID
  public NodeAttributes getAttributes(int nodeID);

  public Set<Character> getColors(byte nodeClass);

  public List<Integer> getPartition(byte nodeClass, char nodeColor);
}