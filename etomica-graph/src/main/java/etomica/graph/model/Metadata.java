/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public interface Metadata extends Comparable<Metadata> {

  public static final char TYPE_EDGE_ANY = 'E';
  public static final char TYPE_NODE_ROOT = 'R';
  public static final char TYPE_NODE_FIELD = 'F';
  public static final char TYPE_NODE_DEFAULT = TYPE_NODE_FIELD;

  // COLOR_CODES is initially populated with A-J.  Further codes can be added
  // dynamically
  public static final char COLOR_CODE_0 = 'A';
  public static final char COLOR_CODE_1 = 'B';
  public static final char COLOR_CODE_2 = 'C';
  public static final char COLOR_CODE_3 = 'D';
  public static final char COLOR_CODE_4 = 'E';
  public static final char COLOR_CODE_5 = 'F';
  public static final char COLOR_CODE_6 = 'G';
  public static final char COLOR_CODE_7 = 'H';
  public static final char COLOR_CODE_8 = 'I';
  public static final char COLOR_CODE_9 = 'J';
  public static final char COLOR_CODE_DEFAULT = COLOR_CODE_0;
  public static List<Character> COLOR_CODES = new ArrayList<Character>(Arrays.asList(new Character[] { COLOR_CODE_0, COLOR_CODE_1,
      COLOR_CODE_2, COLOR_CODE_3, COLOR_CODE_4, COLOR_CODE_5, COLOR_CODE_6, COLOR_CODE_7, COLOR_CODE_8,
      COLOR_CODE_9}));

  // http://www.december.com/html/spec/colorsvg.html

  // COLORS is populated with 13 colors initially, which should accomodate a
  // few custom codes.  More colors can be added dynamically, but should be
  // valid SVG colors as described at the above URL.
  public static String COLOR_0 = "black";
  public static String COLOR_1 = "red";
  public static String COLOR_2 = "green";
  public static String COLOR_3 = "darkviolet";
  public static String COLOR_4 = "royalblue";
  public static String COLOR_5 = "darkorange";
  public static String COLOR_6 = "fuchsia";
  public static String COLOR_7 = "yellow";
  public static String COLOR_8 = "darkblue";
  public static String COLOR_9 = "slategray";
  public static String COLOR_10 = "orange";
  public static String COLOR_11 = "maroon";
  public static String COLOR_12 = "olive";
  public final static ArrayList<String> COLORS = new ArrayList<String>(Arrays.asList(new String[] { COLOR_0, COLOR_1, COLOR_2, COLOR_3, COLOR_4, COLOR_5,
      COLOR_6, COLOR_7, COLOR_8, COLOR_9, COLOR_10, COLOR_11, COLOR_12 }));

  public final static HashMap<Character,String> COLOR_MAP = new HashMap<Character,String>();

  public final static HashMap<Character,Integer> DASH_MAP = new HashMap<Character,Integer>();

  public final static HashMap<Character,Integer> BOND_ORDER_MAP = new HashMap<Character,Integer>();

  public char getColor();

  public char getType();

  public boolean isCompatible(Metadata other);

  public boolean isSameColor(Metadata other);

  public boolean isSameType(Metadata other);
}
