/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import java.awt.Color;
import java.util.*;

import etomica.virial.cluster2.mvc.WizardState;

public class ClusterWizardState implements WizardState {

  // page 1 data: global keys
  private static final String KEY_CLUSTER = "cluster";
  public static final String KEY_NAME = KEY_CLUSTER + ".name";
  public static final String KEY_TOTAL_NODES = KEY_CLUSTER + ".totalNodes";
  public static final String KEY_ROOT_NODES = KEY_CLUSTER + ".rootNodes";
  public static final String KEY_COLOR_SCHEME = KEY_CLUSTER + ".colorScheme";
  public static final String KEY_ISOMORPH_FREE = KEY_CLUSTER + ".isomorphFree";
  // page 1 data: global defaults
  public static final String DEFVAL_NAME = "MyCluster1"; // "NewCluster";
  public static final Integer DEFVAL_TOTAL_NODES = 6; // 4;
  public static final Integer DEFVAL_ROOT_NODES = 3; // 2;
  public static final String DEFVAL_MONOCHROMATIC = "monochromatic";
  public static final String DEFVAL_MULTICOLORED = "multi-colored";
  public static final String DEFVAL_COLOR_SCHEME = DEFVAL_MULTICOLORED; // DEFVAL_MONOCHROMATIC;
  public static final Boolean DEFVAL_ISOMORPH_FREE = false; // true

  // page 2 data: connectivity keys
  private static final String KEY_CONNECTIVITY = "connectivity";
  public static final String KEY_CLASS_ANY = KEY_CONNECTIVITY + ".any";
  public static final String KEY_CLASS_CONNECTED = KEY_CONNECTIVITY + ".connected";
  public static final String KEY_CLASS_BICONNECTED = KEY_CONNECTIVITY + ".biconnected";
  public static final String KEY_CLASS_REEHOOVER = KEY_CONNECTIVITY + ".reehoover";
  public static final String KEY_EXCLUDE_NODAL_POINTS = KEY_CONNECTIVITY + ".nodalPoints";
  public static final String KEY_EXCLUDE_ARTICULATION_POINTS = KEY_CONNECTIVITY + ".articulationPoints";
  public static final String KEY_EXCLUDE_ARTICULATION_PAIRS = KEY_CONNECTIVITY + ".articulationPairs";
  // page 2 data: connectivity defaults
  public static final Boolean DEFVAL_CLASS_ANY = true;
  public static final Boolean DEFVAL_CLASS_CONNECTED = false;
  public static final Boolean DEFVAL_CLASS_BICONNECTED = false;
  public static final Boolean DEFVAL_CLASS_REEHOOVER = false;
  public static final Boolean DEFVAL_EXCLUDE_NODAL_POINTS = false;
  public static final Boolean DEFVAL_EXCLUDE_ARTICULATION_POINTS = false;
  public static final Boolean DEFVAL_EXCLUDE_ARTICULARION_PAIRS = false;

  // page 3 data: color map keys
  private static final String KEY_COLOR = "color";
  public static final String KEY_COLOR1 = KEY_COLOR + ".1";
  public static final String KEY_COLOR2 = KEY_COLOR + ".2";
  public static final String KEY_COLOR3 = KEY_COLOR + ".3";
  public static final String KEY_COLOR4 = KEY_COLOR + ".4";
  public static final String KEY_COLOR5 = KEY_COLOR + ".5";
  public static final String KEY_COLOR6 = KEY_COLOR + ".6";
  public static final String KEY_COLOR7 = KEY_COLOR + ".7";
  public static final String KEY_COLOR8 = KEY_COLOR + ".8";
  public static final String KEY_COLOR9 = KEY_COLOR + ".9";
  public static final String KEY_COLORA = KEY_COLOR + ".A";
  public static final String KEY_COLORB = KEY_COLOR + ".B";
  public static final String KEY_COLORC = KEY_COLOR + ".C";
  public static final String KEY_COLORD = KEY_COLOR + ".D";
  public static final String KEY_COLORE = KEY_COLOR + ".E";
  public static final String KEY_COLORF = KEY_COLOR + ".F";
  public static final List<String> KEY_COLORS = Arrays.asList(new String[] { KEY_COLOR1, KEY_COLOR2,
      KEY_COLOR3, KEY_COLOR4, KEY_COLOR5, KEY_COLOR6, KEY_COLOR7, KEY_COLOR8, KEY_COLOR9, KEY_COLORA,
      KEY_COLORB, KEY_COLORC, KEY_COLORD, KEY_COLORE, KEY_COLORF });
  public static final String KEY_MAPPED_COLORS = KEY_COLOR + ".mapping";
  // page 3 data: color map defaults (by name)
  public static final ColorEntry DEFVAL_BLACK = new ColorEntry(Color.black, "black");
  public static final ColorEntry DEFVAL_DARKGRAY = new ColorEntry(Color.darkGray, "dark gray");
  public static final ColorEntry DEFVAL_SLATEGRAY = new ColorEntry(new Color(112, 128, 144), "slate gray");
  public static final ColorEntry DEFVAL_LIGHTGRAY = new ColorEntry(Color.lightGray, "light gray");
  public static final ColorEntry DEFVAL_NAVY = new ColorEntry(new Color(0, 0, 128), "navy");
  public static final ColorEntry DEFVAL_BLUE = new ColorEntry(Color.blue, "blue");
  public static final ColorEntry DEFVAL_CYAN = new ColorEntry(Color.cyan, "cyan");
  public static final ColorEntry DEFVAL_MAROON = new ColorEntry(new Color(128, 0, 0), "maroon");
  public static final ColorEntry DEFVAL_RED = new ColorEntry(Color.red, "red");
  public static final ColorEntry DEFVAL_ORANGE = new ColorEntry(Color.orange, "orange");
  public static final ColorEntry DEFVAL_YELLOW = new ColorEntry(Color.yellow, "yellow");
  public static final ColorEntry DEFVAL_DARKGREEN = new ColorEntry(new Color(0, 128, 0), "dark green");
  public static final ColorEntry DEFVAL_GREEN = new ColorEntry(Color.green, "green");
  public static final ColorEntry DEFVAL_MAGENTA = new ColorEntry(Color.magenta, "magenta");
  public static final ColorEntry DEFVAL_PINK = new ColorEntry(Color.pink, "pink");
  // page 3 data: color map defaults (by id)
  private static final ColorEntry DEFVAL_COLOR1 = DEFVAL_BLACK;
  private static final ColorEntry DEFVAL_COLOR2 = DEFVAL_DARKGRAY;
  private static final ColorEntry DEFVAL_COLOR3 = DEFVAL_SLATEGRAY;
  private static final ColorEntry DEFVAL_COLOR4 = DEFVAL_LIGHTGRAY;
  private static final ColorEntry DEFVAL_COLOR5 = DEFVAL_NAVY;
  private static final ColorEntry DEFVAL_COLOR6 = DEFVAL_BLUE;
  private static final ColorEntry DEFVAL_COLOR7 = DEFVAL_CYAN;
  private static final ColorEntry DEFVAL_COLOR8 = DEFVAL_DARKGREEN;
  private static final ColorEntry DEFVAL_COLOR9 = DEFVAL_GREEN;
  private static final ColorEntry DEFVAL_COLORA = DEFVAL_MAROON;
  private static final ColorEntry DEFVAL_COLORB = DEFVAL_RED;
  private static final ColorEntry DEFVAL_COLORC = DEFVAL_ORANGE;
  private static final ColorEntry DEFVAL_COLORD = DEFVAL_YELLOW;
  private static final ColorEntry DEFVAL_COLORE = DEFVAL_MAGENTA;
  private static final ColorEntry DEFVAL_COLORF = DEFVAL_PINK;
  public static final List<ColorEntry> DEFVAL_COLORS = Arrays.asList(new ColorEntry[] { DEFVAL_COLOR1,
      DEFVAL_COLOR2, DEFVAL_COLOR3, DEFVAL_COLOR4, DEFVAL_COLOR5, DEFVAL_COLOR6, DEFVAL_COLOR7,
      DEFVAL_COLOR8, DEFVAL_COLOR9, DEFVAL_COLORA, DEFVAL_COLORB, DEFVAL_COLORC, DEFVAL_COLORD,
      DEFVAL_COLORE, DEFVAL_COLORF });
  public static final List<ColorEntry> DEFVAL_MAPPED_COLORS = new ArrayList<ColorEntry>();

  // page 4 data: color mappings
  private static final String KEY_COLOR_MAP = "colormap";
  public static final String KEY_NODE_COLOR1 = KEY_COLOR_MAP + ".1";
  public static final String KEY_NODE_COLOR2 = KEY_COLOR_MAP + ".2";
  public static final String KEY_NODE_COLOR3 = KEY_COLOR_MAP + ".3";
  public static final String KEY_NODE_COLOR4 = KEY_COLOR_MAP + ".4";
  public static final String KEY_NODE_COLOR5 = KEY_COLOR_MAP + ".5";
  public static final String KEY_NODE_COLOR6 = KEY_COLOR_MAP + ".6";
  public static final String KEY_NODE_COLOR7 = KEY_COLOR_MAP + ".7";
  public static final String KEY_NODE_COLOR8 = KEY_COLOR_MAP + ".8";
  public static final String KEY_NODE_COLOR9 = KEY_COLOR_MAP + ".9";
  public static final String KEY_NODE_COLORA = KEY_COLOR_MAP + ".A";
  public static final String KEY_NODE_COLORB = KEY_COLOR_MAP + ".B";
  public static final String KEY_NODE_COLORC = KEY_COLOR_MAP + ".C";
  public static final String KEY_NODE_COLORD = KEY_COLOR_MAP + ".D";
  public static final String KEY_NODE_COLORE = KEY_COLOR_MAP + ".E";
  public static final String KEY_NODE_COLORF = KEY_COLOR_MAP + ".F";
  public static final List<String> KEY_NODE_COLORS = Arrays.asList(new String[] { KEY_NODE_COLOR1,
      KEY_NODE_COLOR2, KEY_NODE_COLOR3, KEY_NODE_COLOR4, KEY_NODE_COLOR5, KEY_NODE_COLOR6, KEY_NODE_COLOR7,
      KEY_NODE_COLOR8, KEY_NODE_COLOR9, KEY_NODE_COLORA, KEY_NODE_COLORB, KEY_NODE_COLORC, KEY_NODE_COLORD,
      KEY_NODE_COLORE, KEY_NODE_COLORF });
  public static final String KEY_ASSIGNED_COLORS = KEY_COLOR_MAP + ".assignment";
  public static final List<ColorEntry> DEFVAL_ASSIGNED_COLORS = new ArrayList<ColorEntry>();

  private Map<String, Object> data;
  private Set<Integer> loadedPageStates;

  public ClusterWizardState() {

    data = new HashMap<String, Object>();
    loadedPageStates = new HashSet<Integer>();
  }

  private Map<String, Object> getData() {

    return data;
  }

  public void clear() {

    getData().clear();
    loadedPageStates.clear();
  }

  public Object getProperty(String key) {

    return getData().get(key);
  }

  public Set<String> getKeys() {

    return getData().keySet();
  }

  public void setProperty(String key, Object value) {

    getData().put(key, value);
  }

  @Override
  public String toString() {

    String result = "";
    for (String key : data.keySet()) {
      result += String.format("%s = %s", key, data.get(key));
    }
    return result;
  }

  public boolean isStateLoaded(int pageId) {

    return loadedPageStates .contains(pageId);
  }

  public void loadDefaultState(int pageId) {

    loadedPageStates.add(pageId);
    switch (pageId) {
      case 1:
        setProperty(KEY_NAME, DEFVAL_NAME);
        setProperty(KEY_TOTAL_NODES, DEFVAL_TOTAL_NODES);
        setProperty(KEY_ROOT_NODES, DEFVAL_ROOT_NODES);
        setProperty(KEY_COLOR_SCHEME, DEFVAL_COLOR_SCHEME);
        setProperty(KEY_ISOMORPH_FREE, DEFVAL_ISOMORPH_FREE);
        break;
      case 2:
        setProperty(KEY_CLASS_ANY, DEFVAL_CLASS_ANY);
        setProperty(KEY_CLASS_CONNECTED, DEFVAL_CLASS_CONNECTED);
        setProperty(KEY_CLASS_BICONNECTED, DEFVAL_CLASS_BICONNECTED);
        setProperty(KEY_CLASS_REEHOOVER, DEFVAL_CLASS_REEHOOVER);
        setProperty(KEY_EXCLUDE_NODAL_POINTS, DEFVAL_EXCLUDE_NODAL_POINTS);
        setProperty(KEY_EXCLUDE_ARTICULATION_POINTS, DEFVAL_EXCLUDE_ARTICULATION_POINTS);
        setProperty(KEY_EXCLUDE_ARTICULATION_PAIRS, DEFVAL_EXCLUDE_ARTICULATION_POINTS);
        break;
      case 3:
        for (int i = 0; i < KEY_COLORS.size(); i++) {
          setProperty(KEY_COLORS.get(i), DEFVAL_COLORS.get(i));
        }
        setProperty(KEY_MAPPED_COLORS, DEFVAL_MAPPED_COLORS);
        break;
      case 4:
        for (int i = 0; i < KEY_NODE_COLORS.size(); i++) {
          setProperty(KEY_NODE_COLORS.get(i), DEFVAL_COLORS.get(0));
        }
        setProperty(KEY_ASSIGNED_COLORS, DEFVAL_ASSIGNED_COLORS);
        break;
      case 5:
        break;
    }
  }
}