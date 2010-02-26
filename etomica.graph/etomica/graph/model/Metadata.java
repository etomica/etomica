package etomica.graph.model;

import java.util.Arrays;
import java.util.List;

public interface Metadata extends Comparable<Metadata> {

  public static final char TYPE_EDGE_ANY = 'E';
  public static final char TYPE_NODE_ROOT = 'R';
  public static final char TYPE_NODE_FIELD = 'F';
  public static final char TYPE_NODE_DEFAULT = TYPE_NODE_FIELD;

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
  public static List<Character> COLOR_CODES = Arrays.asList(new Character[] { COLOR_CODE_0, COLOR_CODE_1,
      COLOR_CODE_2, COLOR_CODE_3, COLOR_CODE_4, COLOR_CODE_5, COLOR_CODE_6, COLOR_CODE_7, COLOR_CODE_8,
      COLOR_CODE_9 });

  // http://www.december.com/html/spec/colorsvg.html

  public static String COLOR_0 = "black"; // Color.BLACK;
  public static String COLOR_1 = "red"; // Color.RED;
  public static String COLOR_2 = "green"; // Color.GREEN;
  public static String COLOR_3 = "darkblue"; // Color.BLUE;
  public static String COLOR_4 = "darkviolet"; // Color.MAGENTA;
  public static String COLOR_5 = "royalblue"; // Color.CYAN;
  public static String COLOR_6 = "yellow"; // Color.YELLOW;
  public static String COLOR_7 = "darkorange"; // Color.ORANGE;
  public static String COLOR_8 = "fuchsia"; // Color.PINK;
  public static String COLOR_9 = "slategray"; // Color.LIGHT_GRAY;
  public static String[] COLORS = new String[] { COLOR_0, COLOR_1, COLOR_2, COLOR_3, COLOR_4, COLOR_5,
      COLOR_6, COLOR_7, COLOR_8, COLOR_9 };

  public char getColor();

  public char getType();

  public boolean isCompatible(Metadata other);

  public boolean isSameColor(Metadata other);

  public boolean isSameType(Metadata other);
}
