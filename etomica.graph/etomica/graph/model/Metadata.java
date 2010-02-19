package etomica.graph.model;

import java.awt.Color;

public interface Metadata extends Comparable<Metadata>{

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

  public char getColor();

  public char getType();

  public boolean isCompatible(Metadata other);

  public boolean isSameColor(Metadata other);

  public boolean isSameType(Metadata other);
}
