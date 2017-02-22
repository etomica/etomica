package org.jmol.g3d;

import java.util.Comparator;

class TextSorter implements Comparator<TextString> {

  public int compare(TextString a, TextString b) {
    return (a == null || b == null ? 0 : a.z > b.z ? -1 : a.z < b.z ? 1 : 0);
  }

}
