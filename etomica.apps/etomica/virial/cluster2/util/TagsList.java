package etomica.virial.cluster2.util;

import java.util.ArrayList;
import java.util.List;

public class TagsList extends ArrayList<String> {

  private static final long serialVersionUID = 5651120322293972959L;

  public TagsList(List<String> tags) {

    super(tags);
  }

  public TagsList() {

    super();
  }

  @Override
  public String toString() {

    String result = "[";
    for (int i = 0; i < size(); i++) {
      result += get(i);
      if (i < size() - 1) {
        result += " > ";
      }
    }
    return result + "]";
  }
}