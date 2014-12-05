/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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