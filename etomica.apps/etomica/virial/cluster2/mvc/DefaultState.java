/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class DefaultState implements State {

  private Map<String, Object> data = new HashMap<String, Object>();

  private Map<String, Object> getData() {

    return data;
  }

  public void clear() {

    getData().clear();
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
    for (String key : data.keySet())  {
      result += String.format("%s = %s", key, data.get(key));
    }
    return result;
  }
}
