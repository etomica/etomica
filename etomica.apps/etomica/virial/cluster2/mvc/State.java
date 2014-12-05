/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc;

import java.util.Set;

public interface State {

  public void clear();

  public Set<String> getKeys();

  public Object getProperty(String key);

  public void setProperty(String key, Object value);
}