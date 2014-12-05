/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.HashMap;
import java.util.Map;

public class ExcludeParameters implements Parameters {

  public Map<Character, char[]> bondMap = new HashMap<Character, char[]>();
}
