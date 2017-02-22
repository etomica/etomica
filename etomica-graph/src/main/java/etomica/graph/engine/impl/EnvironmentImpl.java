/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine.impl;

import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import etomica.graph.engine.Environment;
import etomica.graph.model.Graph;

public class EnvironmentImpl implements Environment {

  private SortedMap<String, String> props = new TreeMap<String, String>();
  private SortedMap<String, Set<Graph>> vars = new TreeMap<String, Set<Graph>>();

  public String getPropertyValue(String propName) {

    return props.get(propName);
  }

  public Set<Graph> getVariableValue(String varName) {

    return vars.get(varName);
  }

  public void setProperty(String propName, String propValue) {

    if (propValue == null) {
      props.remove(propName);
    }
    props.put(propName, propValue);
  }

  public void setVariable(String varName, Set<Graph> varValue) {

    if (varValue == null) {
      vars.remove(varName);
    }
    vars.put(varName, varValue);
  }

  public Set<String> getPropertyNames() {

    return props.keySet();
  }

  public Set<String> getVariableNames() {

    return vars.keySet();
  }
}