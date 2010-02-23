package etomica.graph.engine;

import java.util.Set;

import etomica.graph.model.Graph;

public interface Environment {

  public Set<String> getPropertyNames();

  public String getPropertyValue(String propName);

  public Set<String> getVariableNames();

  public Set<Graph> getVariableValue(String varName);

  public void setProperty(String propName, String value);

  public void setVariable(String varName, Set<Graph> value);
}