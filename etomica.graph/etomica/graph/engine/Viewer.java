package etomica.graph.engine;

import java.util.Set;

import etomica.graph.model.Graph;

public interface Viewer {

  public void close();

  public void open();

  public void update(Set<Graph> graphs);
}