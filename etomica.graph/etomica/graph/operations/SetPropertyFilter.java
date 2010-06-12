package etomica.graph.operations;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.property.Property;

public class SetPropertyFilter implements Unary {

  protected final Property property;

  public SetPropertyFilter(Property property) {
    this.property = property;
  }

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {
    assert (params == null);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      if (property.check(g)) {
        result.add(g);
      }
    }
    return result;
  }

}
