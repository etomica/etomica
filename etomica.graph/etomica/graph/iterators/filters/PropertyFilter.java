package etomica.graph.iterators.filters;


import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;
import etomica.graph.property.Property;

public class PropertyFilter extends LocalFilter {

  private Property property;

  public PropertyFilter(GraphIterator iterator, Property property) {

    super(iterator);
    this.property = property;
  }

  @Override
  protected boolean accept(Graph g) {

    return property.check(g);
  }
}