package etomica.virial.cluster2.graph;

import java.util.NoSuchElementException;

public interface EdgesGenerator extends Tagged {

  public boolean hasNext();

  public Edges next() throws NoSuchElementException;
}