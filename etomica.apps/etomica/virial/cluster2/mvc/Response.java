package etomica.virial.cluster2.mvc;

import java.util.List;

public interface Response {

  public abstract State getData();

  public abstract List<MVCException> getErrors();
}