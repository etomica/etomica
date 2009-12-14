package etomica.virial.cluster2.mvc;

public interface Action {

  public ActionResponse execute(final ViewResponse response, final State state);
}