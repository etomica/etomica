package etomica.virial.cluster2.mvc;

public interface ActionResponse extends Response {

  public Action getAction();

  public ActionStatus getStatus();
}