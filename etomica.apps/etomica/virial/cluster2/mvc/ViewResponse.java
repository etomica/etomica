package etomica.virial.cluster2.mvc;

public interface ViewResponse extends Response {

  public View getView();

  public ViewStatus getStatus();
}