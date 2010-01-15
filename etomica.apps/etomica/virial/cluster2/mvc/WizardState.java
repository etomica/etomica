package etomica.virial.cluster2.mvc;

public interface WizardState extends State {

  // sets the default state for the page
  public void loadDefaultState(int pageId);

  // determines if the state for the given page has been previously loaded
  public boolean isStateLoaded(int pageId);
}