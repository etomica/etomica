package etomica.virial.cluster2.mvc;

public interface WizardView extends View {

  public void attachPageView(WizardPageView page);

  public void detachPageView(WizardPageView page);

  public void close();
}