package etomica.action;

import etomica.Action;
import etomica.Integrator;

 /**
  * Elementary action performed on an integrator.
  */
public interface IntegratorAction extends Action {

    public void setIntegrator(Integrator integrator);
    
    public Integrator getIntegrator();

}