package etomica.action;

import etomica.api.IIntegrator;

 /**
  * Elementary action performed on an integrator.
  */
public interface IntegratorAction extends Action {

    public void setIntegrator(IIntegrator integrator);
    
    public IIntegrator getIntegrator();

}