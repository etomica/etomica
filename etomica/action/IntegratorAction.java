package etomica.action;

import etomica.integrator.IIntegrator;

 /**
  * Elementary action performed on an integrator.
  */
public interface IntegratorAction extends Action {

    public void setIntegrator(IIntegrator integrator);
    
    public IIntegrator getIntegrator();

}