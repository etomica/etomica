package etomica.action;

import etomica.integrator.Integrator;


/**
 * Action that calls the reset method of an integrator.
  */

/*
 * History
 * Created on Feb 5, 2005 by kofke
 */
public class IntegratorReset extends IntegratorActionAdapter {

    public IntegratorReset(Integrator integrator) {
        super("Reset integrator");
        this.integrator = integrator;
    }

    public void actionPerformed() {
        if(integrator != null) integrator.reset();
    }

}
