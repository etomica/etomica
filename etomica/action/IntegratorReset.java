package etomica.action;

import etomica.Integrator;


/**
 * Action that calls the reset method of an integrator.
  */

/*
 * History
 * Created on Feb 5, 2005 by kofke
 */
public class IntegratorReset extends IntegratorActionAdapter {

    /**
     * @param label
     */
    public IntegratorReset(Integrator integrator) {
        super("Reset integrator");
    }

    /* (non-Javadoc)
     * @see etomica.Action#actionPerformed()
     */
    public void actionPerformed() {
        if(integrator != null) integrator.reset();
    }

}
