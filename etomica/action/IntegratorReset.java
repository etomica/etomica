package etomica.action;

import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IIntegrator;


/**
 * Action that calls the reset method of an integrator.
  */
public class IntegratorReset extends IntegratorActionAdapter {

    public IntegratorReset() {
        super();
    }
    
    public IntegratorReset(IIntegrator integrator, boolean ignoreOverlap) {
        super(integrator);
        this.ignoreOverlap = ignoreOverlap;
    }

    public void actionPerformed() {
        if(integrator != null) {
            try {
                integrator.reset();
            }
            catch (ConfigurationOverlapException e) {
                if (!ignoreOverlap) {
                    throw new RuntimeException(e);
                }
            }
        }
    }

    /**
     * @return Returns the ignoreOverlap.
     */
    public boolean isIgnoreOverlap() {
        return ignoreOverlap;
    }

    /**
     * @param ignoreOverlap The ignoreOverlap to set.
     */
    public void setIgnoreOverlap(boolean ignoreOverlap) {
        this.ignoreOverlap = ignoreOverlap;
    }

    private static final long serialVersionUID = 1L;
    private boolean ignoreOverlap;
}
