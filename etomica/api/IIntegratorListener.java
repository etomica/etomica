package etomica.api;

public interface IIntegratorListener {

    /**
     * Invoked when integration begins.
     * @param e
     */
    public void integratorInitialized(IIntegratorEvent e);
   
    /**
     * Invoked at the beginning of each integrator step.
     * @param e
     */
    public void integratorStepStarted(IIntegratorEvent e);
    
    /**
     * Invoked at the end of each integrator step.
     * @param e
     */
    public void integratorStepFinished(IIntegratorEvent e);
    
}
