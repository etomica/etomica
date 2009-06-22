package etomica.api;

public interface IIntegratorListenerMD extends IIntegratorListener {

    /**
     * Invoked after the integrator has computed the forces on all atoms
     * in the system.
     * @param e
     */
    public void integratorForceComputed(IIntegratorEvent e);
    
}
