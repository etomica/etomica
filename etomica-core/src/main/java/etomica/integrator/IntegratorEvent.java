package etomica.integrator;

public class IntegratorEvent {

    private final Integrator integrator;

    public IntegratorEvent(Integrator integrator) {
        this.integrator = integrator;
    }

    public Integrator getIntegrator() {
        return integrator;
    }
}
