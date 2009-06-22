package etomica.api;

public interface ISimulationEventManager {

    public void addListener(ISimulationListener listener);
    
    public void removeListener(ISimulationListener listener);
    
}
