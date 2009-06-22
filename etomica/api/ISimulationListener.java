package etomica.api;

public interface ISimulationListener {

    public void simulationBoxAdded(ISimulationBoxEvent e);
    
    public void simulationBoxRemoved(ISimulationBoxEvent e);
    
    public void simulationSpeciesAdded(ISimulationSpeciesEvent e);
    
    public void simulationSpeciesRemoved(ISimulationSpeciesEvent e);
    
    public void simulationSpeciesIndexChanged(ISimulationSpeciesIndexEvent e);
    
    public void simulationSpeciesMaxIndexChanged(ISimulationIndexEvent e);
    
    public void simulationAtomTypeIndexChanged(ISimulationAtomTypeIndexEvent e);

    public void simulationAtomTypeMaxIndexChanged(ISimulationIndexEvent e);
}
