package etomica.server.dao;

import etomica.meta.SimulationModel;

import java.util.UUID;

public interface SimulationStore {
    public SimulationModel get(UUID id);
    public void put(UUID id, SimulationModel model);
}
