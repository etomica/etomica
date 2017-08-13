package etomica.server.dao;

import etomica.meta.SimulationModel;

import javax.inject.Singleton;
import java.util.UUID;
import java.util.concurrent.ConcurrentHashMap;

@Singleton
public class SimulationStore extends ConcurrentHashMap<UUID, SimulationModel> {

}
