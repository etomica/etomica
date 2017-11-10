package etomica.server.dao;

import etomica.meta.SimulationModel;

import javax.inject.Singleton;
import java.util.UUID;
import java.util.concurrent.ConcurrentHashMap;

@Singleton
public class SimulationStore extends ConcurrentHashMap<UUID, SimulationModel> {
    @Override
    public SimulationModel get(Object o) {
        if (o.equals(UUID.fromString("11111111-1111-1111-1111-111111111111"))) {
            return (SimulationModel) super.values().toArray()[0];
        }
        return super.get(o);
    }
}
