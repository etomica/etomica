package etomica.server.resources;

import etomica.meta.ComponentIndex;
import etomica.server.dao.SimulationStore;
import etomica.server.representations.SimClassInfo;
import etomica.simulation.Simulation;

import javax.inject.Inject;
import javax.ws.rs.GET;
import javax.ws.rs.Path;
import javax.ws.rs.Produces;
import javax.ws.rs.core.MediaType;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

@Path("/simulations")
@Produces(MediaType.APPLICATION_JSON)
public class SimulationsIndexResource {

    private final ComponentIndex<Simulation> simIndex;
    private final SimulationStore simStore;

    @Inject
    public SimulationsIndexResource(ComponentIndex<Simulation> simIndex, SimulationStore simStore) {
        this.simIndex = simIndex;
        this.simStore = simStore;
    }

    @GET
    public List<SimClassInfo> list() {
        return simIndex.getComponentSet().stream()
                .filter(cls -> {
                    try {
                        cls.getConstructor();
                        return true;
                    } catch (NoSuchMethodException e) {
                        return false;
                    }
                })
                .map(SimClassInfo::forClass)
                .collect(Collectors.toList());
    }

    @GET
    @Path("/instances")
    public List<Map<String, Object>> getInstances() {
        return this.simStore.entrySet().stream().map((e) -> {
            Map<String, Object> map = new HashMap<>();
            map.put("id", e.getKey());
            map.put("classInfo", SimClassInfo.forClass(e.getValue().getSimulation().getClass()));
            return map;
        }).collect(Collectors.toList());
    }
}
