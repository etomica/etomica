package etomica.server.resources;

import etomica.meta.ComponentIndex;
import etomica.server.representations.SimClassInfo;
import etomica.simulation.Simulation;

import javax.inject.Inject;
import javax.ws.rs.GET;
import javax.ws.rs.Path;
import javax.ws.rs.Produces;
import javax.ws.rs.core.MediaType;
import java.util.List;
import java.util.stream.Collectors;

@Path("/simulations")
@Produces(MediaType.APPLICATION_JSON)
public class SimulationsIndexResource {

    private final ComponentIndex<Simulation> simIndex;

    @Inject
    public SimulationsIndexResource(ComponentIndex<Simulation> simIndex) {
        this.simIndex = simIndex;
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
}
