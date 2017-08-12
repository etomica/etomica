package etomica.serverdw.resources;

import etomica.meta.ComponentIndex;
import etomica.simulation.Simulation;

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

    public SimulationsIndexResource(ComponentIndex<Simulation> simIndex) {
        this.simIndex = simIndex;
    }

    @GET
    public List<String> list() {
        return simIndex.getComponentSet().stream()
                .filter(cls -> {
                    try {
                        cls.getConstructor();
                        return true;
                    } catch (NoSuchMethodException e) {
                        return false;
                    }
                })
                .map(Class::getCanonicalName).collect(Collectors.toList());
    }
}
