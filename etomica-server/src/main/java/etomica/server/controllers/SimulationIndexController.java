package etomica.server.controllers;

import etomica.meta.ComponentIndex;
import etomica.meta.SimulationModel;
import etomica.server.EtomicaServer;
import etomica.simulation.Simulation;
import org.jooby.Result;
import org.jooby.Results;
import org.jooby.mvc.GET;
import org.jooby.mvc.POST;
import org.jooby.mvc.Path;
import org.jooby.mvc.Produces;

import java.util.UUID;

@Path("/simulations")
public class SimulationIndexController {
    private static final ComponentIndex<Simulation> simIndex = new ComponentIndex<>(Simulation.class);

    @GET
    @Produces("application/json")
    public Result list() {
        return Results.json(simIndex.getComponentSet().stream()
                .filter(cls -> {
                    try {
                        cls.getConstructor();
                        return true;
                    } catch (NoSuchMethodException e) {
                        return false;
                    }
                })
                .map(Class::getCanonicalName).toArray());
    }



}
