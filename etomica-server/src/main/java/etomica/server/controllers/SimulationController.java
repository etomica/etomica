package etomica.server.controllers;

import etomica.meta.SimulationModel;
import etomica.server.EtomicaServer;
import etomica.simulation.Simulation;
import org.jooby.Response;
import org.jooby.Result;
import org.jooby.Results;
import org.jooby.mvc.*;

import java.util.Map;
import java.util.UUID;

@Path("/simulations")
public class SimulationController {

    @GET
    @Path("/:id")
    @Produces("application/json")
    public Result simulationStructure(UUID id) {
        SimulationModel model = EtomicaServer.simulations.get(id);
        return Results.json(model);
    }

    @POST
    @Consumes("application/json")
    public Result createSimulation(@Body SimulationConstructor params) {
        UUID simId = UUID.randomUUID();
        System.out.println(params.className);
        try {
            Simulation sim = (Simulation) Class.forName(params.className).newInstance();
            SimulationModel model = new SimulationModel(sim);
            EtomicaServer.simulations.put(simId, model);
            return Results.with(simId, 201);
        } catch (InstantiationException | IllegalAccessException | ClassNotFoundException e) {
            return Results.with(404);
        }
    }

    private static class SimulationConstructor {
        public String className;
        public Map<String, Object> constructorParams;
    }
}
