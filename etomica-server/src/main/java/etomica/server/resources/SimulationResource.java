package etomica.server.resources;

import etomica.meta.SimulationModel;
import etomica.server.dao.SimulationStore;
import etomica.server.representations.SimulationConstructor;
import etomica.simulation.Simulation;

import javax.inject.Inject;
import javax.validation.constraints.NotNull;
import javax.ws.rs.*;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;
import java.net.URI;
import java.util.Map;
import java.util.UUID;

@Path("/simulations")
@Produces(MediaType.APPLICATION_JSON)
@Consumes(MediaType.APPLICATION_JSON)
public class SimulationResource {
    private final SimulationStore simStore;

    @Inject
    public SimulationResource(SimulationStore simStore) {
        this.simStore = simStore;
    }

    @GET
    @Path("{id}")
    public SimulationModel structure(@PathParam("id") String id) {
        UUID uuid = UUID.fromString(id);
        return simStore.get(uuid);
    }

    @POST
    public UUID createSimulation(@NotNull SimulationConstructor constructionParams) {
        UUID id = UUID.randomUUID();
        try {
            Simulation sim = (Simulation) Class.forName(constructionParams.className).newInstance();
            simStore.put(id, new SimulationModel(sim));
            return id;
        } catch (IllegalAccessException | InstantiationException | ClassNotFoundException e) {
            throw new ServerErrorException("No simulation instance with that id", Response.Status.NOT_FOUND);
        }
    }


}
