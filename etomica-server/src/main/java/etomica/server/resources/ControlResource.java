package etomica.server.resources;

import etomica.action.activity.Controller;
import etomica.meta.SimulationModel;
import etomica.server.dao.SimulationStore;
import etomica.server.representations.StatusAction;

import javax.inject.Inject;
import javax.ws.rs.*;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;
import java.util.Map;
import java.util.UUID;

@Path("/simulations/{id}/control")
@Produces(MediaType.APPLICATION_JSON)
@Consumes(MediaType.APPLICATION_JSON)
public class ControlResource {

    private final SimulationStore simStore;

    @Inject
    public ControlResource(SimulationStore simStore) {
        this.simStore = simStore;
    }

    @PUT
    public void doCommand(@PathParam("id") String id, StatusAction action) {
        Controller controller = simStore.get(UUID.fromString(id)).getSimulation().getController();
        switch (action.getStatus()) {
            case START:
                if(!controller.isActive()) {
                    new Thread(controller::actionPerformed).start();
                } else if(controller.isPaused()) {
                    controller.unPause();
                } else {
                    throw new WebApplicationException(Response.Status.CONFLICT);
                }
                break;
            case PAUSE:
                if(controller.isActive() && !controller.isPaused()) {
                    controller.pause();
                } else {
                    throw new WebApplicationException(Response.Status.CONFLICT);
                }
                break;
            case RESET:
                break;

        }

    }


}
