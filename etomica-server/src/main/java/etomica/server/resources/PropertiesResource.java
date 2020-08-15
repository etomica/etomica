package etomica.server.resources;

import etomica.meta.SimulationModel;
import etomica.meta.properties.Property;
import etomica.meta.wrappers.Wrapper;
import etomica.server.dao.SimulationStore;
import etomica.server.representations.PropertyUpdate;

import javax.inject.Inject;
import javax.ws.rs.*;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;
import java.util.UUID;

@Path("/simulations/{simId}/properties")
@Consumes(MediaType.APPLICATION_JSON)
@Produces(MediaType.APPLICATION_JSON)
public class PropertiesResource {
    private final SimulationStore simStore;

    @Inject
    public PropertiesResource(SimulationStore simStore) {
        this.simStore = simStore;
    }

    @PUT
    public void updateProperty(@PathParam("simId") String simId, PropertyUpdate propUpdate) {
        SimulationModel model = simStore.get(UUID.fromString(simId));
        Wrapper<?> wrapper = model.getWrapperById(propUpdate.getId());
        Property wrapperProp = wrapper.getValueProperties().stream()
                .filter(p -> p.getName().equalsIgnoreCase(propUpdate.getProperty()))
                .findFirst().orElseThrow(() -> new WebApplicationException(Response.Status.BAD_REQUEST));
        model.getSimulation().getController().submitActionInterrupt(() -> {
            wrapperProp.invokeWriter(propUpdate.getNewValue());
        });
    }

}
