package etomica.server.resources;

import etomica.meta.SimulationModel;
import etomica.meta.properties.Property;
import etomica.meta.wrappers.Wrapper;
import etomica.server.dao.SimulationStore;
import etomica.server.representations.PropertyUpdate;
import io.dropwizard.jersey.PATCH;
import org.apache.commons.beanutils.BeanUtils;
import org.apache.commons.beanutils.PropertyUtils;

import javax.inject.Inject;
import javax.ws.rs.*;
import javax.ws.rs.core.MediaType;
import java.util.List;
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
                .findFirst().get();
        model.getSimulation().getController().doActionNow(() -> {
            wrapperProp.invokeWriter(propUpdate.getNewValue());
        });
    }

}
