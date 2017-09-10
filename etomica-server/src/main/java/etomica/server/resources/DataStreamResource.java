package etomica.server.resources;

import etomica.box.Box;
import etomica.data.DataDump;
import etomica.data.DataPump;
import etomica.data.IDataSource;
import etomica.data.IEtomicaDataSource;
import etomica.integrator.Integrator;
import etomica.meta.SimulationModel;
import etomica.server.dao.SimulationStore;
import etomica.server.representations.MeterConstructor;
import etomica.species.Species;
import org.apache.commons.beanutils.PropertyUtils;

import javax.inject.Inject;
import javax.validation.constraints.NotNull;
import javax.ws.rs.POST;
import javax.ws.rs.Path;
import javax.ws.rs.PathParam;
import javax.ws.rs.Produces;
import javax.ws.rs.core.MediaType;
import java.beans.PropertyDescriptor;
import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.UUID;

@Path("/simulations/{simId}/data")
@Produces(MediaType.APPLICATION_JSON)
public class DataStreamResource {
    private final SimulationStore simStore;
    private final SimulationModel simModel;
    private static final Set<Class> SIM_CONTEXT_CLASSES = new HashSet<>(Arrays.asList(
            Box.class,
            Integrator.class,
            Species.class
    ));

    @Inject
    public DataStreamResource(SimulationStore simStore, @PathParam("simId") String simId) {
        this.simStore = simStore;
        this.simModel = simStore.get(UUID.fromString(simId));
    }

    @POST
    public UUID createDataStream(@NotNull MeterConstructor constructionParams) {
        UUID dataId = UUID.randomUUID();

        try {
            IEtomicaDataSource meter = (IEtomicaDataSource) Class.forName(constructionParams.className).newInstance();
            PropertyDescriptor[] descriptors = PropertyUtils.getPropertyDescriptors(meter);
            for(PropertyDescriptor descriptor : descriptors) {
                if(descriptor.getWriteMethod() != null) {
                    Object param = constructionParams.paramsMap.get(descriptor.getName());
                    if(SIM_CONTEXT_CLASSES.stream().anyMatch(cls -> cls.isAssignableFrom(descriptor.getPropertyType()))) {
                        Object simObject = simModel.getWrapperById((Long) param).getWrapped();
                        descriptor.getWriteMethod().invoke(meter, simObject);
                    } else {
                        descriptor.getWriteMethod().invoke(meter, param);
                    }
                }
            }

            DataPump pump = new DataPump(meter, new DataDump());
        } catch(ClassNotFoundException | IllegalAccessException | InstantiationException e) {
            e.printStackTrace();
        } catch(InvocationTargetException e) {
            e.printStackTrace();
        }

        return dataId;
    }
}
