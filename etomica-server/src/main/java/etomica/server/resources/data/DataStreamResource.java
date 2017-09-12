package etomica.server.resources.data;

import etomica.box.Box;
import etomica.data.DataDump;
import etomica.data.DataPumpListener;
import etomica.data.IDataSink;
import etomica.data.IEtomicaDataSource;
import etomica.integrator.Integrator;
import etomica.integrator.mcmove.MCMove;
import etomica.meta.DataSourceIndex;
import etomica.meta.SimulationModel;
import etomica.potential.PotentialMaster;
import etomica.server.core.construction.Construction;
import etomica.server.dao.DataStreamStore;
import etomica.server.dao.SimulationStore;
import etomica.server.representations.ConstructionParams;
import etomica.server.representations.ConstructionInfo;
import etomica.space.Space;
import etomica.species.Species;
import etomica.util.random.IRandom;
import org.apache.commons.beanutils.PropertyUtils;

import javax.inject.Inject;
import javax.validation.constraints.NotNull;
import javax.ws.rs.*;
import javax.ws.rs.core.MediaType;
import java.beans.IndexedPropertyDescriptor;
import java.beans.PropertyDescriptor;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@Path("/simulations/{simId}/data")
@Produces(MediaType.APPLICATION_JSON)
public class DataStreamResource {
    private final SimulationStore simStore;
    private final DataStreamStore dataStore;
    private final DataSourceIndex index;

    @Inject
    public DataStreamResource(SimulationStore simStore, DataStreamStore dataStore, DataSourceIndex index) {
        this.simStore = simStore;
        this.dataStore = dataStore;
        this.index = index;
    }


    @POST
    public UUID createDataStream(@NotNull ConstructionParams constructionParams, @PathParam("simId") String simId) {
        SimulationModel simModel = simStore.get(UUID.fromString(simId));
        UUID dataId = UUID.randomUUID();

        try {
            Constructor constructor = Construction.getEligibleConstructor(Class.forName(constructionParams.className)).orElse(null);
            Object[] params = constructionParams.constructorParams.stream()
                    .map(o -> simModel.getWrapperById((Long) o)).toArray();

            IEtomicaDataSource meter = (IEtomicaDataSource) constructor.newInstance(params);

            constructionParams.paramsMap.forEach((propName, prop) -> {
                try {
                    if(Construction.isContextClass(PropertyUtils.getPropertyType(meter, propName))) {
                        PropertyUtils.setSimpleProperty(
                                meter,
                                propName,
                                simModel.getWrapperById(((Integer)prop).longValue()).getWrapped()
                        );
                    } else {
                        PropertyUtils.setSimpleProperty(meter, propName, prop);
                    }
                } catch (IllegalAccessException | InvocationTargetException | NoSuchMethodException e) {
                    e.printStackTrace();
                }
            });

            DataDump dump = new DataDump();
            DataPumpListener pump = new DataPumpListener(meter, dump, 100);
            this.dataStore.put(dataId, new DataStreamStore.DataPlumbing(pump, dump));

        } catch (ClassNotFoundException | IllegalAccessException | InstantiationException e) {
            e.printStackTrace();
        } catch (InvocationTargetException e) {
            e.printStackTrace();
        }

        return dataId;
    }

    @GET
    public List<ConstructionInfo> list(@PathParam("simId") String simId) {
        SimulationModel model = this.simStore.get(UUID.fromString(simId));

        Set<Class> dataSources = index.getComponentSet().stream()
                .filter(cls -> !IDataSink.class.isAssignableFrom(cls)).collect(Collectors.toSet());

        return dataSources.stream()
                .map(cls -> Construction.getConstructionInfo(cls, model))
                .flatMap(opt -> opt.map(Stream::of).orElseGet(Stream::empty))
                .collect(Collectors.toList());
    }
}
