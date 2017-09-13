package etomica.server.resources.data;

import etomica.data.*;
import etomica.meta.ComponentIndex;
import etomica.meta.DataSourceIndex;
import etomica.meta.SimulationModel;
import etomica.server.core.construction.Construction;
import etomica.server.dao.DataStreamStore;
import etomica.server.dao.SimulationStore;
import etomica.server.representations.ConstructionInfo;
import etomica.server.representations.ConstructionParams;
import org.apache.commons.beanutils.PropertyUtils;

import javax.inject.Inject;
import javax.validation.constraints.NotNull;
import javax.ws.rs.*;
import javax.ws.rs.core.MediaType;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.UUID;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@Path("/simulations/{simId}/data")
@Produces(MediaType.APPLICATION_JSON)
public class DataStreamResource {
    private final SimulationStore simStore;
    private final DataStreamStore dataStore;
    private final ComponentIndex<IEtomicaDataSource> dataSourceIndex;
    private final ComponentIndex<DataAccumulator> accumulatorIndex;

    @Inject
    public DataStreamResource(SimulationStore simStore, DataStreamStore dataStore, DataSourceIndex index) {
        this.simStore = simStore;
        this.dataStore = dataStore;
        this.dataSourceIndex = new ComponentIndex<>(IEtomicaDataSource.class);
        accumulatorIndex = new ComponentIndex<>(DataAccumulator.class);
    }


    @POST
    public UUID createDataStream(@NotNull ConstructionParams constructionParams, @PathParam("simId") String simId) {
        SimulationModel simModel = simStore.get(UUID.fromString(simId));
        UUID dataId = UUID.randomUUID();

        try {
            IEtomicaDataSource meter = Construction.<IEtomicaDataSource>createInstance(constructionParams, simModel)
                    .orElseThrow(() -> new WebApplicationException("Unable to create meter"));

            DataDump dump = new DataDump();
            DataPumpListener pump = new DataPumpListener(meter, dump, 100);
            this.dataStore.put(dataId, new DataStreamStore.DataPlumbing(pump, dump));
            return dataId;

        } catch (ClassNotFoundException e) {
            e.printStackTrace();
            throw new WebApplicationException("Unable to create meter");
        }
    }

    @GET
    @Path("meters")
    public List<ConstructionInfo> listMeters(@PathParam("simId") String simId) {
        SimulationModel model = this.simStore.get(UUID.fromString(simId));

        Set<Class> dataSources = dataSourceIndex.getComponentSet().stream()
                .filter(cls -> !IDataSink.class.isAssignableFrom(cls)).collect(Collectors.toSet());

        return dataSources.stream()
                .map(cls -> Construction.getConstructionInfo(cls, model))
                .flatMap(opt -> opt.map(Stream::of).orElseGet(Stream::empty))
                .collect(Collectors.toList());
    }

    @GET
    @Path("accumulators")
    public List<ConstructionInfo> listAccumulators(@PathParam("simId") String simId) {
        SimulationModel model = this.simStore.get(UUID.fromString(simId));

        Set<Class<? extends DataAccumulator>> accumulators = accumulatorIndex.getComponentSet();

        return accumulators.stream()
                .map(cls -> Construction.getConstructionInfo(cls, model))
                .flatMap(opt -> opt.map(Stream::of).orElseGet(Stream::empty))
                .collect(Collectors.toList());
    }
}
