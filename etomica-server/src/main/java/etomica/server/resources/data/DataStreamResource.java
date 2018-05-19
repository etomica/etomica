package etomica.server.resources.data;

import etomica.data.*;
import etomica.meta.ComponentIndex;
import etomica.meta.DataSourceIndex;
import etomica.meta.SimulationModel;
import etomica.server.core.construction.Construction;
import etomica.server.dao.DataStreamStore;
import static etomica.server.dao.DataStreamStore.DataPlumbing;
import etomica.server.dao.SimulationStore;
import etomica.server.representations.ConstructionInfo;
import etomica.server.representations.ConstructionParams;
import io.dropwizard.jersey.params.UUIDParam;

import javax.inject.Inject;
import javax.validation.constraints.NotNull;
import javax.ws.rs.*;
import javax.ws.rs.core.MediaType;
import java.lang.reflect.Modifier;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.UUID;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@Path("/simulations/{simId}/data")
@Produces(MediaType.APPLICATION_JSON)
public class DataStreamResource {
    private final SimulationStore simStore;
    private final DataStreamStore dataStore;
    private final ComponentIndex<IDataSource> dataSourceIndex;
    private final ComponentIndex<DataAccumulator> accumulatorIndex;

    @Inject
    public DataStreamResource(SimulationStore simStore, DataStreamStore dataStore, DataSourceIndex index) {
        this.simStore = simStore;
        this.dataStore = dataStore;
        this.dataSourceIndex = new ComponentIndex<>(IDataSource.class);
        accumulatorIndex = new ComponentIndex<>(DataAccumulator.class);
    }

    @GET
    public List<UUID> getStreamsForSim(@PathParam("simId") String simId) {
        UUID id = UUID.fromString(simId);
        return dataStore.entrySet().stream()
                .filter(entry -> entry.getValue().getSimId().equals(id))
                .map(Map.Entry::getKey)
                .collect(Collectors.toList());
    }


    @POST
    public UUID createDataStream(@NotNull ConstructionParams constructionParams, @PathParam("simId") String simId) {
        SimulationModel simModel = simStore.get(UUID.fromString(simId));
        UUID dataId = UUID.randomUUID();

        try {
            IDataSource meter = Construction.<IDataSource>createInstance(constructionParams, simModel)
                    .orElseThrow(() -> new WebApplicationException("Unable to create meter"));

            DataDump dump = new DataDump();
            DataPumpListener pump = new DataPumpListener(meter, dump, 100);
            this.dataStore.put(dataId, new DataStreamStore.DataPlumbing(pump, dump, UUID.fromString(simId)));
            return dataId;

        } catch (ClassNotFoundException e) {
            e.printStackTrace();
            throw new WebApplicationException("Unable to create meter");
        }
    }

    @POST
    @Path("{dataId}/accumulator")
    public void addAccumulator(@NotNull ConstructionParams params, @PathParam("simId") String simID, @PathParam("dataId") String dataId) {
        SimulationModel model = simStore.get(UUID.fromString(simID));

        try {
            DataAccumulator accumulator = Construction.<DataAccumulator>createInstance(params, model)
                    .orElseThrow(() -> new WebApplicationException("Unable to create accumulator"));

            DataPlumbing plumbing = this.dataStore.get(UUID.fromString(dataId));
            DataPump pump = plumbing.getPump();
            DataDump dump = plumbing.getDump();
            pump.setDataSink(accumulator);
            accumulator.addDataSink(dump);
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
            throw new WebApplicationException("Unable to create accumulator");
        }
    }

    @GET
    @Path("meters")
    public List<ConstructionInfo> listMeters(@PathParam("simId") String simId) {
        SimulationModel model = this.simStore.get(UUID.fromString(simId));

        Set<Class<?>> dataSources = dataSourceIndex.getComponentSet().stream()
                .filter(cls -> !IDataSink.class.isAssignableFrom(cls))
                .filter(cls -> !Modifier.isInterface(cls.getModifiers()))
                .filter(cls -> !cls.isInterface())
                .collect(Collectors.toSet());

        return dataSources.stream()
                .map(cls -> Construction.getConstructionInfo(cls, model))
                .flatMap(opt -> opt.map(Stream::of).orElseGet(Stream::empty))
                .collect(Collectors.toList());
    }

    @GET
    @Path("accumulators")
    public List<ConstructionInfo> listAccumulators(@PathParam("simId") String simId) {
        SimulationModel model = this.simStore.get(UUID.fromString(simId));

        Set<Class<?>> accumulators = accumulatorIndex.getComponentSet().stream()
                .filter(cls -> !Modifier.isInterface(cls.getModifiers()))
                .filter(cls -> !cls.isInterface())
                .collect(Collectors.toSet());

        return accumulators.stream()
                .map(cls -> Construction.getConstructionInfo(cls, model))
                .flatMap(opt -> opt.map(Stream::of).orElseGet(Stream::empty))
                .collect(Collectors.toList());
    }

    @GET
    @Path("{dataId}")
    public IDataInfo getDataInfo(@PathParam("simId") String simId, @PathParam("dataId") String dataId) {
        //TODO: serializers
        return this.dataStore.get(UUID.fromString(dataId)).getDump().getDataInfo();
    }
}
