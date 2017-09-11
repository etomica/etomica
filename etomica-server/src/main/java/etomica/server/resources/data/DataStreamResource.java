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
import etomica.server.dao.DataStreamStore;
import etomica.server.dao.SimulationStore;
import etomica.server.representations.MeterConstructor;
import etomica.server.representations.MeterInfo;
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

@Path("/simulations/{simId}/data")
@Produces(MediaType.APPLICATION_JSON)
public class DataStreamResource {
    private static final Set<Class> SIM_CONTEXT_CLASSES = new HashSet<>(Arrays.asList(
            Box.class,
            Integrator.class,
            Species.class,
            Space.class,
            PotentialMaster.class,
            MCMove.class,
            IRandom.class
    ));
    private final SimulationStore simStore;
    private final DataStreamStore dataStore;
    private final DataSourceIndex index;

    @Inject
    public DataStreamResource(SimulationStore simStore, DataStreamStore dataStore, DataSourceIndex index) {
        this.simStore = simStore;
        this.dataStore = dataStore;
        this.index = index;
    }

    static Optional<Constructor> getMeterConstructor(Class meterClass) {
        Constructor[] constructors = meterClass.getConstructors();
        List<Constructor> eligibleConstructors = new ArrayList<>();
        constructorLoop:
        for (Constructor constructor : constructors) {
            for (Class paramClass : constructor.getParameterTypes()) {
                if (!containsAssignableFrom(SIM_CONTEXT_CLASSES, paramClass)) {
                    continue constructorLoop;
                }
            }
            eligibleConstructors.add(constructor);
        }

        eligibleConstructors.sort(Comparator.comparingInt(Constructor::getParameterCount));
        if (eligibleConstructors.size() == 0) {
            return Optional.empty();
        } else {
            return Optional.of(eligibleConstructors.get(0));
        }
    }

    private static boolean containsAssignableFrom(Set<Class> classes, Class cls) {
        return classes.stream().anyMatch(c -> c.isAssignableFrom(cls));
    }

    @POST
    public UUID createDataStream(@NotNull MeterConstructor constructionParams, @PathParam("simId") String simId) {
        SimulationModel simModel = simStore.get(UUID.fromString(simId));
        UUID dataId = UUID.randomUUID();

        try {
            Constructor constructor = getMeterConstructor(Class.forName(constructionParams.className)).orElse(null);
            Object[] params = constructionParams.constructorParams.stream()
                    .map(o -> simModel.getWrapperById((Long) o)).toArray();

            IEtomicaDataSource meter = (IEtomicaDataSource) constructor.newInstance(params);

            constructionParams.paramsMap.forEach((propName, prop) -> {
                try {
                    if(containsAssignableFrom(SIM_CONTEXT_CLASSES, PropertyUtils.getPropertyType(meter, propName))) {
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
    public List<MeterInfo> list(@PathParam("simId") String simId) {
        SimulationModel model = this.simStore.get(UUID.fromString(simId));

        Set<Class> dataSources = index.getComponentSet().stream()
                .filter(cls -> !IDataSink.class.isAssignableFrom(cls)).collect(Collectors.toSet());

        List<MeterInfo> meterInfos = new ArrayList<>();
        // use flatMap
        meterInfosLoop:
        for (Class dataSource : dataSources) {
            Constructor constructor = getMeterConstructor(dataSource).orElse(null);
            if (constructor != null) {
                List<String> constructorParamTypes = Arrays.stream(constructor.getParameterTypes())
                        .map(Class::getCanonicalName)
                        .collect(Collectors.toList());

                Map<String, List<Long>> classOptions = new HashMap<>();
                for (Class paramClass : constructor.getParameterTypes()) {
                    List<Long> ids = model.getAllIdsOfType(paramClass);
                    if (ids.isEmpty()) {
                        continue meterInfosLoop;
                    }
                    classOptions.put(
                            paramClass.getCanonicalName(),
                            ids
                    );
                }

                PropertyDescriptor[] propertyDescriptors = PropertyUtils.getPropertyDescriptors(dataSource);
                Map<String, String> propertyMap = new HashMap<>();
                if (propertyDescriptors != null && propertyDescriptors.length != 0) {
                    Arrays.stream(propertyDescriptors)
                            .filter(propertyDescriptor -> {
                                return propertyDescriptor.getWriteMethod() != null &&
                                        !(propertyDescriptor instanceof IndexedPropertyDescriptor);
                            })
                            .filter(propertyDescriptor -> {
                                return propertyDescriptor.getPropertyType().isPrimitive()
                                        || propertyDescriptor.getPropertyType().equals(String.class)
                                        || containsAssignableFrom(SIM_CONTEXT_CLASSES, propertyDescriptor.getPropertyType());
                            })
                            .forEach(propertyDescriptor -> {
                                propertyMap.put(propertyDescriptor.getName(), propertyDescriptor.getPropertyType().getCanonicalName());
                                if(containsAssignableFrom(SIM_CONTEXT_CLASSES, propertyDescriptor.getPropertyType())) {
                                    classOptions.put(
                                            propertyDescriptor.getPropertyType().getCanonicalName(),
                                            model.getAllIdsOfType(propertyDescriptor.getPropertyType())
                                    );
                                }
                            });
                }

                MeterInfo info = new MeterInfo(
                        dataSource.getCanonicalName(),
                        constructorParamTypes,
                        classOptions,
                        propertyMap
                );

                meterInfos.add(info);
            }
        }
        return meterInfos;
    }
}
