package etomica.server.core.construction;

import etomica.box.Box;
import etomica.integrator.Integrator;
import etomica.integrator.mcmove.MCMove;
import etomica.meta.SimulationModel;
import etomica.potential.PotentialMaster;
import etomica.server.representations.ConstructionInfo;
import etomica.server.representations.ConstructionParams;
import etomica.space.Space;
import etomica.species.Species;
import etomica.util.random.IRandom;
import org.apache.commons.beanutils.PropertyUtils;

import java.beans.IndexedPropertyDescriptor;
import java.beans.PropertyDescriptor;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.*;
import java.util.stream.Collectors;

public class Construction {
    public static final Set<Class<?>> SIM_CONTEXT_CLASSES = Collections.unmodifiableSet(new HashSet<>(Arrays.asList(
            Box.class,
            Integrator.class,
            Species.class,
            Space.class,
            PotentialMaster.class,
            MCMove.class,
            IRandom.class
    )));

    public static boolean isContextClass(Class cls) {
        return SIM_CONTEXT_CLASSES.stream().anyMatch(c -> c.isAssignableFrom(cls));
    }

    public static Optional<Constructor<?>> getEligibleConstructor(Class meterClass) {
        Constructor<?>[] constructors = meterClass.getConstructors();
        return Arrays.stream(constructors)
                .filter(constructor -> Arrays.stream(constructor.getParameterTypes()).allMatch(Construction::isContextClass))
                .sorted(Comparator.comparingInt(Constructor::getParameterCount))
                .findFirst();
    }

    public static Optional<ConstructionInfo> getConstructionInfo(Class cls, SimulationModel model) {
        return getEligibleConstructor(cls).flatMap(constructor -> {
            List<String> constructorParamTypes = Arrays.stream(constructor.getParameterTypes())
                    .map(Class::getCanonicalName)
                    .collect(Collectors.toList());

            Map<String, List<Long>> classOptions = new HashMap<>();
            for(Class paramClass : constructor.getParameterTypes()) {
                List<Long> ids = model.getAllIdsOfType(paramClass);
                if(ids.isEmpty()) {
                    return Optional.empty();
                } else {
                    classOptions.put(paramClass.getCanonicalName(), ids);
                }
            }

            PropertyDescriptor[] descriptors = PropertyUtils.getPropertyDescriptors(cls);
            Map<String, String> propertyMap = new HashMap<>();
            Arrays.stream(descriptors)
                    .filter(desc -> !(desc instanceof IndexedPropertyDescriptor))
                    .filter(desc -> desc.getWriteMethod() != null)
                    .filter(desc -> desc.getPropertyType().isPrimitive() || desc.getPropertyType().equals(String.class) || isContextClass(desc.getPropertyType()))
                    .forEach(desc -> {
                        propertyMap.put(desc.getName(), desc.getPropertyType().getCanonicalName());
                        if(isContextClass(desc.getPropertyType())) {
                            classOptions.put(
                                    desc.getPropertyType().getCanonicalName(),
                                    model.getAllIdsOfType(desc.getPropertyType())
                            );
                        }
                    });

            return Optional.of(new ConstructionInfo(
                    cls.getCanonicalName(),
                    constructorParamTypes,
                    classOptions,
                    propertyMap
            ));
        });
    }

    public static <T> Optional<T> createInstance(ConstructionParams params, SimulationModel model) throws ClassNotFoundException {
        return getEligibleConstructor(Class.forName(params.className)).flatMap(constructor -> {
            Object[] constructorParams = params.constructorParams.stream()
                    .map(o -> model.getWrapperById(((Integer) o).longValue()).getWrapped()).toArray();

            try {
                T obj = (T) constructor.newInstance(constructorParams);
                for(Map.Entry<String, Object> prop : params.paramsMap.entrySet()) {

                    if(isContextClass(PropertyUtils.getPropertyType(obj, prop.getKey()))) {
                        PropertyUtils.setSimpleProperty(
                                obj,
                                prop.getKey(),
                                model.getWrapperById(((Integer) prop.getValue()).longValue()).getWrapped()
                        );
                    } else {
                        PropertyUtils.setSimpleProperty(obj, prop.getKey(), prop.getValue());
                    }

                }

                return Optional.of(obj);
            } catch (InstantiationException | IllegalAccessException | InvocationTargetException | NoSuchMethodException e) {
                e.printStackTrace();
                return Optional.empty();
            }
        });
    }


}
