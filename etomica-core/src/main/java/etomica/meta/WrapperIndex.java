package etomica.meta;

import etomica.meta.wrappers.*;
import etomica.space.Vector;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import static etomica.meta.Common.REFLECTIONS;

/**
 * Provides a map between a class and the constructor of its appropriate wrapper type.
 */

//put all this in SimulationModel?
public class WrapperIndex {
    private static final Map<Class, Constructor> wrapperMap;

    static {
        wrapperMap = new HashMap<>();
        REFLECTIONS.getSubTypesOf(Wrapper.class).forEach(wrapperClass -> {
            wrapperMap.put(
                    wrapperClass.getDeclaredConstructors()[0].getParameterTypes()[0],
                    wrapperClass.getDeclaredConstructors()[0]
            );
        });

        wrapperMap.put(Object.class, ObjectWrapper.class.getDeclaredConstructors()[0]);


    }


    public static Wrapper getWrapper(Object o, SimulationModel simModel) {
        if (o == null) {
            return null;
        } else if(o.getClass().isArray()) {
            return new ArrayWrapper((Object[]) o, simModel);
        } else if(Collection.class.isAssignableFrom(o.getClass())) {
            return new CollectionWrapper((Collection) o, simModel);
        } else if(Vector.class.isAssignableFrom(o.getClass())) {    //would like a more extensible way to handle interface wrappers
            return new VectorWrapper((Vector) o, simModel);
        }
        return getWrapper(o, o.getClass(), simModel);
    }

    @SuppressWarnings("unchecked")
    private static Wrapper getWrapper(Object o, Class cls, SimulationModel simModel) {
        if(wrapperMap.containsKey(cls)) {
            try {
                return (Wrapper) wrapperMap.get(cls).newInstance(o, simModel);
            } catch (InstantiationException | IllegalAccessException | InvocationTargetException e) {
                throw new RuntimeException(e);
            }
        } else {
            return getWrapper(o, cls.getSuperclass(), simModel);
        }
    }

/* from SimulationTree
    protected void findChildren(Map<Class, List<InstanceProperty>> classes) {

        classes.putIfAbsent(this.wrapper.getWrappedClass(), this.wrapper.getProperties());

        if(wrapper instanceof CollectionWrapper) {
            ((CollectionWrapper<?>) wrapper).getElements().forEach(wrapper -> {
                children.add(new SimulationTree(wrapper, classes));
            });
        } else {
            wrapper.getChildren().forEach((propName, wrapper) -> {
                children.add(new PropertyTree(wrapper, propName, classes));
            });
        }
    }
*/

    //           System.out.printf("%s -> %s (%s)\n", this.wrappedClass.getSimpleName(), prop.getName(), prop.getPropertyType().getSimpleName());

//            if (prop.isIndexedProperty()) {
//                if (!prop.canCount()) {
//                    continue;
//                }
//
//                int count = prop.invokeCount();
//                Object[] propValues = new Object[count];
//                for (int i = 0; i < count; i++) {
//                    propValues[i] = prop.invokeReader(i);
//                }
//                children.put(prop.getName(), new ArrayWrapper(propValues));
//            } else if (List.class.isAssignableFrom(prop.getPropertyType())) {
//                List childList = (List) prop.invokeReader();
//                for (Object o : childList) {
//                    children.put(prop.getName(), WrapperIndex.getWrapper(o));
//                }
//            } else {
//                getWrapper(prop.invokeReader());
//            }
//   //             children.put(prop.getName(), WrapperIndex.getWrapper(prop.invokeReader()));
//            }

}
