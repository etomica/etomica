package etomica.meta;

import etomica.meta.wrappers.Wrapper;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.Map;

import static etomica.meta.Common.REFLECTIONS;

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

        wrapperMap.put(Object.class, Wrapper.class.getDeclaredConstructors()[0]);

    }


    public static Wrapper getWrapper(Object o) {
        if (o == null) {
            return null;
        }
        return getWrapper(o, o.getClass());
    }

    @SuppressWarnings("unchecked")
    private static Wrapper getWrapper(Object o, Class cls) {
        if(wrapperMap.containsKey(cls)) {
            try {
                return (Wrapper) wrapperMap.get(cls).newInstance(o);
            } catch (InstantiationException | IllegalAccessException | InvocationTargetException e) {
                throw new RuntimeException(e);
            }
        } else {
            return getWrapper(o, cls.getSuperclass());
        }
    }

}
