package etomica.meta;

import org.reflections.Reflections;

import java.util.Set;

import static etomica.meta.Common.REFLECTIONS;

public class ComponentIndex<T> {
    private final Class<T> componentClass;


    public ComponentIndex(Class<T> componentClass) {
        this.componentClass = componentClass;
    }

    public Set<Class<? extends T>> getComponentSet() {
        return REFLECTIONS.getSubTypesOf(componentClass);
    }
}
