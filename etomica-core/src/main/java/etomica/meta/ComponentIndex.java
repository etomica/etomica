package etomica.meta;

import org.reflections.Reflections;

import java.util.Set;

public class ComponentIndex<T> {
    private final Class<T> componentClass;

    private static final Reflections reflections = new Reflections("etomica");

    public ComponentIndex(Class<T> componentClass) {
        this.componentClass = componentClass;
    }

    public Set<Class<? extends T>> getComponentSet() {
        return reflections.getSubTypesOf(componentClass);
    }
}
