package etomica.meta;

import java.util.List;
import java.util.Set;

import static etomica.meta.Common.*;

public class ComponentIndex<T> {
    private final Class<T> componentClass;


    public ComponentIndex(Class<T> componentClass) {
        this.componentClass = componentClass;
    }

    public List<Class<?>> getComponentSet() {
        if (componentClass.isInterface()) {
            return CLASSPATH_SCAN.getClassesImplementing(componentClass.getName()).loadClasses();
        } else {
            return CLASSPATH_SCAN.getSubclasses(componentClass.getName()).loadClasses();
        }

    }
}
