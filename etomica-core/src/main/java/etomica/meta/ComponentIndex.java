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
            return CLASSPATH_SCAN.classNamesToClassRefs(CLASSPATH_SCAN.getNamesOfClassesImplementing(componentClass));
        } else {
            return CLASSPATH_SCAN.classNamesToClassRefs(CLASSPATH_SCAN.getNamesOfSubclassesOf(componentClass));
        }

    }
}
