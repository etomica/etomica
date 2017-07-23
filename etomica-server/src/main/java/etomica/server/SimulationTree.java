package etomica.server;

import etomica.meta.InstanceProperty;
import etomica.meta.wrappers.CollectionWrapper;
import etomica.meta.wrappers.Wrapper;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;


public class SimulationTree {
    private final List<SimulationTree> children = new ArrayList<>();
    protected Wrapper<?> wrapper;

    public SimulationTree(Wrapper<?> wrapper, Map<Class, List<InstanceProperty>> classes) {
        this.wrapper = wrapper;

        if (wrapper != null) {
            findChildren(classes);
        }
    }

    public Wrapper<?> getWrapper() {
        return wrapper;
    }

    public List<SimulationTree> getChildren() {
        return children;
    }

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
}
