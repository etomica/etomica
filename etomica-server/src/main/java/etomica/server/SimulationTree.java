package etomica.server;

import etomica.meta.InstanceProperty;
import etomica.meta.wrappers.Wrapper;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;


public class SimulationTree {
    private final List<SimulationTree> children = new ArrayList<>();
    private Wrapper<?> data;
    private String propName;

    public SimulationTree(Wrapper<?> data, String propName, Map<Class, List<InstanceProperty>> classes) {
        this.data = data;
        this.propName = propName;

        if (data != null) {
            findChildren(classes);
        }
    }

    public String getPropName() {
        return propName;
    }

    public Wrapper<?> getData() {
        return data;
    }

    public List<SimulationTree> getChildren() {
        return children;
    }

    private void findChildren(Map<Class, List<InstanceProperty>> classes) {
        classes.putIfAbsent(this.data.getWrappedClass(), this.data.getProperties());

        data.getChildren().forEach((propName, wrapper) -> {
            if (wrapper != null) {
                classes.putIfAbsent(wrapper.getWrappedClass(), wrapper.getProperties());
            }

            children.add(new SimulationTree(wrapper, propName, classes));
        });
    }
}
