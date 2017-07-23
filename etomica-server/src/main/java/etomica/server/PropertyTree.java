package etomica.server;

import etomica.meta.InstanceProperty;
import etomica.meta.wrappers.Wrapper;

import java.util.List;
import java.util.Map;

public class PropertyTree extends SimulationTree {

    private final String propName;

    public PropertyTree(Wrapper<?> wrapper, String propName, Map<Class, List<InstanceProperty>> classes) {
        super(wrapper, classes);
        this.propName = propName;
    }

    public String getPropName() {
        return propName;
    }
}
