package etomica.server;

import etomica.meta.InstanceProperty;
import etomica.meta.wrappers.SimulationWrapper;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SimulationModel {

    private final Map<Class, List<InstanceProperty>> classes;
    private final SimulationTree tree;

    public SimulationModel(SimulationWrapper base) {
        this.classes = new HashMap<>();
        this.tree = new PropertyTree(base, base.getWrappedClass().getSimpleName(), classes);
    }

    public Map<Class, List<InstanceProperty>> getClasses() {
        return classes;
    }

    public SimulationTree getTree() {
        return tree;
    }

}
