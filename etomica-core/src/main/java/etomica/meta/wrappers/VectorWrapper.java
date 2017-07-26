package etomica.meta.wrappers;

import etomica.meta.SimulationModel;
import etomica.meta.properties.VectorProperty;
import etomica.space.Vector;

/**
 * Created by kofke on 7/24/17.
 */
public class VectorWrapper extends Wrapper<Vector> {

    public VectorWrapper(Vector vector, SimulationModel simModel, boolean doSerialize) {
        super(vector, simModel, doSerialize);

        properties.add(new VectorProperty(vector));
 //       double[] x = new double[vector.getD()];
 //       properties.add(vector.assignTo(x));

        valueProps.addAll(properties);
    }
}
