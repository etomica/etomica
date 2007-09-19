package etomica.virial;

import etomica.box.Box;
import etomica.potential.IPotential;
import etomica.space.Space;

/**
 * Mayer function that wraps another MayerFunctionSpherical and ignores 
 * the temperature passed to f(), asking the wrapped MayerFunction
 * for its value at the MayerFunctionSpherical Temperature's temperature.
 */
public class MayerFunctionSphericalTemperature extends MayerFunctionSpherical {

    public MayerFunctionSphericalTemperature(Space space, MayerFunctionSpherical fWrapped, double temperature) {
        super(space);
        mayerFunction = fWrapped;
        beta = 1/temperature;
    }

    public double f(double r2, double b) {
        return mayerFunction.f(r2,beta);
    }

    public IPotential getPotential() {
        return potential;
    }
    
    public void setBox(Box newBox) {
        potential.setBox(newBox);
    }

    private final MayerFunctionSpherical mayerFunction;
    private final double beta;
	private IPotential potential;

}
