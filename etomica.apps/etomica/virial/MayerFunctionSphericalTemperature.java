package etomica.virial;

import etomica.space.CoordinatePair;
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

    public double f(CoordinatePair cPair, double b) {
        return mayerFunction.f(cPair,beta);
    }

    private final MayerFunctionSpherical mayerFunction;
    private final double beta;

}
