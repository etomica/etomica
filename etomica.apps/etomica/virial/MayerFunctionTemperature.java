package etomica.virial;

import etomica.AtomSet;

/**
 * Mayer function that wraps another MayerFunction and ignores 
 * the temperature passed to f(), asking the wrapped MayerFunction
 * for its value at the MayerFunctionTemperature's temperature.
 */
public class MayerFunctionTemperature implements MayerFunction, java.io.Serializable {

    public MayerFunctionTemperature(MayerFunction fWrapped, double temperature) {
        beta = 1/temperature;
        mayerFunction = fWrapped;
    }

    public double f(AtomSet pair, double b) {
        return mayerFunction.f(pair,beta);
    }

    private final MayerFunction mayerFunction;
    private final double beta;
}
