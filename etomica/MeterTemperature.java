package simulate;

import simulate.*;
import java.beans.*;

public class MeterTemperature extends simulate.Meter
{
    public MeterTemperature()
    {
        super();
        setLabel("Temperature (K)");
    }

    public double currentValue()
    {
        double ke = phaseSpace.kineticEnergy.currentValue();
        return (2./(double)(phaseSpace.atomCount*Simulation.D))*ke*Constants.KE2T* Constants.SCALE * Constants.SCALE;
    }

}