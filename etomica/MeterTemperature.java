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
        // This method is derived from class simulate.Meter
        // to do: code goes here
        return phase.getKineticTemperature();
    }

}