package simulate;

import simulate.*;
import java.beans.*;

public class MeterEnergy extends simulate.Meter
{
    public MeterEnergy()
    {
        super();
        setLabel("Energy");
    }

    public double currentValue()
    {
        // This method is derived from class simulate.Meter
        // to do: code goes here
        return phase.getTotalEnergy()/phase.nMoleculeTotal;
    }

}