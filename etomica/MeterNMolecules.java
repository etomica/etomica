package simulate;

import simulate.*;
import java.beans.*;

public class MeterNMolecules extends simulate.Meter
{
    public MeterNMolecules()
    {
        super();
        setLabel("Molecules");
    }

    public double currentValue()
    {
        return phase.moleculeCount;
    }

}
