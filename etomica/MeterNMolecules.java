package simulate;

import simulate.*;
import java.beans.*;

public class MeterNMolecules extends simulate.Meter
{
    public MeterNMolecules()
    {
        super();
        setLabel("Number of molecules");
    }

    public double currentValue()
    {
        return phase.moleculeCount;
    }

}
