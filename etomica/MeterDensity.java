package simulate;

import simulate.*;
import java.beans.*;

public class MeterDensity extends simulate.Meter
{
    private final double scaleSquared = Constants.SCALE*Constants.SCALE;
    public MeterDensity()
    {
        super();
        setLabel("Density (mol/l)");
    }

    public double currentValue()
    {
        return phase.nMoleculeTotal/(phase.space.volume*scaleSquared*Constants.DEPTH*Constants.MOL_PER_LITER2SIM);
    }

}
