package simulate;

public class MeterKineticEnergy extends simulate.Meter
{
    public MeterKineticEnergy()
    {
        super();
        setLabel("Kinetic Energy");
    }

    public double currentValue()
    {
        double ke = 0.0;
        for(Atom a=phase.firstAtom(); a!=null; a=a.nextAtom()) {
            ke += a.coordinate.kineticEnergy();
        }
        return ke;
    }
}