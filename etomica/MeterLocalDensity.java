package simulate;

public abstract class MeterLocalDensity extends simulate.Meter
{
    private final double scaleSquared = Constants.SCALE*Constants.SCALE;
    public final static int ALL_SPECIES = -100;
    int speciesIndex;
    double volume;
    
    public MeterLocalDensity()
    {
        super();
        setLabel("Local Density (mol/l)");
        setSpeciesIndex(ALL_SPECIES);
    }

    public final void setSpeciesIndex(int i) {speciesIndex = i;}
    public final int getSpeciesIndex() {return speciesIndex;}
    
    public final void setBounds(int x, int y, int width, int height) {
        super.setBounds(x, y, width, height);
        computeVolume();
    }
    
    /**
     Method to compute the volume of the local region where the density is measured, and set associated parameters
     */
    public abstract void computeVolume();
    
    /**
     Method that specifies if a molecule is inside the local region where the density is measured
     */
    public abstract boolean contains(Molecule m);

    public double currentValue()
    {
        int nSum = 0;
        for(Molecule m=phase.firstMolecule(); m!=null; m=m.getNextMolecule()) {
            if(this.contains(m)) nSum++;
        }
        return nSum/(volume*scaleSquared*Constants.DEPTH*Constants.MOL_PER_LITER2SIM);
    }

}
