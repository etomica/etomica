package simulate;

/**
 * Meter for measurement of density or species mole fraction within a specified volume
 * Subclasses define the shape to be rectangular, spherical, etc.
 */

public abstract class MeterLocalDensity extends simulate.Meter
{
    private final double scaleSquared = Constants.SCALE*Constants.SCALE;
    public final static int ALL_SPECIES = -100;
    int speciesIndex;
    double volume;
    boolean moleFractionMode = false;  //if true, compute mole fraction rather than molar density
    String initialLabel;
    
    public MeterLocalDensity()
    {
        super();
        initialLabel = "Local Density (mol/l)";
        setLabel(initialLabel);
        setSpeciesIndex(ALL_SPECIES);
        computeVolume();      
    }

    public final void setMoleFractionMode(boolean  d) {
        moleFractionMode = d;
        if(moleFractionMode && label.equals(initialLabel)) setLabel("Mole fraction");  //don't change label if already set to something else
        if(moleFractionMode && speciesIndex==ALL_SPECIES)  setSpeciesIndex(0);
    }
    public final boolean getMoleFractionMode() {return moleFractionMode;}
    
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
        if(moleFractionMode) {  //compute local mole fraction
            int totalSum = 0, speciesSum = 0;
            for(Molecule m=phase.firstMolecule(); m!=null; m=m.getNextMolecule()) {
                 if(this.contains(m) && !(m.firstAtom() instanceof AtomHardWall)) {
                    totalSum++;
                    if(m.getSpeciesIndex() == speciesIndex) speciesSum++;
                 }
            }
            if(totalSum == 0) return Double.NaN;
            return (double)speciesSum/(double)totalSum;
        }
        else {                 //compute local molar density
            int nSum = 0;
            if(speciesIndex == ALL_SPECIES) {   //total density
              for(Molecule m=phase.firstMolecule(); m!=null; m=m.getNextMolecule()) {
                 if(this.contains(m) && !(m.firstAtom() instanceof AtomHardWall)) nSum++;
              }}
            else {                              //species density
              for(Molecule m=phase.firstMolecule(); m!=null; m=m.getNextMolecule()) {
                 if(this.contains(m) && m.getSpeciesIndex()==speciesIndex) nSum++;
               }
            }       
            return nSum/(volume*scaleSquared*Constants.DEPTH*Constants.MOL_PER_LITER2SIM);
        }
    }
}
