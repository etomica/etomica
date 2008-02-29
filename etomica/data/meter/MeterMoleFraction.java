package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.data.DataSourceScalar;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.species.ISpecies;
import etomica.units.Fraction;

/**
 * Meter for measurement of the species mole fraction in a box.
 * Mole fraction is defined (number of molecules of species)/(number of molecules in box).
 *
 * @author David Kofke
 */
public class MeterMoleFraction extends DataSourceScalar {
    private static final long serialVersionUID = 1L;
    private ISpecies species;
   
    public MeterMoleFraction() {
        super("Mole Fraction",Fraction.DIMENSION);
    }
    
    public void setSpecies(ISpecies s) {
        species = s;
    }
    public ISpecies getSpecies() {return species;}

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Species mole fraction in a box");
        return info;
    }

    public double getDataAsScalar() {
        if (box == null) throw new IllegalStateException("must call setBox before using meter");
    	return (species == null) ? Double.NaN :
         	(double)box.getNMolecules(species)/(double)box.moleculeCount();
     }

    /**
     * @return Returns the box.
     */
    public IBox getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
    public void setBox(Box box) {
        this.box = box;
    }

    private Box box;
}//end of MeterMoleFraction
