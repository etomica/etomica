package etomica.data.meter;

import etomica.Atom;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Phase;
import etomica.Species;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorSpeciesDependent;
import etomica.space2d.Vector2D;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;

/**
 * Meter for measurement of density or species mole fraction within a specified subvolume
 * By selecting appropriate options, measurements can be taken for mole fraction of one species,
 * number density of one species, or total number density in the specified subvolume.
 * Subclasses define the shape to be rectangular, spherical, etc.
 */
 
 //setBounds is no longer implemented.  Class may need repair

public abstract class MeterLocalDensity extends MeterScalar implements EtomicaElement
{
    /**
     * Class variable used to specify that all species are included in number-density calculation
     *
     * @see #setSpeciesIndex
     */
    private Species species;
    private AtomIteratorSpeciesDependent iterator = new AtomIteratorLeafAtoms();
    protected double volume;
    /**
     * If <code>true</code> average mole fraction (n<sub>i</sub>/n<sub>total</sub>), if <code>false</code>, average total number density (n<sub>i</sub>/V)
     * Default is <code>false</code>
     */
    private boolean moleFractionMode = false;
    private String initialLabel;
    
    public MeterLocalDensity() {
        super();
        initialLabel = "Local Density";
        setLabel(initialLabel);
        setSpecies(null);
        computeVolume();      
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Local number density in a subregion of a phase");
        return info;
    }
    public Dimension getDimension() {return new DimensionRatio(Dimension.QUANTITY, Dimension.VOLUME);}
      
    /**
     * Accessor method to set the mole-fraction mode flag
     * Changes label appropriately
     *
     * @see #moleFractionMode
     */
    public final void setMoleFractionMode(boolean d) {
        moleFractionMode = d;
        if(moleFractionMode && getLabel().equals(initialLabel)) setLabel("Mole fraction");  //don't change label if already set to something else
        if(moleFractionMode && species == null)  {}
    }
    /**
     * Accessor method to get the value of the mole-fraction mode flag
     *
     * @see #moleFractionMode
     */
    public final boolean getMoleFractionMode() {return moleFractionMode;}
    
    /**
     * Accessor method to set which species mole-fraction or molar-density is averaged
     * To set to total number density, invoke with static ALL_SPECIES field as argument
     */
    public final void setSpecies(Species s) {species = s;}
    /**
     * Accessor method to get the current value of the species index
     *
     * @see #setSpeciesIndex
     */
    public final Species getSpecies() {return species;}
    
    /**
     * Sets the size and location of the meter, which may affect the definition of the local volume
     */
//    public final void setBounds(int x, int y, int width, int height) {
//        super.setBounds(x, y, width, height);
//        computeVolume();
//    }
    
    /**
     * Method to compute the volume of the local region where the density is measured, and set associated parameters
     */
    public abstract void computeVolume();
    
    /**
     * Method that specifies if a molecule is inside the local region where the density is measured
     */
    public abstract boolean contains(Atom m);

    /**
     * @return the current value of the local density or local mole fraction
     */
    public double getDataAsScalar(Phase p)
    {
        if(moleFractionMode) {  //compute local mole fraction
            int totalSum = 0, speciesSum = 0;
            iterator.setPhase(p);
            iterator.reset();
            while(iterator.hasNext()) {
                Atom m = iterator.nextAtom();
                if(this.contains(m)) {
                    totalSum++;
                    if(m.node.parentSpecies() == species) speciesSum++;
                }
            }
            if(totalSum == 0) return Double.NaN;
            return (double)speciesSum/(double)totalSum;
        }
        //compute local molar density
        int nSum = 0;
        if(species == null) {   //total density
            iterator.setPhase(p);
            iterator.reset();
            while(iterator.hasNext()) {
                if(this.contains(iterator.nextAtom())) nSum++;
            }
        }
        else {                              //species density
            iterator.setPhase(p);
            iterator.setSpecies(new Species[] {species});
            while(iterator.hasNext()) {
                if(this.contains(iterator.nextAtom())) nSum++;
           }
        }       
        return nSum/volume;
    }
    
        /**
        * Meter for measurement of the local-molecule-number density inside of a rectangular volume
        * Assumes a 2D space
        */

        public static class Cube extends etomica.data.meter.MeterLocalDensity
        {
            double xCenter, yCenter, halfWidth, halfHeight;
            
            public Cube()
            {
                super();
//                this.setBounds(0,0,300,300);
            }

            /**
             * Method to compute the volume of the local region where the density is measured
             * Volume is determined converting the size of the component (measured in pixels) to simulation units
             */
            public void computeVolume() {
/* NEEDS REPAIR                Rectangle rectangle = getBounds();
                halfWidth = 0.5*rectangle.width / BaseUnit.Length.Sim.TO_PIXELS;
                halfHeight = 0.5*rectangle.height / BaseUnit.Length.Sim.TO_PIXELS;
                volume = 4 * halfWidth * halfHeight;
                xCenter = rectangle.x/BaseUnit.Length.Sim.TO_PIXELS + halfWidth;
                yCenter = rectangle.y/BaseUnit.Length.Sim.TO_PIXELS + halfHeight;
*/            }
            
            /**
             * Method that specifies if a molecule center-of-mass is inside the local region where the density is measured
             */
            public boolean contains(Atom m) {
                Vector2D r = (Vector2D)m.coord.position();  //molecule center-of-mass
                if(Math.abs(r.x(0)-xCenter) > halfWidth) {return false;}
                else if(Math.abs(r.x(1)-yCenter) > halfHeight) {return false;}
                else {return true;}
            }
            
            /**
             * Method to render the local volume on the screen
             */
 /*           public void draw(Graphics g, int[] origin, double scale) {
/* NEEDS REPAIR                Rectangle rectangle = getBounds();
                int x0 = origin[0]+(int)(rectangle.x*scale);
                int y0 = origin[1]+(int)(rectangle.y*scale);
                int w = (int)(rectangle.width*scale);
                int h = (int)(rectangle.height*scale);
                g.drawRect(x0,y0,w,h);
           }  */
        }//end of Cube
            
}
