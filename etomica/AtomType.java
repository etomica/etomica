package etomica;
import etomica.nbr.NeighborManagerAgent;
import etomica.units.Dimension;
//import etomica.electrostatics.*;

/**
 * AtomType holds fields that are common to many atoms. It serves many
 * functions:
 * <ul>
 * <li>it holds parameters that define the general features of the atom (e.g.,
 * size).
 * <li>it holds an index manager that is used to interpret the atoms index held
 * in its node
 * <li>it holds a class used to manage neighbors, by indicating which
 * potentials apply to the atom
 * </ul>
 * The AtomType of an atom is set by its AtomFactory when it builds a molecule.
 * Each Atom has an instance variable named "type" that holds the AtomType
 * instance.
 */

public class AtomType implements java.io.Serializable {

    public static Parameter.Source[] parameterSource = new Parameter.Source[0];
    AtomFactory creator;//set in constructor of AtomFactory
    public Parameter[] parameter;
    private Parameter.Size sizeParameter = Default.SIZE_PARAMETER;
    private Parameter.Energy energyParameter = Default.ENERGY_PARAMETER;
    private Parameter.Mass massParameter = Default.MASS_PARAMETER;
    protected int speciesIndex = -1;
    private Species species;
    
    //fields for linked list of all instances of AtomType
    public final AtomType previousInstance;
    private static AtomType lastInstance;
    
    private final NeighborManagerAgent neighborManagerAgent;
    
    public double mass, rm;
    private final AtomIndexManager indexManager;
    
//    private Parameter.Electrostatic electroParameter;
    
    public AtomType(AtomIndexManager indexManager) {
        this(indexManager, Default.ATOM_MASS);
    }
    public AtomType(AtomIndexManager indexManager, double mass) {
        
        this.indexManager = indexManager;
        
        //update linked list of instances
        this.previousInstance = lastInstance;
        lastInstance = this;
        
        //set up global parameters
        parameter = new Parameter[parameterSource.length];
        for(int i=0; i<parameter.length; i++) {
            parameter[i] = parameterSource[i].makeParameter();
        }
        setMass(mass);

//        System.out.println("AtomType constructor:"+mass);
        neighborManagerAgent = new NeighborManagerAgent();
    }
    
    public AtomIndexManager getIndexManager() {
        return indexManager;
    }
    
    /**
     * Returns true if an atom of this type is descended from
     * an atom (any atom) having the given type.
     */
    public boolean isDescendedFrom(AtomType type) {
        return indexManager.isDescendedFrom(type.indexManager);
    }
    
    protected void addGlobalParameter(Parameter.Source source) {
        Parameter[] newParameter = new Parameter[parameter.length+1];
        for(int i=0; i<parameter.length; i++) newParameter[i] = parameter[i];
        newParameter[parameter.length] = source.makeParameter();
        parameter = newParameter;
    }
    
    /**
     * Adds given parameter source to parameter-source array and returns index
     * indicating where in atomtype parameter-array the source's parameter will
     * be placed.
     */
    public static int requestParameterIndex(Parameter.Source source) {
        Parameter.Source[] newSource = new Parameter.Source[parameterSource.length+1];
        for(int i=0; i<parameterSource.length; i++) newSource[i] = parameterSource[i];
        int index = parameterSource.length;
        newSource[index] = source;
        parameterSource = newSource;
        
        //make parameter for any existing AtomType instances
        for(AtomType t=lastInstance; t!=null; t=t.previousInstance) {
            t.addGlobalParameter(source);
        }
        return index;
    }
    
    public AtomFactory creator() {return creator;}
    
    
    /**
     * Returns the depth of this atom in the atom hierarchy.  That is, returns
     * the number of parent relations between this atom and the species master.
     */
    public int getDepth() {
        return indexManager.getDepth();
    }
    /**
     * @return Returns the species.
     */
    public Species getSpecies() {
        return species;
    }
    /**
     * @param species The species to set.
     */
    public void setSpecies(Species species) {
        this.species = species;
        speciesIndex = species.getIndex();
    }
    /**
     * @return Returns the speciesIndex.
     */
    public final int getSpeciesIndex() {
        return speciesIndex;
    }

    /**
    * Sets  mass of this atom and updates reciprocal mass accordingly.  Setting
    * mass to largest machine double (Double.MAX_VALUE) causes reciprocal mass 
    * to be set to zero.
    * 
    * @param mass   new value for mass
    */
    public void setMass(double m) {
        massParameter.setMass(m);
        mass = m;
        rm = (m==Double.MAX_VALUE) ? 0.0 : 1.0/mass;
    }
    public final double getMass() {return mass;}
    public final Dimension getMassDimension() {return Dimension.MASS;}
    public final double rm() {return rm;}

    public NeighborManagerAgent getNbrManagerAgent() {
    	return neighborManagerAgent;
    }
    
    public boolean isInteracting() {
        return neighborManagerAgent.getPotentials().length > 0;
    }
    
    //prototype of a real atom type
/*    public final static class Carbon extends Sphere {
        public Carbon() {
            super(12.0, Color.black, 1.1);  //mass, color, diameter  
        }
    }
    public final static class Carbon12 extends Sphere {
        public Carbon12() {
            super(12.0, Color.black, 1.1);
            this.setName("Carbon" + Integer.toString(CarbonID++));
        }
    }
    
    public final static class Hydrogen extends Sphere {
        public Hydrogen() {
            super(1.0, Color.cyan, 0.5);
            this.setName("Hydrogen" + Integer.toString(HydrogenID++));
        }
    }
    
    public final static class Oxygen extends Sphere {
        public Oxygen() {
            super(16.0, Color.red, 1.3);
            this.setName("Oxygen" + Integer.toString(OxygenID++));
        }
    }
    
    public final static class Nitrogen extends Sphere {
        public Nitrogen() {
            super(14.0, Color.blue, 1.2);
            this.setName("Nitrogen" + Integer.toString(NitrogenID++));
        }
    }    
    
*/    
    //interfaces for anisotropic atom types
    public interface Rotator {
        public double[] momentOfInertia(); //diagonal elements of (diagonalized) moment of inertia; should always be a 3-element array
    }
    public interface SphericalTop extends Rotator {} //guarantees Ixx = Iyy = Izz
    public interface CylindricalTop extends Rotator {} //guarantees Ixx = Iyy
    public interface AsymmetricTop extends Rotator {} //all moment-of-inertia elements unequal
    
}
        
