package etomica;

import etomica.units.Dimension;

/**
 * Species in which molecules are made of a single atom of type OrientedSphere
 *
 * @author David Kofke
 * @see AtomType.OrientedSphere
 * 
 */
/* History
 * 08/12/03 (DAK) use sim instead of space in AtomFactoryMono constructor
 */
public class SpeciesSpheresRotating extends Species implements EtomicaElement {
    
    public double mass;
    
    public AtomType.OrientedSphere protoType;
    //static method used to make factory on-the-fly in the constructor
    private static AtomFactoryMono makeFactory(Simulation sim) {
        AtomFactoryMono f = new AtomFactoryMono(sim, sim.iteratorFactory.neighborSequencerFactory());
        AtomType type = new AtomType.OrientedSphere(f, Default.ATOM_MASS, Default.ATOM_SIZE);
        f.setType(type);
        return f;
    }
        
    public SpeciesSpheresRotating() {
        this(Simulation.instance);
    }
    public SpeciesSpheresRotating(int n) {
        this(Simulation.instance, n);
    }
    public SpeciesSpheresRotating(Simulation sim) {
        this(sim, Default.MOLECULE_COUNT);
    }
    public SpeciesSpheresRotating(Simulation sim, int nM) {
        super(sim, makeFactory(sim));
        factory.setSpecies(this);
        protoType = (AtomType.OrientedSphere)((AtomFactoryMono)factory).type();
        nMolecules = nM;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Molecules formed from spheres with an attached rotatable direction");
        return info;
    }
              
    // Exposed Properties
    public final double getMass() {return protoType.getMass();}
    public final void setMass(double m) {
        mass = m;
        allAtoms(new AtomAction() {public void actionPerformed(Atom a) {a.coord.setMass(mass);}});
    }
    public Dimension getMassDimension() {return Dimension.MASS;}
                
    public final double getDiameter() {return protoType.diameter(null);}
    public void setDiameter(double d) {protoType.setDiameter(d);}
    public Dimension getDiameterDimension() {return Dimension.LENGTH;}
    
/*    public static void main(String[] args) {
        P2RoughSphere.main(args);
    }
*/    
}


