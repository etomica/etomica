package etomica;

import java.awt.Color;
import etomica.units.Dimension;

/**
 * Species in which molecules are made of a single atom of type OrientedSphere
 *
 * @author David Kofke
 * @see AtomType.OrientedSphere
 * 
 */
public class SpeciesSpheresRotating extends Species implements EtomicaElement {
    
    public double mass;
    
    public AtomType.OrientedSphere protoType;
    //static method used to make factory on-the-fly in the constructor
    private static AtomFactoryMono makeFactory(Simulation sim) {
        AtomFactoryMono f = new AtomFactoryMono(sim);
        AtomType type = new AtomType.OrientedSphere(f, Default.ATOM_MASS, Default.ATOM_COLOR, Default.ATOM_SIZE);
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
        protoType = (AtomType.OrientedSphere)((AtomFactoryMono)factory).type();
        nMolecules = nM;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Molecules formed from spheres with an attached rotatable direction");
        return info;
    }
              
    // Exposed Properties
    public final double getMass() {return protoType.mass();}
    public final void setMass(double m) {
        mass = m;
        allAtoms(new AtomAction() {public void actionPerformed(Atom a) {a.coord.setMass(mass);}});
    }
    public Dimension getMassDimension() {return Dimension.MASS;}
                
    public final double getDiameter() {return protoType.diameter();}
    public void setDiameter(double d) {protoType.setDiameter(d);}
    public Dimension getDiameterDimension() {return Dimension.LENGTH;}
                    
    public final Color getColor() {return protoType.color();}
    public final void setColor(Color c) {protoType.setColor(c);}
    
    public static void main(String[] args) {
        P2RoughSphere.main(args);
    }
    
}


