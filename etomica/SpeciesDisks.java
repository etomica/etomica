package etomica;
import java.awt.*;
import etomica.units.Dimension;

/**
 * Species in which molecules are made of arbitrary number of disks (same number for all molecules, though) 
 * with each disk having the same mass and size (same type).
 */
public class SpeciesDisks extends Species implements EtomicaElement {

//  The atomType is not declared final here becuase it makes setting up the constructors easier,
//  but effectively it cannot be changed once initialized; the instance is passed to the atoms, where
//  it is declared final
//  Note that the parameters of the type can be changed; only the instance of it is frozen once the atoms are made
//    (this is the same behavior as declaring it final)
    public AtomType.Disk protoType;
    private static AtomFactory makeFactory(Simulation sim, int na) {
        AtomFactoryMono f = new AtomFactoryMono(sim);
        AtomType type = new AtomType.Disk(f, Default.ATOM_MASS, Default.ATOM_COLOR, Default.ATOM_SIZE);
        f.setType(type);
        AtomFactoryMulti fm = new AtomFactoryMulti(sim,f);
        fm.setAtomsPerGroup(na);
  //      return fm;
        return f;
    }
        
    public SpeciesDisks() {
        this(Simulation.instance);
    }
    public SpeciesDisks(int n) {
        this(Simulation.instance, n);
    }
    public SpeciesDisks(Simulation sim) {
        this(sim, Default.MOLECULE_COUNT);
    }
    public SpeciesDisks(Simulation sim, int n) {
        this(sim, n, 1);
  /*      if(sim.space().D() == 1) {
            protoType = new AtomType.Rod(Default.ATOM_MASS,Default.ATOM_COLOR,Default.ATOM_SIZE);
        }
        else {
            protoType = new AtomType.Disk(Default.ATOM_MASS,Default.ATOM_COLOR,Default.ATOM_SIZE);
        }*/
        nMolecules = n;        
    }
    public SpeciesDisks(Simulation sim, int nM, int nA) {
        super(sim, makeFactory(sim, nA));
        nMolecules = nM;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Species with molecules composed of one or more spherical atoms");
        return info;
    }
              
    // Exposed Properties
    public final double getMass() {return protoType.mass();}
    public final void setMass(double mass) {protoType.setMass(mass);}
    public Dimension getMassDimension() {return Dimension.MASS;}
                
    public final double getDiameter() {return protoType.diameter();}
    public void setDiameter(double d) {protoType.setDiameter(d);}
    public Dimension getDiameterDimension() {return Dimension.LENGTH;}
                    
    public final Color getColor() {return protoType.color();}
    public final void setColor(Color c) {protoType.setColor(c);}
}


