package etomica.dimer;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.box.Box;
import etomica.chem.elements.Copper;
import etomica.chem.elements.Tin;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisBetaSnA5;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.meam.ParameterSetMEAM;
import etomica.meam.PotentialMEAM;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.IVector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Simulation using Henkelman's Dimer method to find a saddle point for
 * an adatom of Sn on a surface, modeled with MEAM.
 * 
 * @author msellers
 *
 */

public class EnergyMapMEAMadatomSnCluster extends Simulation{

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "MEAM Md3D";
    public IntegratorEnergyMap integratorMAP;
    public Box box;
    public SpeciesSpheresMono sn, snFix, snAdatom, movable;
    public PotentialMEAM potential;
    
    public static void main(String[] args){
        double height;
        height = Double.parseDouble(args[0]);
        String fileTail = ""+height;
        
    	final String APP_NAME = "EnergyMapMEAMadatomSn";
    	final EnergyMapMEAMadatomSn sim = new EnergyMapMEAMadatomSn(height, fileTail);
    	
    	
    	sim.activityIntegrateMAP.setMaxSteps(1);
        sim.getController().actionPerformed();

    	
    }
}