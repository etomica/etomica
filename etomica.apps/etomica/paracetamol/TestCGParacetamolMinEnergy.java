package etomica.paracetamol;

import etomica.action.PDBWriter;
import etomica.atom.IAtomPositioned;
import etomica.atom.IMolecule;
import etomica.space.Space;
import etomica.units.ElectronVolt;
import etomica.units.Kelvin;
import etomica.util.numerical.ConjugateGradientMultiDimensional;
import etomica.util.numerical.FiniteDifferenceDerivative;

/**
 * To minimize the configuration free energy by invoking non-linear conjugate gradient.
 */
public class TestCGParacetamolMinEnergy{

    public TestCGParacetamolMinEnergy() {
    	
    }

    public static void main(String[] args) {

        int nA = 192;
        double temperature = Kelvin.UNIT.toSim(100);
  
        // parse arguments
        if (args.length > 0) {
            temperature = Kelvin.UNIT.toSim(Double.parseDouble(args[0]));
        }

        System.out.println(ElectronVolt.UNIT.toSim(1));
        
        System.out.println("Running CG Energy Minimization Orthorhombic Paracetamol simulation");
        System.out.println(nA + " atoms " +" and temperature "+ Kelvin.UNIT.fromSim(temperature) +"K");

        // construct simulation
        SimulationOrthorhombicParacetamol sim = new SimulationOrthorhombicParacetamol(Space.getInstance(3), nA, temperature);
       
        // Energy Minimization
        NumericalDerivativeEnergyParacetamol function = new NumericalDerivativeEnergyParacetamol(sim.box, sim.potentialMaster);
        ConjugateGradientMultiDimensional conjugateGradient = new ConjugateGradientMultiDimensional();
       
        function.setCoordinateDefinition(sim.coordinateDefinition);
        
       // function.species = sim.getSpeciesManager().getSpecies()[0];
        
        
        System.out.println("CG Coordinate of first atom that overlaps: "
        		+((IAtomPositioned)((IMolecule)sim.box.getMoleculeList(sim.getSpeciesManager().getSpecies()[0]).getAtom(5))
        				.getChildList().getAtom(7)).getPosition());
        System.out.println("CG Coordinate of second atom that overlaps: "
        		+((IAtomPositioned)((IMolecule)sim.box.getMoleculeList(sim.getSpeciesManager().getSpecies()[0]).getAtom(7))
        				.getChildList().getAtom(18)).getPosition());
        
        
        double[] u = sim.coordinateDefinition.calcU(sim.coordinateDefinition.getBasisCells()[0].molecules);
        
 
        double[] num = new double[48];
        for (int i=0; i<48;i++){
        	if (i%6==3){
        		num[i]=-1e-17;
        	} else {
        	num[i] =1e-9;
        	}//System.out.print(u[i]+" ");
        }
        
        System.out.println("CG Coordinate of first atom after setToU 0: "
        		+((IAtomPositioned)((IMolecule)sim.box.getMoleculeList(sim.getSpeciesManager().getSpecies()[0]).getAtom(5))
        				.getChildList().getAtom(7)).getPosition());
//        System.out.println("CG Coordinate of second atom after setToU 0: "
//        		+((IAtomPositioned)((IAtomGroup)sim.box.getMoleculeList(sim.getSpeciesManager().getSpecies()[0]).getAtom(7))
//        				.getChildList().getAtom(18)).getPosition());
        
        sim.coordinateDefinition.setToU(sim.coordinateDefinition.getBasisCells()[0].molecules, num);
        System.out.println("CG Coordinate of first atom after setToU num: "
        		+((IAtomPositioned)((IMolecule)sim.box.getMoleculeList(sim.getSpeciesManager().getSpecies()[0]).getAtom(5))
       				.getChildList().getAtom(7)).getPosition());
        
        sim.coordinateDefinition.setToU(sim.coordinateDefinition.getBasisCells()[0].molecules, u);
        System.out.println("CG Coordinate of first atom after setToU u: "
        		+((IAtomPositioned)((IMolecule)sim.box.getMoleculeList(sim.getSpeciesManager().getSpecies()[0]).getAtom(5))
        				.getChildList().getAtom(7)).getPosition());
//        System.out.println("CG Coordinate of second atom after setToU u: "
//        		+((IAtomPositioned)((IAtomGroup)sim.box.getMoleculeList(sim.getSpeciesManager().getSpecies()[0]).getAtom(7))
//        				.getChildList().getAtom(18)).getPosition());
        
        
        
        conjugateGradient.conjugateGradient(u, 0.001, function);
        System.out.println("The minimum energy is: "+ conjugateGradient.getFunctionMinimimumValue());
        
        PDBWriter pdbWriter = new PDBWriter(sim.box);
        pdbWriter.setFileName("CG_Energy_Min_Ortho_"+ Kelvin.UNIT.fromSim(temperature)+"_K.pdb");
        pdbWriter.actionPerformed();
        
    }

    private static final long serialVersionUID = 1L;
    public CoordinateDefinitionParacetamol coordinateDefinition;
    public ConjugateGradientMultiDimensional nonLinearCG;
    protected FiniteDifferenceDerivative finiteDifferenceDerivative;
}