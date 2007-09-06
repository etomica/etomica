package etomica.paracetamol;

import etomica.action.PDBWriter;
import etomica.space.Space;
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
  
        long simSteps = 10000;

        // parse arguments
        String filename = "Normal_Modes_Paracetamol_FormII_Min_"+ Kelvin.UNIT.fromSim(temperature)+"K";
        if (args.length > 0) {
            filename = args[0];
        }
    
        if (args.length > 1) {
            simSteps = Long.parseLong(args[1]);
        }
        if (args.length > 2) {
            temperature = Kelvin.UNIT.toSim(Double.parseDouble(args[2]));
        }

        System.out.println("Running "+ " Orthorhombic Paracetamol simulation");
        System.out.println(nA + " atoms " +" and temperature "+ Kelvin.UNIT.fromSim(temperature) +"K");
        System.out.println(simSteps+ " steps");
        System.out.println("output data to " + filename);

        // construct simulation
      
        
        SimulationOrthorhombicParacetamol sim = new SimulationOrthorhombicParacetamol(Space.getInstance(3), nA, temperature);
       
        // Energy Minimization
        EnergyFunctionParacetamol function = new EnergyFunctionParacetamol(sim.box, sim.potentialMaster);
        ConjugateGradientMultiDimensional conjugateGradient = new ConjugateGradientMultiDimensional();
       
        function.setCoordinateDefinition(sim.coordinateDefinition);
        double[] u = sim.coordinateDefinition.calcU(sim.coordinateDefinition.getBasisCells()[0].molecules);
        
        conjugateGradient.conjugateGradient(u, 0.001, function);
        System.out.println("The minimum energy is: "+ conjugateGradient.getFunctionMinimimumValue());
        
        PDBWriter pdbWriter = new PDBWriter(sim.box);
        pdbWriter.setFileName("Energy_Min_Ortho_"+ Kelvin.UNIT.fromSim(temperature)+"_K.pdb");
        pdbWriter.actionPerformed();
        
    }

    private static final long serialVersionUID = 1L;
    public CoordinateDefinitionParacetamol coordinateDefinition;
    public ConjugateGradientMultiDimensional nonLinearCG;
    protected FiniteDifferenceDerivative finiteDifferenceDerivative;
}