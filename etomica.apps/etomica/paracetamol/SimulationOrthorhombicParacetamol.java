package etomica.paracetamol;

import etomica.atom.AtomType;
import etomica.atom.AtomTypeGroup;
import etomica.box.Box;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.units.ElectronVolt;

/**
 * Simulation of 3D Paracetamol molecules in Form II (Orthorhombic)
 */
public class SimulationOrthorhombicParacetamol extends Simulation {

    public SimulationOrthorhombicParacetamol(Space space, int numAtoms, double temperature) {
        super(space, true);
        potentialMaster = new PotentialMaster(space);
        
        BasisOrthorhombicParacetamol basis = new BasisOrthorhombicParacetamol();;
        primitive = new PrimitiveOrthorhombic(space, 17.248, 12.086, 7.382);

        ConformationParacetamolOrthorhombic conformation = new ConformationParacetamolOrthorhombic(space);
        SpeciesParacetamol species = new SpeciesParacetamol(this);
        ((AtomTypeGroup)species.getMoleculeType()).setConformation(conformation);
        getSpeciesManager().addSpecies(species);

        box = new Box(this);
        addBox(box);
        box.setDimensions(Space.makeVector(new double[] {25,25,25}));
        box.setNMolecules(species, numAtoms);
        
            /*
             * Intermolecular Potential
             */
        
            double truncationRadiusCC   = 3.0* 3.395524116;
            double truncationRadiusCHy  = 3.0* 2.670105986;
            double truncationRadiusHyHy = 3.0* 2.099665865;
            double truncationRadiusCN   = 3.0* 3.237739512;
            double truncationRadiusNO   = 3.0* 3.035146951;
            double truncationRadiusNN   = 3.0* 3.087286897;
            double truncationRadiusHyN  = 3.0* 2.546030404;
            double truncationRadiusHyO  = 3.0* 2.503031484;
            double truncationRadiusOO   = 3.0* 2.983887553;
            double truncationRadiusCO   = 3.0* 3.183058614;
            double truncationRadiusHpHp = 3.0* 1.543178334;
            double truncationRadiusCHp  = 3.0* 2.289082817;
            double truncationRadiusHpN  = 3.0* 2.182712943;
            double truncationRadiusOHp  = 3.0* 2.145849932;
            double truncationRadiusHyHp = 3.0* 1.800044389;
            
            P2ElectrostaticDreiding potentialCC   = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(3832.14700),
            		0.277778, ElectronVolt.UNIT.toSim(25.286949));
            P2ElectrostaticDreiding potentialCHy  = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim( 689.53672),
            		0.272480, ElectronVolt.UNIT.toSim( 5.978972));
            P2ElectrostaticDreiding potentialHyHy = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim( 124.07167), 
            		0.267380, ElectronVolt.UNIT.toSim( 1.413698));
            P2ElectrostaticDreiding potentialCN   = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(3179.51460),
            		0.271003, ElectronVolt.UNIT.toSim(19.006710));
            P2ElectrostaticDreiding potentialNO   = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(2508.04480),
            		0.258398, ElectronVolt.UNIT.toSim(12.898341));
            P2ElectrostaticDreiding potentialNN   = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(2638.02850),
            		0.264550, ElectronVolt.UNIT.toSim(14.286224));
            P2ElectrostaticDreiding potentialHyN  = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim( 572.10541),
            		0.265957, ElectronVolt.UNIT.toSim( 4.494041));
            P2ElectrostaticDreiding potentialHyO  = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim( 543.91604),
            		0.259740, ElectronVolt.UNIT.toSim( 4.057452));
            P2ElectrostaticDreiding potentialOO   = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(2384.46580),
            		0.252525, ElectronVolt.UNIT.toSim(11.645288));
            P2ElectrostaticDreiding potentialCO   = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(3022.85020), 
            		0.264550, ElectronVolt.UNIT.toSim(17.160239));
            P2ElectrostaticDreiding potentialHpHp = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(  52.12899), 
            		0.214592, ElectronVolt.UNIT.toSim( 0.222819));
            P2ElectrostaticDreiding potentialCHp  = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim( 446.95185), 
            		0.242131, ElectronVolt.UNIT.toSim( 2.373693));
            P2ElectrostaticDreiding potentialHpN  = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim( 370.83387), 
            		0.236967, ElectronVolt.UNIT.toSim( 1.784166));
            P2ElectrostaticDreiding potentialOHp  = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim( 352.56176), 
            		0.232019, ElectronVolt.UNIT.toSim( 1.610837));
            P2ElectrostaticDreiding potentialHyHp = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(  80.42221), 
            		0.238095, ElectronVolt.UNIT.toSim( 0.561248));
            
            // CA-CA
            if(truncationRadiusCC > 0.5*box.getBoundary().getDimensions().x(0)) {
                throw new RuntimeException("Truncation radius too large.  " +
                		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
                }
            P2SoftSphericalTruncated interpotentialCC = new P2SoftSphericalTruncated (potentialCC, truncationRadiusCC); 
            potentialMaster.addPotential(interpotentialCC, new AtomType[]{(
            		(AtomFactoryParacetamol)species.getMoleculeFactory()).cType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).cType} );
            
            // CA-HY
            if(truncationRadiusCHy > 0.5*box.getBoundary().getDimensions().x(0)) {
                throw new RuntimeException("Truncation radius too large.  " +
                		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
                }
            P2SoftSphericalTruncated interpotentialCHy = new P2SoftSphericalTruncated (potentialCHy, truncationRadiusCHy); 
            potentialMaster.addPotential(interpotentialCHy, new AtomType[]{(
            		(AtomFactoryParacetamol)species.getMoleculeFactory()).cType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).hyType} );
            
            // HY-HY
            if(truncationRadiusHyHy > 0.5*box.getBoundary().getDimensions().x(0)) {
                throw new RuntimeException("Truncation radius too large.  " +
                		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
                }
            P2SoftSphericalTruncated interpotentialHyHy = new P2SoftSphericalTruncated (potentialHyHy, truncationRadiusHyHy); 
            potentialMaster.addPotential(interpotentialHyHy, new AtomType[]{(
            		(AtomFactoryParacetamol)species.getMoleculeFactory()).hyType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).hyType} );
                   
            // CA-NI
            if(truncationRadiusCN > 0.5*box.getBoundary().getDimensions().x(0)) {
                throw new RuntimeException("Truncation radius too large.  " +
                		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
                }
            P2SoftSphericalTruncated interpotentialCN = new P2SoftSphericalTruncated (potentialCN, truncationRadiusCN); 
            potentialMaster.addPotential(interpotentialCN, new AtomType[]{(
            		(AtomFactoryParacetamol)species.getMoleculeFactory()).cType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).nType} );
            
            // NI-OX
            if(truncationRadiusNO > 0.5*box.getBoundary().getDimensions().x(0)) {
                throw new RuntimeException("Truncation radius too large.  " +
                		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
                }
            P2SoftSphericalTruncated interpotentialNO = new P2SoftSphericalTruncated (potentialNO, truncationRadiusNO); 
            potentialMaster.addPotential(interpotentialNO, new AtomType[]{(
            		(AtomFactoryParacetamol)species.getMoleculeFactory()).nType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).oType} );
            
            //NI-NI
            if(truncationRadiusNN > 0.5*box.getBoundary().getDimensions().x(0)) {
                throw new RuntimeException("Truncation radius too large.  " +
                		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
                }
            P2SoftSphericalTruncated interpotentialNN = new P2SoftSphericalTruncated (potentialNN, truncationRadiusNN); 
            potentialMaster.addPotential(interpotentialNN, new AtomType[]{(
            		(AtomFactoryParacetamol)species.getMoleculeFactory()).nType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).nType} );
            
            // HY-NI
            if(truncationRadiusHyN > 0.5*box.getBoundary().getDimensions().x(0)) {
                throw new RuntimeException("Truncation radius too large.  " +
                		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
                }
            P2SoftSphericalTruncated interpotentialHyN = new P2SoftSphericalTruncated (potentialHyN, truncationRadiusHyN); 
            potentialMaster.addPotential(interpotentialHyN, new AtomType[]{(
            		(AtomFactoryParacetamol)species.getMoleculeFactory()).hyType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).nType} );
            
            // HY-OX
            if(truncationRadiusHyO > 0.5*box.getBoundary().getDimensions().x(0)) {
                throw new RuntimeException("Truncation radius too large. " +
                		" Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
                }
            P2SoftSphericalTruncated interpotentialHyO = new P2SoftSphericalTruncated (potentialHyO, truncationRadiusHyO); 
            potentialMaster.addPotential(interpotentialHyO, new AtomType[]{(
            		(AtomFactoryParacetamol)species.getMoleculeFactory()).hyType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).oType} );
                 
            // OX-OX
            if(truncationRadiusOO > 0.5*box.getBoundary().getDimensions().x(0)) {
                throw new RuntimeException("Truncation radius too large. " +
                		" Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
                }
            P2SoftSphericalTruncated interpotentialOO = new P2SoftSphericalTruncated (potentialOO, truncationRadiusOO); 
            potentialMaster.addPotential(interpotentialOO, new AtomType[]{(
            		(AtomFactoryParacetamol)species.getMoleculeFactory()).oType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).oType} );
            
            // CA-OX
            if(truncationRadiusCO > 0.5*box.getBoundary().getDimensions().x(0)) {
                throw new RuntimeException("Truncation radius too large. " +
                		" Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
                }
            P2SoftSphericalTruncated interpotentialCO = new P2SoftSphericalTruncated (potentialCO, truncationRadiusCO); 
            potentialMaster.addPotential(interpotentialCO, new AtomType[]{(
            		(AtomFactoryParacetamol)species.getMoleculeFactory()).cType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).oType} );
            
            // HP-HP
            if(truncationRadiusHpHp > 0.5*box.getBoundary().getDimensions().x(0)) {
                throw new RuntimeException("Truncation radius too large. " +
                		" Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
                }
            P2SoftSphericalTruncated interpotentialHpHp = new P2SoftSphericalTruncated (potentialHpHp, truncationRadiusHpHp); 
            potentialMaster.addPotential(interpotentialHpHp, new AtomType[]{(
            		(AtomFactoryParacetamol)species.getMoleculeFactory()).hpType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).hpType} );
            
            // CA-HP
            if(truncationRadiusCHp > 0.5*box.getBoundary().getDimensions().x(0)) {
                throw new RuntimeException("Truncation radius too large. " +
                		" Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
                }
            P2SoftSphericalTruncated interpotentialCHp = new P2SoftSphericalTruncated (potentialCHp, truncationRadiusCHp); 
            potentialMaster.addPotential(interpotentialCHp, new AtomType[]{(
            		(AtomFactoryParacetamol)species.getMoleculeFactory()).cType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).hpType} );
                   
            // HP-NI
            if(truncationRadiusHpN > 0.5*box.getBoundary().getDimensions().x(0)) {
                throw new RuntimeException("Truncation radius too large.  " +
                		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
                }
            P2SoftSphericalTruncated interpotentialHpN = new P2SoftSphericalTruncated (potentialHpN, truncationRadiusHpN); 
            potentialMaster.addPotential(interpotentialHpN, new AtomType[]{(
            		(AtomFactoryParacetamol)species.getMoleculeFactory()).hpType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).nType} );
            
            // OX-HP
            if(truncationRadiusOHp > 0.5*box.getBoundary().getDimensions().x(0)) {
                throw new RuntimeException("Truncation radius too large. " +
                		" Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
                }
            P2SoftSphericalTruncated interpotentialOHp = new P2SoftSphericalTruncated (potentialOHp, truncationRadiusOHp); 
            potentialMaster.addPotential(interpotentialOHp, new AtomType[]{(
            		(AtomFactoryParacetamol)species.getMoleculeFactory()).oType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).hpType} );
            
            // HY-HP
            if(truncationRadiusHyHp > 0.5*box.getBoundary().getDimensions().x(0)) {
                throw new RuntimeException("Truncation radius too large.  " +
                		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
                }
            P2SoftSphericalTruncated interpotentialHyHp = new P2SoftSphericalTruncated (potentialHyHp, truncationRadiusHyHp); 
            potentialMaster.addPotential(interpotentialHyHp, new AtomType[]{(
            		(AtomFactoryParacetamol)species.getMoleculeFactory()).hyType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).hpType} );
      
            potentialMaster.lrcMaster().setEnabled(false);
       
            
        boundary = new BoundaryRectangularPeriodic(space, getRandom(), 1);
        boundary.setDimensions(Space.makeVector(new double[] {2*17.248, 3*12.086, 4*7.382}));
        box.setBoundary(boundary);

        coordinateDefinition = new CoordinateDefinitionParacetamol(box, primitive, basis);
        coordinateDefinition.setBasisOrthorhombic();
        coordinateDefinition.initializeCoordinates(new int []{2,3,4});
        
    }
    
    private static final long serialVersionUID = 1L;
    public PotentialMaster potentialMaster;
    public Box box;
    public BoundaryRectangularPeriodic boundary;
    public PrimitiveOrthorhombic primitive;
    public CoordinateDefinitionParacetamol coordinateDefinition;
    
}