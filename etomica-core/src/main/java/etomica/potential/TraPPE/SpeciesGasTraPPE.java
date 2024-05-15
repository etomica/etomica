package etomica.potential.TraPPE;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.*;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.*;
import etomica.units.Degree;
import etomica.units.Kelvin;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class SpeciesGasTraPPE {

    public static SpeciesGasTraPPE.ChemForm ChemForm;
    protected ChemForm chemForm;
    protected Space space;
    public ISpecies getMethane(){
        double cMass = Carbon.INSTANCE.getMass();
        double hMass = Hydrogen.INSTANCE.getMass();
        species = SpeciesGeneral.monatomic(Space3D.getInstance(), AtomType.simple("CH4", cMass + 4 * hMass));
        return species;
    }

    public P2LennardJones speciesMethaneLJ(){
        double epsilon = Kelvin.UNIT.toSim(148);
        double sigma = 3.73;
        return new P2LennardJones(sigma, epsilon);
    }

    public ISpecies speciesAlkane(int nSpheres, PotentialMasterBonding pmBonding, Space space){
        ISpecies species = SpeciesAlkane.makeBuilder(nSpheres)
                .build();
        P3BondAngle p3 = new P3BondAngle(Math.PI*114.0/180.0, Kelvin.UNIT.toSim(62500));
        List<int[]> triplets = new ArrayList<>();
        for (int i=0; i<nSpheres-2; i++) {
            triplets.add(new int[]{i,i+1,i+2});
        }
        pmBonding.setBondingPotentialTriplet(species, p3, triplets);

        P4BondTorsion p4 = new P4BondTorsion(space, 0, Kelvin.UNIT.toSim(355.03), Kelvin.UNIT.toSim(-68.19), Kelvin.UNIT.toSim(791.32));
        List<int[]> quads = new ArrayList<>();
        for (int i=0; i<nSpheres-3; i++) {
            quads.add(new int[]{i,i+1,i+2,i+3});
        }
        pmBonding.setBondingPotentialQuad(species, p4, quads);
        return species;
    }

    public PotentialMasterCell setAlkaneLJ (ISpecies species, double rc, SpeciesManager sm, Box box, int cellRange, PotentialMasterBonding pmBonding, BondingInfo bondingInfo) {
        PotentialMasterCell potentialMaster = new PotentialMasterCell(sm, box, cellRange, pmBonding.getBondingInfo());
        AtomType typeCH3 = species.getAtomType(0);
        AtomType typeCH2 = species.getAtomType(1);
        double epsilonCH2 = Kelvin.UNIT.toSim(46.0);
        double epsilonCH3 = Kelvin.UNIT.toSim(98.0);
        double epsilonCH2CH3 = Math.sqrt(epsilonCH2*epsilonCH3);
        double sigmaCH2 = 3.95;
        double sigmaCH3 = 3.75;
        double sigmaCH2CH3 = (sigmaCH2+sigmaCH3)/2;
        TruncationFactory tf = new TruncationFactoryForceShift(rc);
        P2LennardJones p2CH2LJ = new P2LennardJones(sigmaCH2, epsilonCH2);
        P2LennardJones p2CH3LJ = new P2LennardJones(sigmaCH3, epsilonCH3);
        P2LennardJones p2CH2CH3LJ = new P2LennardJones(sigmaCH2CH3, epsilonCH2CH3);
        IPotential2 p2CH2 = tf.make(p2CH2LJ);
        IPotential2 p2CH3 = tf.make(p2CH3LJ);
        IPotential2 p2CH2CH3 = tf.make(p2CH2CH3LJ);
        potentialMaster.setPairPotential(typeCH2, typeCH2, p2CH2);
        potentialMaster.setPairPotential(typeCH2, typeCH3, p2CH2CH3);
        potentialMaster.setPairPotential(typeCH3, typeCH3, p2CH3);
        return potentialMaster;
    }

    public PotentialMasterCell setAlkeneLJ (ISpecies species, double rc, SpeciesManager sm, Box box, int cellRange, PotentialMasterBonding pmBonding, BondingInfo bondingInfo, String confName) {
        PotentialMasterCell potentialMaster = new PotentialMasterCell(sm, box, cellRange, pmBonding.getBondingInfo());
        if(confName.equals("ethene")){
            AtomType typeCH2 = species.getTypeByName("ch2");
            double epsilonCH2 = Kelvin.UNIT.toSim(85.0);
            double sigmaCH2 = 3.675;
            TruncationFactory tf = new TruncationFactoryForceShift(rc);
            P2LennardJones p2CH2LJ = new P2LennardJones(sigmaCH2, epsilonCH2);
            IPotential2 p2CH2 = tf.make(p2CH2LJ);
            potentialMaster.setPairPotential(typeCH2, typeCH2, p2CH2);
        } else if (confName.equals("propene")) {
            AtomType typeCH3 = species.getTypeByName("ch3");
            AtomType typeCH2 = species.getTypeByName("ch2");
            AtomType typeCH = species.getTypeByName("ch");
            double epsilonCH2 = Kelvin.UNIT.toSim(85.0);
            double epsilonCH3 = Kelvin.UNIT.toSim(98.0);
            double epsilonCH = Kelvin.UNIT.toSim(47.0);
            double epsilonCHCH3 = Math.sqrt(epsilonCH*epsilonCH3);
            double epsilonCHCH2 = Math.sqrt(epsilonCH*epsilonCH2);
            double epsilonCH2CH3 = Math.sqrt(epsilonCH2*epsilonCH3);
            double sigmaCH2 = 3.675;
            double sigmaCH3 = 3.75;
            double sigmaCH = 3.73;
            double sigmaCH2CH3 = (sigmaCH2+sigmaCH3)/2;
            double sigmaCHCH3 = (sigmaCH+sigmaCH3)/2;
            double sigmaCHCH2 = (sigmaCH+sigmaCH2)/2;
            P2LennardJones p2CH2LJ = new P2LennardJones(sigmaCH2, epsilonCH2);
            P2LennardJones p2CH3LJ = new P2LennardJones(sigmaCH3, epsilonCH3);
            P2LennardJones p2CHLJ = new P2LennardJones(sigmaCH, epsilonCH);
            P2LennardJones p2CH2CH3LJ = new P2LennardJones(sigmaCH2CH3, epsilonCH2CH3);
            P2LennardJones p2CHCH3LJ = new P2LennardJones(sigmaCHCH3, epsilonCHCH3);
            P2LennardJones p2CHCH2LJ = new P2LennardJones(sigmaCHCH2, epsilonCHCH2);
            TruncationFactory tf = new TruncationFactoryForceShift(rc);
            IPotential2 p2CH2 = tf.make(p2CH2LJ);
            IPotential2 p2CH3 = tf.make(p2CH3LJ);
            IPotential2 p2CH2CH3 = tf.make(p2CH2CH3LJ);
            IPotential2 p2CH = tf.make(p2CHLJ);
            IPotential2 p2CHCH3 = tf.make(p2CHCH3LJ);
            IPotential2 p2CHCH2 = tf.make(p2CHCH2LJ);
            potentialMaster.setPairPotential(typeCH2, typeCH2, p2CH2);
            potentialMaster.setPairPotential(typeCH, typeCH, p2CH);
            potentialMaster.setPairPotential(typeCH3, typeCH3, p2CH3);
            potentialMaster.setPairPotential(typeCH2, typeCH3, p2CH2CH3);
            potentialMaster.setPairPotential(typeCH, typeCH3, p2CHCH3);
            potentialMaster.setPairPotential(typeCH, typeCH2, p2CHCH2);
        }
        return potentialMaster;
    }


    public enum ChemForm{
        CH4, C2H6, C3H8, C2H4, C3H6
    }
    public AtomType[] atomTypes;
    public double[] sigma;
    public double[] epsilon;
    public double[] charge;
    public ISpecies species;

    public ISpecies speciesGasTraPPE(Space space, ChemForm chemForm, boolean isDynamic){
        this.space = space;
        this.chemForm = chemForm;
        double cMass = Carbon.INSTANCE.getMass();
        double hMass = Hydrogen.INSTANCE.getMass();
       // AtomType ch2Type = AtomType.simple("ch2", cMass + 2 * hMass);
       // AtomType chType = AtomType.simple("ch", cMass +  hMass);
        //AtomType ch3Type = AtomType.simple("ch3", cMass + 3*hMass);
        if(chemForm == ChemForm.CH4){

            AtomType typeCH4 = AtomType.simple("ch4", cMass + 4*hMass);
            atomTypes = new AtomType[]{typeCH4};

            double sigmaCH4 = 3.73; // Angstrom
            double epsilonCH4 = Kelvin.UNIT.toSim(148);

            sigma = new double[]{sigmaCH4};
            epsilon = new double[]{epsilonCH4};
            charge = new double[]{0};

            Vector3D posCH4 = new Vector3D(new double[]{0, 0, 0});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH4, posCH4, "CH4")
                    .setDynamic(isDynamic)
                    .build();


        } else if (chemForm == ChemForm.C2H6) {

            AtomType typeCH3 = AtomType.simple("ch3", cMass + 3*hMass);
            atomTypes = new AtomType[]{typeCH3};

            double sigmaCH3 = 3.75;
            double epsilonCH3 = Kelvin.UNIT.toSim(98);
            double bondLength = 1.54;

            sigma = new double[]{sigmaCH3};
            epsilon = new double[]{epsilonCH3};
            charge = new double[]{0.0};

            Vector3D posnC1 = new Vector3D(new double[]{0, 0, 0});
            Vector3D posnC2 = new Vector3D(new double[]{bondLength, 0, 0});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH3, posnC1, "C1")
                    .addAtom(typeCH3, posnC2, "C2")
                    .setDynamic(isDynamic)
                    .build();

        } else if (chemForm == ChemForm.C3H8) {

            AtomType typeCH3 = AtomType.simple("ch3", cMass + 3*hMass);
            AtomType typeCH2 = AtomType.simple("ch2", cMass + 2*hMass);
            atomTypes = new AtomType[]{typeCH3, typeCH2};

            double sigmaCH3 = 3.75;
            double sigmaCH2 = 3.95;
            double epsilonCH3 = Kelvin.UNIT.toSim(98);
            double epsilonCH2 = Kelvin.UNIT.toSim(46);
            double bondLength = 1.54;
            double angle = Degree.UNIT.toSim(114);

            sigma = new double[]{sigmaCH3, sigmaCH2};
            epsilon = new double[]{epsilonCH3, epsilonCH2};
            charge = new double[]{0.0};

            Vector3D posnC1 = new Vector3D(new double[]{bondLength, 0, 0});
            Vector3D posnC2 = new Vector3D(new double[]{0, 0, 0});
            Vector3D posnC3 = new Vector3D(new double[]{bondLength* Math.cos(angle), bondLength*Math.sin(angle), 0});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH3, posnC1, "C1")
                    .addAtom(typeCH2, posnC2, "C2")
                    .addAtom(typeCH3, posnC3, "C3")
                    .setDynamic(isDynamic)
                    .build();

        } else if (chemForm == ChemForm.C2H4) {

            AtomType typeCH2ene = AtomType.simple("ch2ene", cMass + 2*hMass);
            atomTypes = new AtomType[]{typeCH2ene};

            double sigmaCH2 = 3.675;
            double epsilonCH2 = Kelvin.UNIT.toSim(85);
            double bondLength = 1.33;

            sigma = new double[]{sigmaCH2};
            epsilon = new double[]{epsilonCH2};
            charge = new double[]{0.0};

            Vector3D posnC1 = new Vector3D(new double[]{0, 0, 0});
            Vector3D posnC2 = new Vector3D(new double[]{bondLength, 0, 0});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH2ene, posnC1, "C1")
                    .addAtom(typeCH2ene, posnC2, "C2")
                    .setDynamic(isDynamic)
                    .build();

        } else if (chemForm == ChemForm.C3H6) {

            AtomType typeCH3 = AtomType.simple("ch3", cMass + 3*hMass);
            AtomType typeCH2ene = AtomType.simple("ch2ene", cMass + 2*hMass);
            AtomType typeCH1ene = AtomType.simple("ch1ene", cMass + 2*hMass);
            atomTypes = new AtomType[]{typeCH3,typeCH1ene, typeCH2ene};

            double sigmaCH3 = 3.75;
            double epsilonCH3 = Kelvin.UNIT.toSim(98);
            double sigmaCH2 = 3.675;
            double epsilonCH2 = Kelvin.UNIT.toSim(85);
            double sigmaCH = 3.73;
            double epsilonCH = Kelvin.UNIT.toSim(47);
            double bondLengthOne = 1.54;
            double bondLengthTwo = 1.33;
            double angle = 119.70;

            sigma = new double[]{sigmaCH3, sigmaCH2, sigmaCH};
            epsilon = new double[]{epsilonCH3, epsilonCH2, epsilonCH};
            charge = new double[]{0.0};

            Vector3D posnC1 = new Vector3D(new double[]{bondLengthOne, 0, 0});
            Vector3D posnC2 = new Vector3D(new double[]{0, 0, 0});
            Vector3D posnC3 = new Vector3D(new double[]{bondLengthTwo*Math.cos(angle), bondLengthTwo*Math.sin(angle), 0});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH3, posnC1, "C1")
                    .addAtom(typeCH1ene, posnC2, "C2")
                    .addAtom(typeCH2ene, posnC3, "C3")
                    .setDynamic(isDynamic)
                    .build();

        }else {
            throw new RuntimeException("unrecognized chem form");
        }
        return species;
    }


    public double [] atomicPot (String atomtype){
        HashMap<String, double[]> atomicConstant = new HashMap<>();
        atomicConstant.put("ch4", new double[]{3.73, 148});
        atomicConstant.put("ch3", new double[]{3.75, 98});
        atomicConstant.put("ch2", new double[]{3.95, 46});
        atomicConstant.put("ch2ene", new double[]{3.675,85});
        atomicConstant.put("ch1ene", new double[]{3.73, 47});
        return atomicConstant.get(atomtype);
    }

}
