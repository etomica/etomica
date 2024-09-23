package etomica.potential.TraPPE;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.*;
import etomica.graph.operations.Int;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.*;
import etomica.units.Degree;
import etomica.units.Electron;
import etomica.units.Kelvin;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SpeciesGasTraPPE {

    public static ChemForm ChemForm;
    protected boolean polar;
    protected static Element elementM = new ElementSimple("M", 0.0);
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
        CH4, C2H6, C3H8, C2H4, C3H6, CO2, N2, O2, NH3
    }
    public AtomType[] atomTypes;
    public double[] sigma;
    public double[] epsilon;
    public double[] charge;
    public Map<Integer, Double> chargeMap;
    public ISpecies species;
    public Map<Integer, String> atomMap;

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
            atomMap.put(0, typeCH4.getName());

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
            atomMap.put(0, typeCH3.getName());
            atomMap.put(1, typeCH3.getName());

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
            atomMap.put(0, typeCH3.getName());
            atomMap.put(1, typeCH2.getName());
            atomMap.put(2, typeCH2.getName());

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
            atomMap.put(0, typeCH2ene.getName());
            atomMap.put(1, typeCH2ene.getName());

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
            atomMap.put(0, typeCH3.getName());
            atomMap.put(1, typeCH2ene.getName());
            atomMap.put(2, typeCH1ene.getName());

            Vector3D posnC1 = new Vector3D(new double[]{bondLengthOne, 0, 0});
            Vector3D posnC2 = new Vector3D(new double[]{0, 0, 0});
            Vector3D posnC3 = new Vector3D(new double[]{bondLengthTwo*Math.cos(angle), bondLengthTwo*Math.sin(angle), 0});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH3, posnC1, "C1")
                    .addAtom(typeCH1ene, posnC2, "C2")
                    .addAtom(typeCH2ene, posnC3, "C3")
                    .setDynamic(isDynamic)
                    .build();

        }else if(chemForm == ChemForm.N2) {

            //Atoms in Compound
            AtomType typeM = new AtomType(elementM);
            AtomType typeN = new AtomType(Nitrogen.INSTANCE);

            atomTypes = new AtomType[]{typeM,typeN};

            //TraPPE Parameters
            double bondLength = 1.10; // Angstrom
            double sigmaN = 3.31; // Angstrom
            double epsilonN = Kelvin.UNIT.toSim(36.0);
            double qN = Electron.UNIT.toSim(-0.482);
            double sigmaM = 0.0; // Angstrom
            double epsilonM = Kelvin.UNIT.toSim(0.0);
            double qM = Electron.UNIT.toSim(0.964);

            //Construct Arrays
            sigma = new double[]{sigmaM, sigmaN};
            epsilon = new double[]{epsilonM, epsilonN};
            charge = new double[]{qM, qN};
            atomMap.put(0, typeN.getName());
            atomMap.put(1, typeN.getName());
            chargeMap.put(0,qN);
            chargeMap.put(1, qN);

            //Get Coordinates
            Vector3D posM = new Vector3D(new double[]{0, 0, 0});
            Vector3D posN1 = new Vector3D(new double[]{-bondLength / 2, 0, 0});
            Vector3D posN2 = new Vector3D(new double[]{+bondLength / 2, 0, 0});

            //Set Geometry
            species = new SpeciesBuilder(space)
                    .addAtom(typeM, posM, "M")
                    .addAtom(typeN, posN1, "N1")
                    .addAtom(typeN, posN2, "N2")
                    .build();
        } else if (chemForm == ChemForm.O2) {

            //Atoms in Compound
            AtomType typeM = new AtomType(elementM);
            AtomType typeO = new AtomType(Oxygen.INSTANCE);

            atomTypes = new AtomType[]{typeM,typeO};

            //TraPPE Parameters
            double bondLength = 1.210; // Angstrom
            double sigmaO = 3.020; // Angstrom
            double epsilonO = Kelvin.UNIT.toSim(49.0);
            double qO = Electron.UNIT.toSim(-0.113);
            double sigmaM = 0.0; // Angstrom
            double epsilonM = Kelvin.UNIT.toSim(0.0);
            double qM = Electron.UNIT.toSim(0.226);

            //Construct Arrays
            sigma = new double[]{sigmaM, sigmaO};
            epsilon = new double[]{epsilonM, epsilonO};
            charge = new double[]{qM, qO};
            atomMap.put(0, typeO.getName());
            atomMap.put(1, typeO.getName());
            chargeMap.put(0, qO);
            chargeMap.put(1, qO);

            //Get Coordinates
            Vector3D posM = new Vector3D(new double[]{0, 0, 0});
            Vector3D posO1 = new Vector3D(new double[]{-bondLength / 2, 0, 0});
            Vector3D posO2 = new Vector3D(new double[]{+bondLength / 2, 0, 0});

            //Set Geometry
            species = new SpeciesBuilder(space)
                    .addAtom(typeM, posM, "M")
                    .addAtom(typeO, posO1, "O1")
                    .addAtom(typeO, posO2, "O2")
                    .build();
        } else if (chemForm == ChemForm.CO2) {

            //Atoms in Compound
            AtomType typeC = new AtomType(Carbon.INSTANCE);
            AtomType typeO = new AtomType(Oxygen.INSTANCE);

            atomTypes = new AtomType[]{typeC,typeO};

            //TraPPE Parameters
            double bondLengthCO = 1.160; // Angstrom
            double sigmaC = 2.800; // Angstrom
            double epsilonC = Kelvin.UNIT.toSim(27.0);
            double qC = Electron.UNIT.toSim(0.700);
            double sigmaO = 3.050; // Angstrom
            double epsilonO = Kelvin.UNIT.toSim(79.0);
            double qO = Electron.UNIT.toSim(-0.350);

            //Construct Arrays
            sigma = new double[]{sigmaC, sigmaO};
            epsilon = new double[]{epsilonC, epsilonO};
            charge = new double[]{qC, qO};
            atomMap.put(0, typeC.getName());
            atomMap.put(1, typeO.getName());
            atomMap.put(2, typeO.getName());

            //Get Coordinates
            Vector3D posC = new Vector3D(new double[]{0, 0, 0});
            Vector3D posO1 = new Vector3D(new double[]{-bondLengthCO, 0, 0});
            Vector3D posO2 = new Vector3D(new double[]{+bondLengthCO, 0, 0});

            //Set Geometry
            species = new SpeciesBuilder(space)
                    .addAtom(typeC, posC, "C")
                    .addAtom(typeO, posO1, "O1")
                    .addAtom(typeO, posO2, "O2")
                    .build();
        } else if (chemForm == ChemForm.NH3) {

            //Atom in Compound
            AtomType typeN = new AtomType(Nitrogen.INSTANCE);
            AtomType typeH = new AtomType(Hydrogen.INSTANCE);
            AtomType typeM = new AtomType(elementM);

            atomTypes = new AtomType[]{typeN,typeH,typeM};

            polar = true;

            //TraPPE Parameters
            double bondLengthNH = 1.012; // Angstrom
            double bondLengthNM = 0.080; // Angstrom
            double thetaHNM = Degree.UNIT.toSim(67.9) ;
            double thetaHNH = Degree.UNIT.toSim(106.7);
            double thetaHNHxy = Degree.UNIT.toSim(60);
            double sigmaN = 3.420; // Angstrom
            double epsilonN = Kelvin.UNIT.toSim(185.0);
            double qN = Electron.UNIT.toSim(0.0);
            double sigmaH = 0.0; // Angstrom
            double epsilonH = Kelvin.UNIT.toSim(0.0);
            double qH = Electron.UNIT.toSim(0.410);
            double sigmaM = 0.0; // Angstrom
            double epsilonM = Kelvin.UNIT.toSim(0.0);
            double qM = Electron.UNIT.toSim(-1.230);

            //Construct Arrays
            sigma = new double[] {sigmaN,sigmaH,sigmaM};
            epsilon = new double[] {epsilonN,epsilonH,epsilonM};
            charge = new double[]{qN, qH, qM};

            //Get Coordinates
            Vector3D posN = new Vector3D(new double[]{0, 0, 0});
            Vector3D posH1 = new Vector3D(new double[]{bondLengthNH * Math.sin(thetaHNM), 0, -bondLengthNH * Math.cos(thetaHNM)});
            Vector3D posH2 = new Vector3D(new double[]{-bondLengthNH * Math.sin(thetaHNM) * Math.cos(thetaHNHxy), bondLengthNH * Math.sin(thetaHNM) * Math.sin(thetaHNHxy), -bondLengthNH * Math.cos(thetaHNM)});
            Vector3D posH3 = new Vector3D(new double[]{-bondLengthNH * Math.sin(thetaHNM) * Math.cos(thetaHNHxy), -bondLengthNH * Math.sin(thetaHNM) * Math.sin(thetaHNHxy), -bondLengthNH * Math.cos(thetaHNM)});
            Vector3D posM = new Vector3D(new double[]{0, 0, -bondLengthNM});

            //Set Geometry
            species = new SpeciesBuilder(space)
                    .addAtom(typeN, posN, "N")
                    .addAtom(typeH, posH1, "H1")
                    .addAtom(typeH, posH2, "H2")
                    .addAtom(typeH, posH3, "H3")
                    .addAtom(typeM, posM, "M")
                    .build();
        }else {
            throw new RuntimeException("unrecognized chem form");
        }
        return species;
    }

    public Map<Integer,int[]> getBonding(ChemForm chemForm){
        Map<Integer, int[]> getBonding = new HashMap<>();
        if(chemForm == ChemForm.C2H4){
            int[] bond = new int[]{0,1};
            getBonding.put(0,bond);
        } else if (chemForm == ChemForm.C2H6) {
            int[] bond = new int[]{0, 1};
            getBonding.put(0, bond);
        }
        return getBonding;
    }
    public void setCharge (double[] charge){
        this.charge = charge;
    }
    public double[] getCharge (){return charge;}

    public void setMap(Map<Integer, String > atomMap){this.atomMap = atomMap;}


    public double [] atomicPot (String atomtype){
        HashMap<String, double[]> atomicConstant = new HashMap<>();
        atomicConstant.put("ch4", new double[]{3.73, 148});
        atomicConstant.put("ch3", new double[]{3.75, 98});
        atomicConstant.put("ch2", new double[]{3.95, 46});
        atomicConstant.put("ch2ene", new double[]{3.675,85});
        atomicConstant.put("ch1ene", new double[]{3.73, 47});
        atomicConstant.put("O", new double[]{3.050, 79});
        atomicConstant.put("C", new double[]{2.80, 27});
        atomicConstant.put("N", new double[]{3.310, 36});
        atomicConstant.put("M", new double[]{0, 0});
        return atomicConstant.get(atomtype);
    }


}
