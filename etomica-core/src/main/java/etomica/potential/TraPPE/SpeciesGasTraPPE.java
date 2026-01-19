package etomica.potential.TraPPE;

import etomica.atom.AtomType;
import etomica.chem.elements.*;
import etomica.space.Space;
import etomica.space3d.Vector3D;
import etomica.species.*;
import etomica.units.Degree;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.util.collections.IntArrayList;

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

    public enum ChemForm {
        CH4, C2H6, C3H8, C2H4, C3H6, CO2, N2, O2, NH3, C4H8ene_1, C4H8ene_13, C4H10, methylPropene_2, methylPropane_2, transButene, cisButene;
    }

    public AtomType[] atomTypes;
    public double[] sigma;
    public double[] epsilon;
    public double[] charge;
    public Map<Integer, Double> chargeMap;
    public ISpecies species;
    public Map<Integer, String> atomMap;

    public ISpecies speciesGasTraPPE(Space space, ChemForm chemForm, boolean isDynamic) {
        this.space = space;
        this.chemForm = chemForm;
        double cMass = Carbon.INSTANCE.getMass();
        double hMass = Hydrogen.INSTANCE.getMass();

        AtomType typeCH4 = AtomType.simple("ch4", cMass + 4 * hMass);
        AtomType typeCH3 = AtomType.simple("ch3", cMass + 3 * hMass);

        AtomType typeCH2ane = AtomType.simple("ch2ane", cMass + 3 * hMass);
        AtomType typeCH2ene = AtomType.simple("ch2ene", cMass + 3 * hMass);

        AtomType typeCHane = AtomType.simple("chane", cMass + hMass);
        AtomType typeCHene = AtomType.simple("chene", cMass + hMass);
        AtomType typeCHdiene = AtomType.simple("chdiene", cMass + hMass);

        AtomType typeC = AtomType.simple("c", cMass);

        if (chemForm == ChemForm.CH4) {
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
            Vector3D posnC3 = new Vector3D(new double[]{bondLength * Math.cos(angle), bondLength * Math.sin(angle), 0});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH3, posnC1, "C1")
                    .addAtom(typeCH2ane, posnC2, "C2")
                    .addAtom(typeCH3, posnC3, "C3")
                    .setDynamic(isDynamic)
                    .build();

        } else if (chemForm == ChemForm.C2H4) {

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

            Vector3D posnC1 = new Vector3D(new double[]{0, 0, 0});
            Vector3D posnC2 = new Vector3D(new double[]{1.54, 0, 0});
            Vector3D posnC3 = new Vector3D(new double[]{2.198, 1.15, 0});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH3, posnC1, "C1")
                    .addAtom(typeCHene, posnC2, "C2")
                    .addAtom(typeCH2ene, posnC3, "C3")
                    .setDynamic(isDynamic)
                    .build();

        } else if (chemForm == ChemForm.N2) {
            chargeMap = new HashMap<>();
            //Atoms in Compound
            AtomType typeM = new AtomType(elementM);
            AtomType typeN = new AtomType(Nitrogen.INSTANCE);

            atomTypes = new AtomType[]{typeM, typeN};

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
            chargeMap.put(0, qM);
            chargeMap.put(1, qN);
            chargeMap.put(2, qN);

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

            atomTypes = new AtomType[]{typeM, typeO};

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
            AtomType typeO = new AtomType(Oxygen.INSTANCE);

            atomTypes = new AtomType[]{typeC, typeO};

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

            atomTypes = new AtomType[]{typeN, typeH, typeM};

            polar = true;

            //TraPPE Parameters
            double bondLengthNH = 1.012; // Angstrom
            double bondLengthNM = 0.080; // Angstrom
            double thetaHNM = Degree.UNIT.toSim(67.9);
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
            sigma = new double[]{sigmaN, sigmaH, sigmaM};
            epsilon = new double[]{epsilonN, epsilonH, epsilonM};
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
        } else if (chemForm == ChemForm.C4H8ene_1) { //1 butene typeCH2aned = typeCH - typeCH2ane - typeCH3

            Vector3D pos1 = new Vector3D(new double[]{0.000000, 0.000000, 0.000000});
            Vector3D pos2 = new Vector3D(new double[]{1.330000, 0.000000, 0.000000});
            Vector3D pos3 = new Vector3D(new double[]{2.093006, 1.337693, 0.000000});
            Vector3D pos4 = new Vector3D(new double[]{3.625392, 1.184740, 0.000000});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH2ene, pos1, "ch2ene")
                    .addAtom(typeCHene, pos2, "chene")
                    .addAtom(typeCH2ane, pos3, "ch2ane")
                    .addAtom(typeCH3, pos4, "ch3")
                    .build();

        } else if (chemForm == ChemForm.C4H8ene_13) { //13 butadiene

            Vector3D pos1 = new Vector3D(new double[]{0.000000, 0.000000, 0.000000});
            Vector3D pos2 = new Vector3D(new double[]{1.330000, 0.000000, 0.000000});
            Vector3D pos3 = new Vector3D(new double[]{2.093006, 1.337693, 0.000000});
            Vector3D pos4 = new Vector3D(new double[]{3.423006, 1.337693, 0.000000});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH2ene, pos1, "ch2ene0")
                    .addAtom(typeCHdiene, pos2, "chdiene0")
                    .addAtom(typeCHdiene, pos3, "chdiene1")
                    .addAtom(typeCH2ene, pos4, "ch2ene1")
                    .build();

        } else if (chemForm == ChemForm.C4H10) { // butane

            Vector3D pos1 = new Vector3D(new double[]{0.000000, 0.000000, 0.000000});
            Vector3D pos2 = new Vector3D(new double[]{1.540000, 0.000000, 0.000000});
            Vector3D pos3 = new Vector3D(new double[]{2.166374, 1.406860, 0.000000});
            Vector3D pos4 = new Vector3D(new double[]{3.706374, 1.406860, 0.000000});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH3, pos1, "ch30")
                    .addAtom(typeCH2ane, pos2, "ch2ane0")
                    .addAtom(typeCH2ane, pos3, "ch2ane1")
                    .addAtom(typeCH3, pos4, "ch31")
                    .build();

        } else if (chemForm == ChemForm.methylPropane_2) {

            Vector3D pos1 = new Vector3D(new double[]{0.000, 0.000, 1.540});
            Vector3D pos2 = new Vector3D(new double[]{0.000, 0.000, 0.000});
            Vector3D pos3 = new Vector3D(new double[]{1.427, 0.000, -0.576});
            Vector3D pos4 = new Vector3D(new double[]{-0.855,1.143, -0.576});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH3, pos1, "ch30")
                    .addAtom(typeCHane, pos2, "chane1")
                    .addAtom(typeCH3, pos3, "ch31")
                    .addAtom(typeCH3, pos4, "ch32")
                    .build();

        } else if (chemForm == ChemForm.methylPropene_2) {

            Vector3D pos1 = new Vector3D(new double[]{0.000000, 0.000000, 0.000000});
            Vector3D pos2 = new Vector3D(new double[]{1.330000, 0.000000, 0.000000});
            Vector3D pos3 = new Vector3D(new double[]{-0.763, 1.337693, 0.000000});
            Vector3D pos4 = new Vector3D(new double[]{-0.763, -1.313, 0.252677});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH3, pos1, "ch30")
                    .addAtom(typeC, pos2, "c1")
                    .addAtom(typeCHene, pos3, "chene0")
                    .addAtom(typeCH3, pos4, "ch31")
                    .build();

        } else if (chemForm == ChemForm.cisButene) {

            Vector3D pos1 = new Vector3D(new double[]{-0.763006, 1.337693, 0.000000});
            Vector3D pos2 = new Vector3D(new double[]{0.000000, 0.000000, 0.000000});
            Vector3D pos3 = new Vector3D(new double[]{1.330000, 0.000000, 0.000000});
            Vector3D pos4 = new Vector3D(new double[]{2.093006, 1.337693, 0.000000});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH3, pos1, "ch30")
                    .addAtom(typeCHene, pos2, "chene0")
                    .addAtom(typeCHene, pos3, "chene1")
                    .addAtom(typeCH3, pos4, "ch31")
                    .build();

        } else if (chemForm == ChemForm.transButene) {

            Vector3D pos1 = new Vector3D(new double[]{-0.763006, 1.337693, 0.000000});
            Vector3D pos2 = new Vector3D(new double[]{0.000000, 0.000000, 0.000000});
            Vector3D pos3 = new Vector3D(new double[]{1.330000, 0.000000, 0.000000});
            Vector3D pos4 = new Vector3D(new double[]{2.093006, -1.337693, 0.000000});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH3, pos1, "ch30")
                    .addAtom(typeCHene, pos2, "chene0")
                    .addAtom(typeCHene, pos3, "chene1")
                    .addAtom(typeCH3, pos4, "ch31")
                    .build();

        } else {
            throw new RuntimeException("unrecognized chem form");
        }
        return species;
    }

    public List<double[]> getAngle(String molName) {
        List<double[]> doubleList = new ArrayList<>();
        double[] bondAngle = new double[2];
        if (molName.equals("butene")) {
            // typeCH2aned = typeCH - typeCH2ane - typeCH3
            bondAngle = new double[]{119.70, 70420};
            doubleList.add(bondAngle);
            bondAngle = new double[]{114, 62500};
            doubleList.add(bondAngle);
        } else if (molName.equals("propane")) {
            bondAngle = new double[]{114, 62500};
            doubleList.add(bondAngle);
        } else if (molName.equals("methylpropane")) {
            bondAngle = new double[]{112, 62500};
            doubleList.add(bondAngle);
        } else if (molName.equals("nbutane")) {
            bondAngle = new double[]{114, 62500};
            doubleList.add(bondAngle);
        } else if (molName.equals("propene")) {
            bondAngle = new double[]{119.70, 70420};
            doubleList.add(bondAngle);
        } else if (molName.equals("13butadiene")) {
            bondAngle = new double[]{119.70, 70420};
            doubleList.add(bondAngle);
        } else if (molName.equals("methylpropene")) {
            bondAngle = new double[]{119.70, 70420};
            doubleList.add(bondAngle);
        } else if (molName.equals("cisbutene")) {
            bondAngle = new double[]{119.70, 70420};
            doubleList.add(bondAngle);
        } else if (molName.equals("transbutene")) {
            bondAngle = new double[]{119.70, 70420};
            doubleList.add(bondAngle);
        }
        return doubleList;
    }

    public double[] getTorsion(String molName) {
        double[] bondTorsion = new double[4];
        if (molName.equals("butene")) { // c0 + c1 ( 1 + cos phi) + c2 ( 1 - cos 2phi)  + c3 ( 1 + cos 3phi)
            bondTorsion = new double[]{688.50, 86.36, -109.77, -282.24};
        } else if (molName.equals("C4H10")) {// c0 + c1 ( 1 + cos phi) + c2 ( 1 - cos 2phi)  + c3 ( 1 + cos 3phi)
            bondTorsion = new double[]{0.00, 355.03, -68.19, 791.32};
        } else if (molName.equals("13butadiene")) {// c0 + c1 cosphi + c2 cos2phi + c3 cos 3phi
            bondTorsion = new double[]{2034.58, 531.57, -1239.35, 460.04, 196.38};
        } else if (molName.equals("cisButene")) { // d0/2 (spi - spi0)^2
            bondTorsion = new double[]{180, 24800};
        } else if (molName.equals("transButene")) {// d0/2 (spi - spi0)^2
            bondTorsion = new double[]{0, 26800};
        }
        return bondTorsion;
    }

    public ISpecies speciesGasTraPPETwo(Space space, ChemForm chemForm, boolean isDynamic) {
        this.space = space;
        this.chemForm = chemForm;
        double cMass = Carbon.INSTANCE.getMass();
        double hMass = Hydrogen.INSTANCE.getMass();
        // AtomType ch2Type = AtomType.simple("ch2", cMass + 2 * hMass);
        // AtomType chType = AtomType.simple("ch", cMass +  hMass);
        //AtomType ch3Type = AtomType.simple("ch3", cMass + 3*hMass);
        if (chemForm == ChemForm.CH4) {

            AtomType typeCH4 = AtomType.simple("ch4", cMass + 4 * hMass);
            atomTypes = new AtomType[]{typeCH4};

            double sigmaCH4 = 3.73; // Angstrom
            double epsilonCH4 = Kelvin.UNIT.toSim(148);

            sigma = new double[]{sigmaCH4};
            epsilon = new double[]{epsilonCH4};
            charge = new double[]{0};
            //atomMap.put(0, typeCH4.getName());

            Vector3D posCH4 = new Vector3D(new double[]{0, 0, 0});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH4, posCH4, "CH4")
                    .setDynamic(isDynamic)
                    .build();

        } else if (chemForm == ChemForm.C2H6) {

            AtomType typeCH3n = AtomType.simple("ch3n", cMass + 3 * hMass);
            atomTypes = new AtomType[]{typeCH3n};

            double sigmaCH3n = 3.75;
            double epsilonCH3n = Kelvin.UNIT.toSim(98);
            double bondLength = 1.54;

            sigma = new double[]{sigmaCH3n};
            epsilon = new double[]{epsilonCH3n};
            charge = new double[]{0.0};
            //  atomMap.put(0, typeCH3.getName());
            // atomMap.put(1, typeCH3.getName());

            Vector3D posnC1 = new Vector3D(new double[]{0, 0, 0});
            Vector3D posnC2 = new Vector3D(new double[]{bondLength, 0, 0});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH3n, posnC1, "C1")
                    .addAtom(typeCH3n, posnC2, "C2")
                    .setDynamic(isDynamic)
                    .build();

        } else if (chemForm == ChemForm.C3H8) {

            AtomType typeCH3 = AtomType.simple("ch3ane", cMass + 3 * hMass);
            AtomType typeCH2 = AtomType.simple("ch2ane", cMass + 2 * hMass);
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
            // atomMap.put(0, typeCH3.getName());
            // atomMap.put(1, typeCH2.getName());
            // atomMap.put(2, typeCH2.getName());

            Vector3D posnC1 = new Vector3D(new double[]{bondLength, 0, 0});
            Vector3D posnC2 = new Vector3D(new double[]{0, 0, 0});
            Vector3D posnC3 = new Vector3D(new double[]{bondLength * Math.cos(angle), bondLength * Math.sin(angle), 0});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH3, posnC1, "C1")
                    .addAtom(typeCH2, posnC2, "C2")
                    .addAtom(typeCH3, posnC3, "C3")
                    .setDynamic(isDynamic)
                    .build();

        } else if (chemForm == ChemForm.C2H4) {

            AtomType typeCH2ene = AtomType.simple("ch2ene", cMass + 2 * hMass);
            atomTypes = new AtomType[]{typeCH2ene};

            double sigmaCH2n = 3.675;
            double epsilonCH2n = Kelvin.UNIT.toSim(85);
            double bondLength = 1.33;

            sigma = new double[]{sigmaCH2n};
            epsilon = new double[]{epsilonCH2n};
            charge = new double[]{0.0};
            //  atomMap.put(0, typeCH2ene.getName());
            //  atomMap.put(1, typeCH2ene.getName());

            Vector3D posnC1 = new Vector3D(new double[]{0, 0, 0});
            Vector3D posnC2 = new Vector3D(new double[]{bondLength, 0, 0});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH2ene, posnC1, "C1n")
                    .addAtom(typeCH2ene, posnC2, "C2n")
                    .setDynamic(isDynamic)
                    .build();

        } else if (chemForm == ChemForm.C3H6) {

            AtomType typeCH3 = AtomType.simple("ch3", cMass + 3 * hMass);
            AtomType typeCH2ene = AtomType.simple("ch2ene", cMass + 2 * hMass);
            AtomType typeCH1ene = AtomType.simple("ch1ene", cMass + hMass);
            atomTypes = new AtomType[]{typeCH3, typeCH1ene, typeCH2ene};

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

            Vector3D posnC1 = new Vector3D(new double[]{0, 0, 0});
            Vector3D posnC2 = new Vector3D(new double[]{1.54, 0, 0});
            Vector3D posnC3 = new Vector3D(new double[]{2.198, 1.15, 0});

            species = new SpeciesBuilder(space)
                    .addAtom(typeCH3, posnC1, "C1")
                    .addAtom(typeCH1ene, posnC2, "C2")
                    .addAtom(typeCH2ene, posnC3, "C3")
                    .setDynamic(isDynamic)
                    .build();

        } else if (chemForm == ChemForm.N2) {

            //Atoms in Compound
            AtomType typeM = new AtomType(elementM);
            AtomType typeN = new AtomType(Nitrogen.INSTANCE);

            atomTypes = new AtomType[]{typeM, typeN};

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
            chargeMap.put(0, qN);
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

            atomTypes = new AtomType[]{typeM, typeO};

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

            atomTypes = new AtomType[]{typeC, typeO};

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
            // atomMap.put(0, typeC.getName());
            // atomMap.put(1, typeO.getName());
            // atomMap.put(2, typeO.getName());

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

            atomTypes = new AtomType[]{typeN, typeH, typeM};

            polar = true;

            //TraPPE Parameters
            double bondLengthNH = 1.012; // Angstrom
            double bondLengthNM = 0.080; // Angstrom
            double thetaHNM = Degree.UNIT.toSim(67.9);
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
            sigma = new double[]{sigmaN, sigmaH, sigmaM};
            epsilon = new double[]{epsilonN, epsilonH, epsilonM};
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
        } else {
            throw new RuntimeException("unrecognized chem form");
        }
        return species;
    }

    public void setCharge(double[] charge) {
        this.charge = charge;
    }

    public double[] getCharge() {
        return charge;
    }


    public double atomicCharge(String gasName, String atomName) {
        double charge = 0.0;
        if (gasName.equals("n2")) {
            if (atomName.equals("N")) {
                charge = Electron.UNIT.toSim(-0.482);
            } else if (atomName.equals("M")) {
                charge = Electron.UNIT.toSim(0.964);
            }
        } else if (gasName.equals("co2")) {
            if (atomName.equals("C")) {
                charge = Electron.UNIT.toSim(0.700);
            } else if (atomName.equals("O")) {
                charge = Electron.UNIT.toSim(-0.350);
            } else if (atomName.equals("M")) {
                charge = 0.0;
            }
        }
        return charge;
    }

    public void setMap(Map<Integer, String> atomMap) {
        this.atomMap = atomMap;
    }
    public static IntArrayList[] getIntArrayList(String atomName){
        //check for larger alkanes-alkenes
        IntArrayList[] bonding;
        if (atomName.equals("ch4")){
            bonding = new IntArrayList[1];
            return bonding;
        } else if (atomName.equals("ethane") || atomName.equals("ethene")) {
            bonding = new IntArrayList[2];
            bonding[0] = new IntArrayList(new int[]{1});
            bonding[1] = new IntArrayList(new int[]{0});
            return bonding;
        } else if (atomName.equals("propane") || atomName.equals("propene")){
            bonding = new IntArrayList[3];
            bonding[0] = new IntArrayList(new int[]{1});
            bonding[1] = new IntArrayList(new int[]{0, 2});
            bonding[2] = new IntArrayList(new int[]{1});
            return bonding;
        } else if (atomName.equals("nbutane") || atomName.equals("butene") || atomName.equals("butadiene") ||
                atomName.equals("cisbutene") || atomName.equals("transbutene")) {
            bonding = new IntArrayList[4];
            bonding[0] = new IntArrayList(new int[]{1});
            bonding[1] = new IntArrayList(new int[]{0, 2});
            bonding[2] = new IntArrayList(new int[]{1, 3});
            bonding[3] = new IntArrayList(new int[]{2});
            return bonding;
        } else if (atomName.equals("methylpropane") || atomName.equals("methylpropene")) {
            bonding = new IntArrayList[4];
            bonding[0] = new IntArrayList(new int[]{1});
            bonding[1] = new IntArrayList(new int[]{0, 2, 3});
            bonding[2] = new IntArrayList(new int[]{1});
            bonding[3] = new IntArrayList(new int[]{1});
            return bonding;
        }else {
            throw new RuntimeException("Bonding for flexible angle not present");
        }
    }


    /*public double[] atomicPot(String atomtype) {
        HashMap<String, double[]> atomicConstant = new HashMap<>();
        // atomicConstant.put("ch4", new double[]{3.73, 148});
        atomicConstant.put("ch4n", new double[]{3.73, 148});
        atomicConstant.put("ch3", new double[]{3.75, 98});
        atomicConstant.put("ch3n", new double[]{3.75, 98});
        atomicConstant.put("ch2", new double[]{3.95, 46});
        atomicConstant.put("ch2n", new double[]{3.95, 46});
        atomicConstant.put("ch3ane", new double[]{3.75, 98});
        atomicConstant.put("ch2ane", new double[]{3.95, 46});
        atomicConstant.put("ch2ene", new double[]{3.675, 85});
        atomicConstant.put("ch1ene", new double[]{3.73, 47});
        atomicConstant.put("O", new double[]{3.050, 79});
        atomicConstant.put("C", new double[]{2.80, 27});
        atomicConstant.put("N", new double[]{3.310, 36});
        atomicConstant.put("M", new double[]{0, 0});
        atomicConstant.put("ch4", new double[]{3.72, 158.5});
        // Check if the atomtype exists in the map
        if (!atomicConstant.containsKey(atomtype)) {
            // If atomtype doesn't exist, throw a RuntimeException with a custom message
            throw new RuntimeException("Atom type '" + atomtype + "' does not exist in atomic constants.");
        }

        // If atomtype exists, return the corresponding values
        return atomicConstant.get(atomtype);
    }*/
    public double[] atomicPot(String atomtype) {
        HashMap<String, double[]> atomicConstant = new HashMap<>();
        //atomicConstant.put("ch4", new double[]{3.73, 148});
        atomicConstant.put("ch4", new double[]{3.73, 148});
       // atomicConstant.put("ch4", new double[]{3.72, 158.5});
        atomicConstant.put("ch3", new double[]{3.75, 98});
        atomicConstant.put("ch2ane", new double[]{3.95, 46});
        atomicConstant.put("chane", new double[]{4.68, 10});
        atomicConstant.put("ch2ene", new double[]{3.675, 85});
        atomicConstant.put("chene", new double[]{3.73, 47});
        atomicConstant.put("chdiene", new double[]{3.710, 52});
      //  atomicConstant.put("c", new double[]{3.75, 20});
        atomicConstant.put("M", new double[]{0.0, 0.0});
        atomicConstant.put("N", new double[]{3.310, 36});
        // Check if the atomtype exists in the map
        if (!atomicConstant.containsKey(atomtype)) {
            // If atomtype doesn't exist, throw a RuntimeException with a custom message
            throw new RuntimeException("Atom type '" + atomtype + "' does not exist in atomic constants.");
        }
        // If atomtype exists, return the corresponding values
        return atomicConstant.get(atomtype);
    }
}