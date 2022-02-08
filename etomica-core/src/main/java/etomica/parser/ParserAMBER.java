package etomica.parser;

import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.box.Box;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.config.ConformationFile;
import etomica.config.ConformationGeneric;
import etomica.config.IConformation;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.potential.compute.NeighborManagerSimple;
import etomica.potential.compute.PotentialComputePair;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesManager;
import etomica.units.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
/*
 * AMBER Parameter Topology file (.prmtop) parser
 * @author Arpit Bansal
 */

public class ParserAMBER {
    public static Stuff makeStuffFromLines(String smiles, List<String> lines, Options opts) {

        System.out.println("Smiles " + smiles);
        String heading = null;
        AtomType[] atomTypes = null;
        P2Harmonic[] p2Bonds = null;
        P2LennardJones[] p2LennardJones = null;
        P3BondAngle[] p3Bonds = null;
        P4BondTorsion[] p4Bonds = null;
        int[] atomCounts = null;
        int[] atomTypeId = null;
        double[] charges = null;
        Vector[] coords = null;
        List<int[]>[] bondedPairs = null;
        List<int[]>[] bondedTriplets = null;
        List<int[]>[] bondedQuads = null;
        char symbol = 'A';
        double[] sigma = null;
        double[] epsilon = null;
        double[] A = null;
        double[] B = null;
        int uniqueAtoms = 0;
        int uniqueBonds = 0;
        int uniqueAngles = 0;
        int uniqueTorsions = 0;
        double[] w = null;
        double[] r0 = null;
        double[] k_epsilon = null;
        double[] angle = null;
        double[] k = null;
        double[] n = null;
        double[] phi = null;
        int iter = 0;
        int countC = 0;
        int countH = 0;
        int numBondsIH = 0;
        int numBondsWH = 0;
        int numAnglesIH = 0;
        int numAnglesWH = 0;
        int numTorsionsIH = 0;
        int numTorsionsWH = 0;
        int[] atomicNumber = null;
        double[] mass = null;
        // let's assume we have "real" AMBER units -- Angstroms, kCal/mol
        Unit eUnit = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT), Mole.UNIT);
        SpeciesBuilder builder = new SpeciesBuilder(Space3D.getInstance());

        Iterator<String> linesIterator = lines.iterator();
        while (linesIterator.hasNext()) {
            String a = linesIterator.next();
            if (a.matches("%FLAG POINTERS")) {
                linesIterator.next();
                a = linesIterator.next();
                int i = 0;
                while (!a.contains("%")) {
                    String[] fields = a.split("\\s+");
                    if (i == 0) {
                        int numAtoms = Integer.parseInt(fields[0]);
                        uniqueAtoms = Integer.parseInt(fields[1]);
                        atomTypeId = new int[numAtoms];
                        atomicNumber = new int[numAtoms];
                        mass = new double[numAtoms];
                        atomTypes = new AtomType[uniqueAtoms];
                        atomCounts = new int[numAtoms];
                        coords = new Vector[numAtoms];
                        sigma = new double[uniqueAtoms + (uniqueAtoms * (uniqueAtoms - 1) / 2)];
                        epsilon = new double[uniqueAtoms + (uniqueAtoms * (uniqueAtoms - 1) / 2)];
                        A = new double[uniqueAtoms + (uniqueAtoms * (uniqueAtoms - 1) / 2)];
                        B = new double[uniqueAtoms + (uniqueAtoms * (uniqueAtoms - 1) / 2)];
                        p2LennardJones = new P2LennardJones[uniqueAtoms + (uniqueAtoms * (uniqueAtoms - 1) / 2)];
                        //System.out.println(String.format("Num Atoms %d", numAtoms));
                        numBondsIH = Integer.parseInt(fields[2]);
                        //System.out.println(String.format("Num Bonds including Hydrogen %d", numBondsIH));
                        numBondsWH = Integer.parseInt(fields[3]);
                        //System.out.println(String.format("Num Bonds without Hydrogen %d", numBondsWH));
                        numAnglesIH = Integer.parseInt(fields[4]);
                        //System.out.println(String.format("Num Angles including Hydrogen %d", numAnglesIH));
                        numAnglesWH = Integer.parseInt(fields[5]);
                        //System.out.println(String.format("Num Angles without Hydrogen %d", numAnglesWH));
                        numTorsionsIH = Integer.parseInt(fields[6]);
                        //System.out.println(String.format("Num Torsions including Hydrogen %d", numTorsionsIH));
                        numTorsionsWH = Integer.parseInt(fields[7]);
                        //System.out.println(String.format("Num Torsions without Hydrogen %d", numTorsionsWH));
                    }
                    if (i == 1) {
                        uniqueBonds = Integer.parseInt(fields[5]);
                        p2Bonds = new P2Harmonic[uniqueBonds];
                        bondedPairs = new ArrayList[p2Bonds.length];
                        for (int j = 0; j < p2Bonds.length; j++) {
                            bondedPairs[j] = new ArrayList<int[]>();
                        }
                        w = new double[uniqueBonds];
                        r0 = new double[uniqueBonds];
                        //System.out.println(String.format("Num Unique Bonds %d", uniqueBonds));
                        uniqueAngles = Integer.parseInt(fields[6]);
                        p3Bonds = new P3BondAngle[uniqueAngles];
                        bondedTriplets = new ArrayList[p3Bonds.length];
                        for (int j = 0; j < p3Bonds.length; j++) {
                            bondedTriplets[j] = new ArrayList<int[]>();
                        }
                        k_epsilon = new double[uniqueAngles];
                        angle = new double[uniqueAngles];
                        //System.out.println(String.format("Num Unique Angles %d", uniqueAngles));
                        uniqueTorsions = Integer.parseInt(fields[7]);
                        p4Bonds = new P4BondTorsion[uniqueTorsions];
                        bondedQuads = new ArrayList[p4Bonds.length];
                        for (int j = 0; j < p4Bonds.length; j++) {
                            bondedQuads[j] = new ArrayList<int[]>();
                        }
                        k = new double[uniqueTorsions];
                        n = new double[uniqueTorsions];
                        phi = new double[uniqueTorsions];
                        //System.out.println(String.format("Num Unique Torsions %d", uniqueTorsions));
                    }
                    ++i;
                    a = linesIterator.next();
                }
            }
            if (a.matches("%FLAG ATOMIC_NUMBER")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Atomic Number");
                iter = 0;
                while (!a.contains("%")) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; ++i) {
                        atomicNumber[10 * iter + i] = Integer.parseInt(fields[i]);
                        //System.out.println(atomicNumber[10 * iter + i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
            }
            if (a.matches("%FLAG MASS")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Mass");
                iter = 0;
                while (!a.contains("%")) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; i++) {
                        mass[5 * iter + i] = Double.parseDouble(fields[i]);
                        //System.out.println(mass[5 * iter + i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
            }
            if (a.matches("%FLAG ATOM_TYPE_INDEX")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Atom Type");
                iter = 0;
                while (!a.contains("%")) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; ++i) {
                        atomTypeId[10 * iter + i] = Integer.parseInt(fields[i]);
                        //atomTypes[10 * iter + i] = new AtomType(new ElementSimple(String.format("%d", atomTypeId[10 * iter + i])));
                        //System.out.println(atomTypeId[10 * iter + i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
            }
            if (a.matches("%FLAG BOND_FORCE_CONSTANT")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Bond Force Constant (in kcal/mol/A^2)");
                iter = 0;
                while (!a.contains("%") && a.length()>0) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; ++i) {
                        w[10 * iter + i] = Double.parseDouble(fields[i]);
                        //System.out.println(w[10 * iter + i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
            }
            if (a.matches("%FLAG BOND_EQUIL_VALUE")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Bond Equilibrium Distance (in A)");
                iter = 0;
                while (!a.contains("%") && a.length()>0) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; ++i) {
                        r0[10 * iter + i] = Double.parseDouble(fields[i]);
                        p2Bonds[10 * iter + i] = new P2Harmonic(opts.space, eUnit.toSim(w[10 * iter + i]), r0[10 * iter + i]);
                        //System.out.println(r0[10 * iter + i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
            }
            if (a.matches("%FLAG ANGLE_FORCE_CONSTANT")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Angle Force Constant (in kcal*rad^2/mol)");
                iter = 0;
                while (!a.contains("%") && a.length()>0) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; ++i) {
                        k_epsilon[10 * iter + i] = Double.parseDouble(fields[i]);
                        //System.out.println(k_epsilon[i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
            }
            if (a.matches("%FLAG ANGLE_EQUIL_VALUE")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Angle Equilibrium Angles (in radians)");
                iter = 0;
                while (!a.contains("%") && a.length()>0) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; ++i) {
                        angle[10 * iter + i] = Double.parseDouble(fields[i]);
                        p3Bonds[10 * iter + i] = new P3BondAngle(opts.space);
                        p3Bonds[10 * iter + i].setEpsilon(eUnit.toSim(k_epsilon[10 * iter + i]));
                        p3Bonds[10 * iter + i].setAngle(angle[10 * iter + i]);
                        //System.out.println(angle[10 * iter + i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
            }
            if (a.matches("%FLAG DIHEDRAL_FORCE_CONSTANT")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Dihedral Force Constant (in kcal/mol)");
                iter = 0;
                while (!a.contains("%") && a.length()>0) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    //System.out.println("|"+fields[0]+"|");
                    for (int i = 0; i < fields.length; ++i) {
                        k[5 * iter + i] = Double.parseDouble(fields[i]);
                        //System.out.println(k[5 * iter + i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
            }
            if (a.matches("%FLAG DIHEDRAL_PERIODICITY")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Dihedral Periodicity");
                iter = 0;
                while (!a.contains("%") && a.length()>0) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; ++i) {
                        n[5 * iter + i] = Double.parseDouble(fields[i]);
                        //System.out.println(n[10 * iter + i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
            }
            if (a.matches("%FLAG DIHEDRAL_PHASE")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Dihedral Phase Shift (in radians)");
                iter = 0;
                while (!a.contains("%") && a.length()>0) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; ++i) {
                        phi[5 * iter + i] = Double.parseDouble(fields[i]);
                        double foo = eUnit.toSim(k[5 * iter + i]);
                        if (phi[5 * iter + i] > 3.1) {
                            foo = -foo;
                        }
                        if (n[5 * iter + i] == 0) {
                            p4Bonds[5 * iter + i] = new P4BondTorsion(opts.space, foo, 0, 0, 0);
                            //System.out.println("Phi = Pi, n=0");
                        }
                        if (n[5 * iter + i] == 1) {
                            p4Bonds[5 * iter + i] = new P4BondTorsion(opts.space, 0, foo, 0, 0);
                            //System.out.println("Phi = Pi, n=1");
                        }
                        if (n[5 * iter + i] == 2) {
                            p4Bonds[5 * iter + i] = new P4BondTorsion(opts.space, 0, 0, -foo, 0);
                            //System.out.println("Phi = Pi, n=2");
                        }
                        if (n[5 * iter + i] == 3) {
                            p4Bonds[5 * iter + i] = new P4BondTorsion(opts.space, 0, 0, 0, foo);
                            //System.out.println("Phi = Pi, n=3");
                        }
                        //System.out.println("Phi = Pi");
                        //System.out.println(phi[5 * iter + i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
                for (double i =0; i < 3.2; i+=0.01) {
                    //System.out.println(i + "\t" + p4Bonds[0].u(Math.cos(i)) + "\t" + p4Bonds[1].u(Math.cos(i)) + "\t" + p4Bonds[2].u(Math.cos(i)));
                }
            }
            if (a.matches("%FLAG LENNARD_JONES_ACOEF")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Lennard Jones A coefficient (4*epsilon*sigma^12)");
                iter = 0;
                while (!a.contains("%")) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; ++i) {
                        A[5 * iter + i] = Double.parseDouble(fields[i]);
                        //System.out.println(A[5 * iter + i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
            }
            if (a.matches("%FLAG LENNARD_JONES_BCOEF")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Lennard Jones B coefficient (4*epsilon*sigma^6)");
                iter = 0;
                while (!a.contains("%")) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; ++i) {
                        B[5 * iter + i] = Double.parseDouble(fields[i]);
                        sigma[5 * iter + i] = Math.pow(A[5 * iter + i] / B[5 * iter + i], 1.0 / 6.0);
                        epsilon[5 * iter + i] = (B[5 * iter + i] * B[5 * iter + i]) / (4 * A[5 * iter + i]);
                        p2LennardJones[5 * iter + i] = new P2LennardJones(opts.space, sigma[5 * iter + i], eUnit.toSim(epsilon[5 * iter + i]));
                        //System.out.println(B[5 * iter + i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
            }
            if (a.matches("%FLAG BONDS_INC_HYDROGEN")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Bonds including Hydrogen");
                iter = 0;
                int[] elements = new int[3 * numBondsIH];
                while (!a.contains("%") && a.length()>0) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; ++i) {
                        elements[10 * iter + i] = Integer.parseInt(fields[i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
                for (int i = 0; i < elements.length; i = i + 3) {
                    bondedPairs[elements[i+2]-1].add(new int[]{elements[i] / 3, elements[i + 1] / 3});
                }
            }
            if (a.matches("%FLAG BONDS_WITHOUT_HYDROGEN")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Bonds without Hydrogen");
                iter = 0;
                int[] elements = new int[3 * numBondsWH];
                while (!a.contains("%") && a.length()>0) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; ++i) {
                        elements[10 * iter + i] = Integer.parseInt(fields[i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
                for (int i = 0; i < elements.length; i = i + 3) {
                    bondedPairs[elements[i+2]-1].add(new int[]{elements[i] / 3, elements[i + 1] / 3});
                }
            }
            if (a.matches("%FLAG ANGLES_INC_HYDROGEN")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Angles including Hydrogen");
                iter = 0;
                int[] elements = new int[4 * numAnglesIH];
                while (!a.contains("%") && a.length()>0) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; ++i) {
                        elements[10 * iter + i] = Integer.parseInt(fields[i]);
                        //System.out.println(elements[10*iter + i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
                for (int i = 0; i < elements.length; i = i + 4) {
                    bondedTriplets[elements[i + 3] - 1].add(new int[]{elements[i] / 3, elements[i + 1] / 3, elements[i + 2] / 3});
                }
            }
            if (a.matches("%FLAG ANGLES_WITHOUT_HYDROGEN")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Angles without Hydrogen");
                iter = 0;
                int[] elements = new int[4 * numAnglesWH];
                while (!a.contains("%") && a.length()>0) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; ++i) {
                        elements[10 * iter + i] = Integer.parseInt(fields[i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
                for (int i = 0; i < elements.length; i = i + 4) {
                    bondedTriplets[elements[i + 3] - 1].add(new int[]{elements[i] / 3, elements[i + 1] / 3, elements[i + 2] / 3});
                }
            }
            if (a.matches("%FLAG DIHEDRALS_INC_HYDROGEN")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Dihedrals including Hydrogen");
                iter = 0;
                int[] elements = new int[5 * numTorsionsIH];
                while (!a.contains("%") && a.length()>0) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; ++i) {
                        elements[10 * iter + i] = Integer.parseInt(fields[i]);
                        //System.out.println(elements[10*iter + i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
                for (int i = 0; i < elements.length; i = i + 5) {
                    elements[i+2] = Math.abs(elements[i+2]);
                    bondedQuads[elements[i + 4] - 1].add(new int[]{elements[i] / 3, elements[i + 1] / 3, elements[i + 2] / 3, elements[i + 3] / 3});
                }
            }
            if (a.matches("%FLAG DIHEDRALS_WITHOUT_HYDROGEN")) {
                linesIterator.next();
                a = linesIterator.next();
                //System.out.println("Dihedrals without Hydrogen");
                iter = 0;
                int[] elements = new int[5 * numTorsionsWH];
                while (!a.contains("%") && a.length()>0) {
                    //System.out.println(a);
                    String[] fields = a.split("\\s+");
                    for (int i = 0; i < fields.length; ++i) {
                        elements[10 * iter + i] = Integer.parseInt(fields[i]);
                    }
                    ++iter;
                    a = linesIterator.next();
                }
                for (int i = 0; i < elements.length; i = i + 5) {
                    elements[i+2] = Math.abs(elements[i+2]);
                    bondedQuads[elements[i + 4] - 1].add(new int[]{elements[i] / 3, elements[i + 1] / 3, elements[i + 2] / 3, elements[i + 3] / 3});
                    //System.out.println(bondedQuads[elements[i + 4] - 1]);
                }
            }
            //System.out.println(bondedPairs[0]);
        }
        //System.out.println(bondedPairs[0]);
        Element C = Carbon.INSTANCE;
        Element H = Hydrogen.INSTANCE;
        AtomType atomTypeC = new AtomType(C);
        AtomType atomTypeH = new AtomType(H);
        for (int i = 0; i < atomTypeId.length; ++i){
            if(atomTypeId[i] == 1)
                ++countC;
            else if(atomTypeId[i] == 2)
                ++countH;
        }
        //builder.addAtom(atomTypeC);
        //builder.addAtom(atomTypeH);
        builder.addCount(atomTypeC, countC);
        builder.addCount(atomTypeH, countH);
        builder.withConformation(ConformationFile.makeConformation(opts.space, smiles+".xyz"));
        builder.setDynamic(true);
        SpeciesGeneral speciesGeneral = builder.build();

        SpeciesManager speciesManager = SpeciesManager.builder().addSpecies(speciesGeneral).build();
        Box box = new Box(opts.space);
        box.addSpeciesNotify(speciesGeneral);
        PotentialMasterBonding.FullBondingInfo fullBondingInfo = new PotentialMasterBonding.FullBondingInfo(speciesManager);
        for (int i = 0; i < p2Bonds.length; ++i) {
            fullBondingInfo.setBondingPotentialPair(speciesGeneral, p2Bonds[i], bondedPairs[i]);
        }
        for (int i = 0; i < p3Bonds.length; ++i) {
            fullBondingInfo.setBondingPotentialTriplet(speciesGeneral, p3Bonds[i], bondedTriplets[i]);
        }
        for (int i = 0; i < p4Bonds.length; ++i) {
            fullBondingInfo.setBondingPotentialQuad(speciesGeneral, p4Bonds[i], bondedQuads[i]);
        }
        PotentialComputePair potentialComputePair = new PotentialComputePair(speciesManager, box, new NeighborManagerSimple(box));
        potentialComputePair.setPairPotential(atomTypeC, atomTypeC, p2LennardJones[0]);
        potentialComputePair.setPairPotential(atomTypeC, atomTypeH, p2LennardJones[1]);
        potentialComputePair.setPairPotential(atomTypeH, atomTypeH, p2LennardJones[2]);

        Potential2Soft [][] potential2Softs = new Potential2Soft[uniqueAtoms][uniqueAtoms];
        potential2Softs[0][0] = p2LennardJones[0];
        potential2Softs[0][1] = p2LennardJones[1];
        potential2Softs[1][0] = p2LennardJones[1];
        potential2Softs[1][1] = p2LennardJones[2];

        return new ParserAMBER.Stuff(fullBondingInfo, potential2Softs, speciesManager, box);
    }

    public static ParserAMBER.Stuff makeStuff(String smiles, Options opts) {
        FileReader fileReader;
        try {
            fileReader = new FileReader(smiles+".prmtop");
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+smiles+", caught IOException: " + e.getMessage());
        }
        List<String> lines = new ArrayList<String>();
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            String line = null;
            while ((line = bufReader.readLine()) != null) {
                lines.add(line.trim());
            }
            bufReader.close();
            //System.out.println(lines);
            return makeStuffFromLines(smiles, lines, opts);
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    public static void main(String[] args) {
        ParserAMBER.makeStuff("CCCCCC", new Options());
    }

    public enum Truncation {NONE, TRUNCATED, SHIFTED, FORCE_SHIFTED, SWITCHED}

    public static class Stuff {
        public final PotentialMasterBonding.FullBondingInfo fullBondingInfo;
        public final Potential2Soft [][] potential2Soft;
        public final SpeciesManager speciesManager;
        public final Box box;
        public Stuff(PotentialMasterBonding.FullBondingInfo fullBondingInfo, Potential2Soft [][] potential2Soft, SpeciesManager speciesManager, Box box) {
            this.fullBondingInfo = fullBondingInfo;
            this.potential2Soft = potential2Soft;
            this.speciesManager = speciesManager;
            this.box = box;
        }
    }

    public static class Options {
        public Space space = Space3D.getInstance();
        public Truncation truncation = Truncation.NONE;
        public double rc = Double.POSITIVE_INFINITY;
    }
}
