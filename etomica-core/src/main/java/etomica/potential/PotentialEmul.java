/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

import etomica.atom.IAtom;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.api.IMoleculeList;
import etomica.space.Vector;
import etomica.atom.MoleculeArrayList;
import etomica.chem.elements.Argon;
import etomica.chem.elements.Helium;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Calorie;
import etomica.units.Kelvin;
import etomica.units.Mole;
import etomica.units.Prefix;
import etomica.units.PrefixedUnit;
import etomica.units.Unit;
import etomica.units.UnitRatio;

/**
 * PotentialEmul class invokes the emul external program to compute energies
 * between molecules.  This class handles only mono-atomic molecules.  The
 * atomic positions are written to an emul input file constructed from a
 * template. 
 */
public class PotentialEmul extends PotentialMolecular {

    protected String templateName;
    public static final Unit emulEnergyUnit = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT), Mole.UNIT);
    protected double r2Core;
    protected final Vector r12Vec;
    
	public PotentialEmul(Space space, String templateName){
	    this(space, templateName, 2.5);
	}
	
	public PotentialEmul(Space space, String templateName, double rCore) {
	    this(space, templateName, countMolecules(templateName), rCore);
	}

	protected PotentialEmul(Space space, String templateName, int nBody, double rCore) {
		super(nBody, space);
		this.r2Core = rCore*rCore;
        r12Vec = space.makeVector();
		this.templateName = templateName;
	}
	
	protected static int countMolecules(String templateName) {
        int nBody = 0;
        try {
            FileReader fileReader = new FileReader(templateName);
            BufferedReader bufReader = new BufferedReader(fileReader);
            String line = null;
            
            while ((line = bufReader.readLine()) != null) {
                if (line.matches("^%molecule$")) nBody++;
            }
            bufReader.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
        return nBody;
	}
	
	public double energy(IMoleculeList molecules) {
        
    	for (int i=0; i<molecules.getMoleculeCount(); i++) {
    	    IAtom atomi = molecules.getMolecule(i).getChildList().getAtom(0);
            for (int j=i+1; j<molecules.getMoleculeCount(); j++) {
                IAtom atomj = molecules.getMolecule(j).getChildList().getAtom(0);
                r12Vec.Ev1Mv2(atomi.getPosition(), atomj.getPosition());
                boundary.nearestImage(r12Vec);
                
            	if (r12Vec.squared() < r2Core) {
            		//System.out.println(""+r12Mag);
            		return Double.POSITIVE_INFINITY;
            	}
            }
    	}

		makeInputFile(molecules);
		runEmul();

		double energy = 0;
		try {
		    energy = readOutputFile();
		}
		catch (RuntimeException ex) {
		    System.err.println("failed to get energy.");
		    System.err.println("Here are the atomic coordinates:");
	        for (int i=0; i<molecules.getMoleculeCount(); i++) {
	            IAtom atomi = molecules.getMolecule(i).getChildList().getAtom(0);
	            System.out.println(i+" "+atomi.getPosition());
	        }
	        throw new RuntimeException(ex);
		}
		
		return energy;
	}
	
	protected double readOutputFile() {
	    File file = new File("calc.out");   
	    if (!file.exists()) {
	        throw new RuntimeException("could not find calc.out");
	    }
	        
	    double energy = Double.NaN;

	    try{    
	        FileReader fileReader = new FileReader(file);           
	        BufferedReader bufReader = new BufferedReader(fileReader);
	            
	        String line;

	        while ((line = bufReader.readLine()) != null) {
	        	// System.out.println(line);
	            if (line.indexOf("kcal/mol") > -1) {
	                String[] strings = line.split(" +");
	                String stringVal = strings[strings.length-1];
	                energy = Double.parseDouble(stringVal); 
	                break;  
	            }
	        }
	        bufReader.close();
	    }
	    catch (IOException e){
	        throw new RuntimeException(e);
	    }
	    if (Double.isNaN(energy)) {
	        System.err.println("could not find energy in calc.out");
	        System.err.println("this is was calc.out looked like:");
            System.err.println("==============");
	        try{    
	            FileReader fileReader = new FileReader(file);           
	            BufferedReader bufReader = new BufferedReader(fileReader);
	                
	            String line;

	            while ((line = bufReader.readLine()) != null) {
	                System.out.println(line);
	            }
	            bufReader.close();
	        }
	        catch (IOException e){
	            throw new RuntimeException(e);
	        }
            System.err.println("==============");

	        throw new RuntimeException();
	    }
	    return emulEnergyUnit.toSim(energy);
	}
	
	public void makeInputFile(IMoleculeList molecules) {
	    
	    try {
    	    FileReader fileReader = new FileReader(templateName);
            BufferedReader bufReader = new BufferedReader(fileReader);
    	    FileWriter fw = new FileWriter("calc.in", false);
    	    String line = null;
    	    int nAtomsWritten = 0;
    	    
            while ((line = bufReader.readLine()) != null) {
                fw.write(line+"\n");
                if (line.matches("^%molecule$")) {
                    // now "0 1"
                    fw.write(bufReader.readLine()+"\n");
                    // the next line has the element and we add the x-y-z coordinate
                    // XXX we only handle mono-atomic molecules
                    IAtom atom = molecules.getMolecule(nAtomsWritten).getChildList().getAtom(0);
                    Vector p = atom.getPosition();
                    fw.write(bufReader.readLine()+" "+p.getX(0)+" "+p.getX(1)+" "+p.getX(2)+"\n");
                    nAtomsWritten++;
                }
                else if (line.matches("^%cluster$")) {
                    // now "0 1"
                    fw.write(bufReader.readLine()+"\n");
                    // the following lines have the element and we add the x-y-z coordinate
                	for (int i=0; i<molecules.getMoleculeCount(); i++) {
                        IAtom atom = molecules.getMolecule(i).getChildList().getAtom(0);
                	    Vector p = atom.getPosition();
                        fw.write(bufReader.readLine()+" "+p.getX(0)+" "+p.getX(1)+" "+p.getX(2)+"\n");
                	}
                }
            }
            bufReader.close();
            fw.close();
	    }
	    catch (IOException e) {
	        throw new RuntimeException(e);
	    }
	}
	
	public static String execPath = "/opt/emul/emul";
	
	protected void runEmul() {
		try{
	        File file = new File("calc.out");   
	        if (file.exists()) {
	            if (!file.delete()) {
	                throw new RuntimeException("unable to clear out calc.out");
	            }
	        }
            ProcessBuilder pb = new ProcessBuilder(execPath, "calc.in");
            pb.redirectErrorStream(true);
            Process proc = pb.start();
            BufferedReader reader = new BufferedReader(new InputStreamReader(proc.getInputStream()));

            boolean unhappy = false;
            String line = reader.readLine();
            if (line != null) {
                unhappy = true;
                System.err.println("emul produced unexpected output:");
            }
            while (line!= null) {
                System.err.println(line);
                line = reader.readLine();
            }
            
            proc.waitFor();
            if (unhappy) throw new RuntimeException("emul was unhappy");
		}
		catch (IOException e){
			System.out.println("Problem running Emul.");
			throw new RuntimeException(e);
		}
		catch (InterruptedException err){
			System.out.println("Problem running Emul.");
			throw new RuntimeException(err);
		}
	}

	public double getRange() {
		return Double.POSITIVE_INFINITY;
	}


	public void setBox(Box box) {
        this.boundary = box.getBoundary();
    }
	
	protected Boundary boundary;
	
	public static void main(String[] args) {
	    Space space = Space3D.getInstance();
	    if (false) {
    	    PotentialEmul p2 = new PotentialEmul(space, "template.in");
    	    Simulation sim = new Simulation(space);
    	    Box box = new Box(new BoundaryRectangularNonperiodic(space), space);
    	    sim.addBox(box);
    	    SpeciesSpheresMono species = new SpeciesSpheresMono(space, Helium.INSTANCE);
    	    sim.addSpecies(species);
    	    box.setNMolecules(species, 2);
    	    p2.setBox(box);
    	    P2LennardJones p2LJ = new P2LennardJones(space, 3.1, Kelvin.UNIT.toSim(98));
    	    p2LJ.setBox(box);
    	    P2QChem p2Q = new P2QChem(space);
    	    p2Q.setBox(box);
    	    P2HePCKLJS p2He = new P2HePCKLJS(space);
    	    p2He.setBox(box);
    	    for (int i=10; i<100; i++) {
    	    	double r = i*0.1;
    		    box.getLeafList().getAtom(1).getPosition().setX(0, r);
    		    double uEmul = p2.energy(box.getMoleculeList(species));
    		    double uLJ = p2LJ.energy(box.getLeafList());
    		    double uQ = 0; //p2Q.energy(box.getMoleculeList(species));
    		    double uHe = p2He.energy(box.getLeafList());
    		    System.out.println(r+" "+uEmul+" "+uLJ+" "+uQ+" "+uHe);
    	    }
	    }

	    if (true) {
	        PotentialEmul p3 = new PotentialEmul(space, "3body_template.in");
            PotentialEmul p2 = new PotentialEmul(space, "2body_template.in");
	        Simulation sim = new Simulation(space);
	        Box box = new Box(new BoundaryRectangularNonperiodic(space), space);
	        sim.addBox(box);
	        SpeciesSpheresMono species = new SpeciesSpheresMono(space, Argon.INSTANCE);
	        sim.addSpecies(species);
	        box.setNMolecules(species, 3);
            p2.setBox(box);
	        p3.setBox(box);
            box.getLeafList().getAtom(0).getPosition().setX(0, -4);
            box.getLeafList().getAtom(1).getPosition().setX(0, -4);
            box.getLeafList().getAtom(0).getPosition().setX(1, -3);
            box.getLeafList().getAtom(1).getPosition().setX(1, 3);
            MoleculeArrayList molecules02 = new MoleculeArrayList(2);
            molecules02.add(box.getMoleculeList().getMolecule(0));
            molecules02.add(box.getMoleculeList().getMolecule(1));
            double u01 = p2.energy(molecules02);
            molecules02.clear();
            molecules02.add(box.getMoleculeList().getMolecule(0));
            molecules02.add(box.getMoleculeList().getMolecule(2));
            MoleculeArrayList molecules12 = new MoleculeArrayList(2);
            molecules12.add(box.getMoleculeList().getMolecule(1));
            molecules12.add(box.getMoleculeList().getMolecule(2));
            
	        for (int i=0; i<100; i++) {
	            double r = -3 + i*0.1;
	            box.getLeafList().getAtom(2).getPosition().setX(0, r);
	            double uEmul = p3.energy(box.getMoleculeList(species));
                double u02 = p2.energy(molecules02);
	            double u12 = p2.energy(molecules12);
	            System.out.printf("%5.2f ", r);
	            System.out.println(uEmul+" "+(uEmul - (u01+u02+u12)));
//	            System.out.println("  "+u01+" "+u02+" "+u12);
	        }
	    }
	}
}

