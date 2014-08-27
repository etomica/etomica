package etomica.potential;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.api.IAtom;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.chem.elements.Helium;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.ISpace;
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
    protected final IVectorMutable r12Vec;
    protected int nBody;
    
	public PotentialEmul(ISpace space, String templateName){
	    this(space, templateName, 2.5);
	}
	
	public PotentialEmul(ISpace space, String templateName, double rCore) {
	    this(space, templateName, countMolecules(templateName), rCore);
	    
	}
	protected PotentialEmul(ISpace space, String templateName, int nBody, double rCore) {
		super(nBody, space);
		this.r2Core = rCore*rCore;
        r12Vec = space.makeVector();
		this.templateName = templateName;
	}
	
	private static int countMolecules(String templateName) {
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

		double energy = readOutputFile();
		
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
	        	System.out.println(line);
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
	        throw new RuntimeException("could not find energy in calc.out");
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
                    fw.write("0 1\n");
                    // XXX we only handle mono-atomic atoms
                    IAtom atom = molecules.getMolecule(nAtomsWritten).getChildList().getAtom(0);
                    IVector p = atom.getPosition();
                    fw.write(atom.getType().getElement().getSymbol()+" "+p.getX(0)+" "+p.getX(1)+" "+p.getX(2)+"\n");
                    nAtomsWritten++;
                }
                else if (line.matches("^%cluster$")) {
                	fw.write("0 1\n");
                	for (int i=0; i<molecules.getMoleculeCount(); i++) {
                        IAtom atom = molecules.getMolecule(i).getChildList().getAtom(0);
                	    IVector p = atom.getPosition();
                	    fw.write(atom.getType().getElement().getSymbol()+" "+p.getX(0)+" "+p.getX(1)+" "+p.getX(2)+"\n");
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
	
	public static String execPath = "/home/fenglai/emul/trunk/src/debug_test/virial/emul";
	
	protected void runEmul() {
		try{
			Runtime rt0 = Runtime.getRuntime();
			String[] args = new String[]{execPath,"calc.in"};
			Process proc0 = rt0.exec(args);
			proc0.waitFor();
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


	public void setBox(IBox box) {
        this.boundary = box.getBoundary();
    }
	
	protected IBoundary boundary;
	
	public static void main(String[] args) {
	    ISpace space = Space3D.getInstance();
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
}

