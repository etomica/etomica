/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Formatter;

/** 
 * 
 * This class requires a bash script called ./runQChem to load the appropriate Q-Chem module and run it.  
 * 
 * This class cannot be used within Eclipse (or on our network in general) because it does not have access to Q-Chem.
 * 
 * Currently, it writes an input file appropriate for particular description of argon.  
 * Some writing to the Q-Chem input file will always be required because the positions of the molecules must be updated.
 * 
 * @author kate (adapted from Tai Boon Tan's P2DLPOLY)
 *
 */


public class P2QChem extends PotentialMolecular {
	
	public P2QChem (Space space){
		super(2, space);
    
	}
	
	public double energy(IMoleculeList atoms) {

		//System.out.println("Made it to P2QChem.energy()");
        
    	IMolecule molecule1 = atoms.get(0);
    	IMolecule molecule2 = atoms.get(1);
    	
    	Vector atomPos1 = molecule1.getChildList().get(0).getPosition();
    	Vector atomPos2  = molecule2.getChildList().get(0).getPosition();
      
    	Vector r12Vec = space.makeVector();
    	r12Vec.Ev1Mv2(atomPos1, atomPos2);
    	
    	double r12Mag = Math.sqrt(r12Vec.squared());
        	
    	if (r12Mag < 2.5) {
    		//System.out.println(""+r12Mag);
    		return 1e7;
    		
    	} else  {
    		
    		makeInputFile(atomPos1, atomPos2);
    		runQChem();
    		
    		File file = new File("argon.out");	
    		while (!file.exists()) {
    			
    			try {
				
    				Thread.sleep(1000);
				
    			} catch (InterruptedException err){
					
				}
				
    		}
    		
    		//System.out.println("argon.out created");
    		
    		double energy = Double.NaN;
    		
    		while (Double.isNaN(energy)) {
    			
    		
	    		try{	
	    			
	    			FileReader fileReader = new FileReader(file);			
	    			BufferedReader bufReader = new BufferedReader(fileReader);
	    			
	    			
	    			String line;
	    			
	    			while ((line = bufReader.readLine()) != null) {
	    				if (line.indexOf("Convergence criterion met") > -1) {
	    					String string = line.split(" +")[2];
	    					energy = Double.parseDouble(string); 
	    				    break;	
	    				}
	    				
	    			}
	    		

	    		}catch (IOException e){
	    			
	    			throw new RuntimeException(e);
	    		}
	    		
	    		if (Double.isNaN(energy)) {
	    			try {
						
	    				Thread.sleep(1000);
					
	    			} catch (InterruptedException err){
						
					}
	    		}
	    		
    		}
    		
    		//System.out.println("Found energy");
	    		
			File input = new File ("argon.in");
			input.renameTo(new File ("argon.in.old"));
			
			File output = new File ("argon.out");
			output.renameTo(new File ("argon.out.old"));
			
			//System.out.println("r12 (Angstroms) = " + r12Mag);
			//System.out.println("ut0t (Hartrees) = " + energy);
			
			energy = energy - 2*(-529.7492732249);
			
			energy = energy*2625.5; //convert to kJ/mol
			energy = energy*1000; //convert to J/mol
			double energyOverR = energy/8.314; // Kelvin
			//System.out.println("u12/k (K) = " + energyOverR);
			return energyOverR;
    	}
			
		
	}
	
	public void makeInputFile(Vector atomPos1, Vector atomPos2) {
		
		String fileName = "argon.in";
		
		try {
    	        	
			Formatter formatter = new Formatter(new File(fileName));
	        
        	formatter.format("%s\n", "$molecule");
        	formatter.format("%s\n", " +0 1");
			
			formatter.format("%s\t", new Object[]{" Ar"});
			formatter.format("%20.12f%20.12f%20.12f\n", new Object[]{atomPos1.getX(0), atomPos1.getX(1), atomPos1.getX(2)});
			
			formatter.format("%s\t", new Object[]{" Ar"});
			formatter.format("%20.12f%20.12f%20.12f\n", new Object[]{atomPos2.getX(0), atomPos2.getX(1), atomPos2.getX(2)});
        
            formatter.format("%s\n\n", "$end");
            
    		formatter.format("%s\n", "$rem");
    		
    		// contents of Mareks arg-33.in file
    		formatter.format("%s\n", "jobtype                 sp");
    		formatter.format("%s\n", "ideriv                  1");
    		formatter.format("%s\n", "exchange                PW86");
    		formatter.format("%s\n", "correlation             PBE");
    		formatter.format("%s\n", "basis                   6-31+G*");
    		formatter.format("%s\n", "basis_lin_dep_thresh    5");
    		formatter.format("%s\n", "scf_guess               sad");
    		formatter.format("%s\n", "scf_final_print         0");
    		formatter.format("%s\n", "scf_convergence         8");
    		formatter.format("%s\n", "symmetry                false");
    		formatter.format("%s\n", "sym_ignore              1");
    		formatter.format("%s\n", "thresh                  14");
    		formatter.format("%s\n", "incdft                  0");
    		formatter.format("%s\n", "xc_grid                 000120000302");
    		formatter.format("%s\n", "dftvdw_jobnumber        1");
    		formatter.format("%s\n", "dftvdw_method           2");
    		formatter.format("%s\n", "dftvdw_use_ele_drv      1");
    		formatter.format("%s\n", "dftvdw_d2xmethod        1");
    		formatter.format("%s\n", "dftvdw_d2xcompute       1");
    		formatter.format("%s\n", "dftvdw_print            1");
    		formatter.format("%s\n", "dftvdw_alpha1           75");
    		formatter.format("%s\n", "dftvdw_alpha2           125");
    		formatter.format("%s\n", "mem_static              516");
    		formatter.format("%s\n", "mem_total               2000");
    		
    		 
    		formatter.format("%s\n", "$end");
            formatter.close();

    		/*formatter.format("%s\n", "   exchange                b3lyp");
    		formatter.format("%s\n", "   basis                   6-311+G(d,p)");
    		formatter.format("%s\n", "   thresh                  14");
    		formatter.format("%s\n", "   xc_grid                 1");
    		formatter.format("%s\n", "   incdft                  0");
    		formatter.format("%s\n", "   dftvdw_jobnumber        2");
    		formatter.format("%s\n", "   dftvdw_method           1");*/
    		

        	
        } catch(IOException e) {
            System.out.println("Problem writing to "+fileName+", caught IOException: " + e.getMessage());
        }
	}
	
public void runQChem() {
		
		try{
			
			Runtime rt0 = Runtime.getRuntime();
			
			String[] args = new String[]{"/san/user/schadel2/u2/Q-Chem/argon/DirectSampling/dft/PW86PBE/631G/runQChem"};
			//System.out.println("qchem argon.in argon.out");
			//String[] args0 = new String[]{"qchem", "argon.in", "argon.out"};
			
			Process proc0 = rt0.exec(args);
			
			proc0.waitFor();
			
			
		}catch (IOException e){
			System.out.println("Problem running Q-Chem.");
			throw new RuntimeException(e);
		}catch (InterruptedException err){
			System.out.println("Problem running Q-Chem.");
			throw new RuntimeException(err);
		}
		
	}


	public double getRange() {
		return Double.POSITIVE_INFINITY;
	}


	public void setBox(Box box) {
        this.box = box;
    }
	
	protected Boundary boundary;

	private Box box;
	

}

