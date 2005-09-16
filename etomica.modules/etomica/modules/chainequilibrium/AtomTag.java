/*
 * Created on Apr 27, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.modules.chainequilibrium;


/**
 * @author Matt Moynihan
 *	This class contains a boolean which marks an atom if it has been 
 *counted as part of a molecule, it contains a  Three states, White, Gray & black,
 * white is unprocessed, gray is being processed, black has been processed. all atoms 
 * start as white 
 */
public class AtomTag{ //implements Atom.AgentSource{

	public boolean gray;	// This is if the atom has been processed
	public boolean black; 	// This is if the atom is BEING processed
	
	public AtomTag(){
		gray = false;
		black = false;
	}
	
// This function returns true if the atom is white
	public boolean white(){
		if((!(gray)) && (!(black))){
			return true;
		}else{
			return false; 
		}
	}

//	 This function returns true if the atom is gray
	public boolean gray(){
		if((gray)&&(!(black))){
			return true;
		}else{
			return false; 
		}
	}

//	 This function returns true if the atom is black	
	public boolean black(){
		if(black){
			return true;
		}else{
			return false;
		}
	}

//	processes the atom, if white --> gray, if gray --> black
	public void process(){	
		if((!(gray)) && (!(black))){
			gray = true;
		}else{
		if((gray)&&(!(black))){
		    black = true;
		}}
	}

//	resets the atom to white
	public void reset(){
		gray = false;
		black = false;
	}
}
