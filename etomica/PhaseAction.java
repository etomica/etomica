package etomica;

//import java.awt.event.ActionEvent;

 /**
  * Superclass of classes that apply some elementary action (transformation) to a phase.
  * 
  */
 
 /* History
  * 08/27/03 (DAK) modifications to ImposePBC to handle molecules versus atoms
  */
public abstract class PhaseAction extends etomica.Action implements PhaseListener {

    public static String getVersion() {return "PhaseAction:01.03.28/"+Action.VERSION;}
    private static final AtomIteratorMolecule moleculeIterator = new AtomIteratorMolecule();

    protected Phase phase;
    public PhaseAction() {this(null);}
    public PhaseAction(Phase p) {
        super();
        setPhase(p);
    }
        
    public void setPhase(Phase p) {
    	if(p == null) throw new IllegalArgumentException("Cannot set null phase in PhaseAction");
    	phase = p;
    }
    public Phase getPhase() {return phase;}
    
    public void actionPerformed(SimulationEvent evt) {
        actionPerformed((PhaseEvent)evt);
    }
    public void actionPerformed(PhaseEvent pe) {
        actionPerformed(pe.phase());
    }
    
    public void actionPerformed() {actionPerformed(phase);}
    
    public abstract void actionPerformed(Phase p);

        
    //************* End of fields and methods for PhaseAction ************//
    
    //************* The remainder defines some subclasses of PhaseAction **********//
    
    public static final class Translate extends PhaseAction implements Action.Undoable {
        
        private Space.Vector translationVector;
        
        public Translate(Phase p) {
            super(p);
        }
        
        public void setTranslationVector(Space.Vector v) {translationVector.E(v);}
        public Space.Vector getTranslationVector() {return translationVector;}
        
        public void actionPerformed(Phase p) {doAction(p, translationVector);}
        
        public void actionPerformed(Space.Vector v) {doAction(phase, v);}
        
        public static void doAction(Phase p, Space.Vector v) {
            if(v == null || p == null) return;
            p.speciesMaster().coord.translateBy(v);
        }
        public void attempt() {doAction(phase, translationVector);}
        public void undo() {
            translationVector.TE(-1);
            doAction(phase, translationVector);
            translationVector.TE(-1);
        }
    }//end of Translate
    
    public static final class ImposePbc extends PhaseAction implements Integrator.IntervalListener.ImposePbc {
        
        
        public ImposePbc(Phase p) {
            super(p);
        }
        
        public void actionPerformed() {
			Space.Boundary boundary = phase.boundary();
			iterator.reset();
			while(iterator.hasNext()) {
				Atom a = iterator.nextAtom();
				//centralImage returns true if it causes the atom to be moved
//				if(boundary.centralImage(a.coord)) a.seq.moveNotify();
				boundary.centralImage(a.coord); //08/27/03 - don't need to notify sequencer because coord.translateBy does this
			}        	
        }
               
        public void actionPerformed(Phase p) {
        	setPhase(p);
        	actionPerformed();
		}
        
        public void intervalAction(Integrator.IntervalEvent evt) {
            actionPerformed();
        }
        
        private AtomIterator iterator;
        private boolean applyToMolecules = false;
        
        public void setPhase(Phase phase) {
        	super.setPhase(phase);
        	if(applyToMolecules) iterator = phase.makeMoleculeIterator();
        	else iterator = phase.makeAtomIterator();
        }
		/**
		 * Returns the iterator that gives the atoms to which central imaging
		 * is applied.
		 * @return AtomIteratorList
		 */
		public AtomIterator getIterator() {
			return iterator;
		}

		/**
		 * Sets the iterator the gives the atoms for central imaging.  Normally
		 * this does not need to be set, but if central imaging scheme is
		 * desired for application at a level other than the molecules or the
		 * leaf atoms, or if it is to be applied only to a subet of the atoms
		 * in a phase, this can be invoked by setting this iterator
		 * appropriately.
		 * @param iterator The iterator to set
		 */
		public void setIterator(AtomIterator iterator) {
			this.iterator = iterator;
		}

		/**
		 * Returns the applyToMolecules.
		 * @return boolean
		 */
		public boolean isApplyToMolecules() {
			return applyToMolecules;
		}

		/**
		 * Sets a flag indicating whether periodic boundaries are applied to the
		 * molecules (true), or to the atoms (false).  If applied to the atoms
		 * (the default case), then central imaging is done to each atom
		 * individually, which could cause a molecule to be split, with some of
		 * its atoms on one edge of the simulation box, and others on the other
		 * edge.  If applied to molecules, the entire molecule will be shifted
		 * as a whole when enforcing central imaging.
		 * @param applyToMolecules The new value of the flag.
		 */
		public void setApplyToMolecules(boolean applyToMolecules) {
			this.applyToMolecules = applyToMolecules;
			setPhase(phase);
		}

    }//end of ImposePBC
    
    
    /**
     * Performs actions that cause volume of system to expand, with particle positions
     * scaled to keep them in the same relative positions.  Inflation is done isotropically,
     * equal in all coordinate directions.
     */
    public static final class Inflate extends PhaseAction implements Action.Undoable {
        
        double scale = 1.0;
        
        public Inflate() {
            super();
        }
        public Inflate(Phase p) {
            super(p);
        }
                
        public void setScale(double s) {scale = s;}
        public double getScale() {return scale;}

        public void actionPerformed(double s) {
            setScale(s);
            doAction(phase, scale);
        }
        public void actionPerformed(Phase p) {doAction(p, scale);}

        /**
         * Performs isotropic inflation.
         */
        public static void doAction(Phase phase, double scale) {
            phase.boundary().inflate(scale);
            moleculeIterator.setPhase(phase);
            moleculeIterator.reset();
            while(moleculeIterator.hasNext()) {
                moleculeIterator.nextAtom().coord.inflate(scale);
            }
        }
        
        public void attempt() {
            doAction(phase, scale);
        }
        public void undo() {
            doAction(phase, 1.0/scale);
        }
    }//end of Inflate 

    /**
     * Performs actions that cause volume of system to expand, with particle positions
     * scaled to keep them in the same relative positions.  Inflation is done anisotropically,
     * so that each dimension can be scaled differently.
     */
    public static class InflateAnisotropic extends PhaseAction implements Action.Undoable {
        
        protected Space.Vector scale;
        protected transient Space.Vector temp;
        private transient Space.Vector oldDimensions;
        
        public InflateAnisotropic(Phase p) {
            super(p);
            if(p != null) {
                scale = p.simulation().space().makeVector();
                temp = p.simulation().space().makeVector();
                oldDimensions = p.simulation().space().makeVector();
            }
        }
                
        public void setScale(Space.Vector s) {scale = s;}
        public Space.Vector getScale() {return scale;}

        public void setPhase(Phase p) {
            super.setPhase(p);
            if(p != null && temp == null) temp = p.simulation().space().makeVector();
        }
 
        public void actionPerformed(Phase p) {doAction(p, scale, temp);}

        /**
         * Performs anisotropic inflation.
         * @param scale Vector describing scaling to be performed in each coordinate
         *              direction.  scale(i) == 1 indicates no change in size in direction i.
         */
        public static void doAction(Phase phase, Space.Vector scale, Space.Vector work) {
            phase.boundary().inflate(scale);
      //      scale.PE(-1.0);
            moleculeIterator.setPhase(phase);
            moleculeIterator.reset();
            while(moleculeIterator.hasNext()) {
                moleculeIterator.nextAtom().coord.inflate(scale);
              //  Atom m = iterator.next();
              //  work.E(m.coord.position());
              //  work.TE(scale);
              //  m.coord.displaceBy(work);
            }
            //scale.PE(1.0);
        }

        public void attempt() {
    //        oldDimensions.E(phase.boundary().dimensions());
            doAction(phase, scale, temp);
        }
        public void undo() {
            temp.E(1.0);
            temp.DE(scale);
            doAction(phase, temp, null);
     //       phase.boundary().setDimensions(oldDimensions);
       //     iterator.setBasis(phase);
       //     iterator.reset();
       //     while(iterator.hasNext()) {
       //         Atom m = iterator.next();
    //            m.coord.replace();
           // }
        }
    }//end of InflateAnisotropic 
    
    /**
     * Toggles the state of enabling the long-range-corrections to the 
     * potential truncation in the phase.
     */
    public static class ToggleLrc extends PhaseAction {
        
        public void actionPerformed(Phase phase) {
            phase.setLrcEnabled(!phase.isLrcEnabled());
        }
    }

}//end of PhaseAction