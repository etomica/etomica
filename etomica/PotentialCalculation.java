package etomica;

public abstract class PotentialCalculation implements AtomSetActive,
											 PhaseActive,
											 AtomActive,
											 AtomPairActive,
											 Atom3Active {

	protected Potential0 potential0;
	protected Potential1 potential1;
	protected Potential2 potential2;
	protected Potential3 potential3;
	public final AtomPairActive.InnerWrapper wrapper;
	public final AtomPairActive.OuterWrapper outerWrapper;

	public PotentialCalculation() {
		this(Simulation.instance.space);
	}
	public PotentialCalculation(Space space) {
		wrapper = new AtomPairActive.InnerWrapper(this, new AtomPair(space));
		outerWrapper = new AtomPairActive.OuterWrapper(this);
	}
 	
	public PotentialCalculation set(Potential0 p0) {
		potential0 = p0;
		return this;
	}
	public PotentialCalculation set(Potential1 p1) {
		potential1 = p1;
		return this;
	}
	public PotentialCalculation set(Potential2 p2) {
		potential2 = p2;
		return this;
	}
	public PotentialCalculation set(Potential3 p3) {
		potential3 = p3;
		return this;
	}
	
//	/**
//	 * Performs the atom-action on the species master of the phase.
//	 * @see etomica.PhaseActive#actionPerformed(Phase)
//	 */
	public void actionPerformed(Phase phase) {
//		actionPerformed(phase.speciesMaster);//problem with lrc in soft md if doing this
	}
	
	public void actionPerformed(AtomSet atomSet) {
		switch(atomSet.nBody()) {
			case 1: actionPerformed((Atom)atomSet); break;
			case 2: actionPerformed((AtomPair)atomSet); break;
			case 3: actionPerformed((Atom3)atomSet); break;
		}
	}
 	   
	public abstract void actionPerformed(Atom atom);
    
	public abstract void actionPerformed(AtomPair pair);
	
	public abstract void actionPerformed(Atom3 atom3);

	/**
	 * @see etomica.AtomPairActive#innerWrapper()
	 */
	public final AtomPairActive.InnerWrapper innerWrapper() {return wrapper;}
	/**
	 * @see etomica.AtomPairActive#outerWrapper()
	 */
	public final AtomPairActive.OuterWrapper outerWrapper() {return outerWrapper;}

	
	public final PotentialGroupWrapper wrapper() {
		if(firstWrapper == null) {
			return new PotentialGroupWrapper(this);
		} else {
			PotentialGroupWrapper wrapper = firstWrapper;
			firstWrapper = wrapper.next;
			wrapper.potentialCalculation = this;
			return wrapper;
		}		
	}
	
	private PotentialGroupWrapper firstWrapper = null;
	
	public static final class PotentialGroupWrapper extends PotentialCalculation {
		private PotentialGroupWrapper next;
		public final IteratorDirective localDirective = new IteratorDirective();
		private PotentialGroup potential;
		private PotentialCalculation potentialCalculation;
		
		PotentialGroupWrapper(PotentialCalculation calculation) {
			potentialCalculation = calculation;
		}
		
		public PotentialGroupWrapper set(PotentialGroup pg) {
			potential = pg;
			return this;
		}
		
		public void release() {
			next = potentialCalculation.firstWrapper;
			potentialCalculation.firstWrapper = this;
		}
		public void actionPerformed(AtomSet atomSet) {
			if(potential.potentialTruncation != null && potential.potentialTruncation.isZero(atomSet)) return;                
	            
			//if the atom of the pair is the one specified for calculation, then
			//it becomes the basis for the sub-potential iterations, and is no longer
			//specified to them via the iterator directive
			if(atomSet.contains(localDirective.atom1())) localDirective.set();
	            
			//loop over sub-potentials
			for(PotentialLinker link=potential.first; link!=null; link=link.next) {
				if(localDirective.excludes(link.potential)) continue; //see if potential is ok with iterator directive
				link.potential.calculate(atomSet, localDirective, potentialCalculation);
			}//end for	    	
		}

		/**
		 * @see etomica.AtomActive#actionPerformed(etomica.Atom)
		 */
		public void actionPerformed(Atom atom) {this.actionPerformed((AtomSet)atom);}

		/**
		 * @see etomica.Atom3Active#actionPerformed(etomica.Atom3)
		 */
		public void actionPerformed(Atom3 atom3) {this.actionPerformed((AtomSet)atom3);}

		/**
		 * @see etomica.AtomPairActive#actionPerformed(etomica.AtomPair)
		 */
		public void actionPerformed(AtomPair pair) {this.actionPerformed((AtomSet)pair);}

	}
	    
    public interface Summable {
        public double sum();
        public Summable reset();
    }
}//end of PotentialCalculation