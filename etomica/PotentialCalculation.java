package etomica;

public class PotentialCalculation implements AtomSetActive,
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
	public void actionPerformed(Phase phase) {}
	
	public void actionPerformed(AtomSet atomSet) {
		switch(atomSet.nBody()) {
			case 1: actionPerformed((Atom)atomSet); break;
			case 2: actionPerformed((AtomPair)atomSet); break;
			case 3: actionPerformed((Atom3)atomSet); break;
		}
	}
 	   
	public void actionPerformed(Atom atom) {}
    
	public void actionPerformed(AtomPair pair) {}
	
	public void actionPerformed(Atom3 atom3) {}

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
			return new PotentialGroupWrapper();
		} else {
			PotentialGroupWrapper wrapper = firstWrapper;
			firstWrapper = wrapper.next;
			return wrapper;
		}		
	}
	
	private PotentialGroupWrapper firstWrapper = null;
	
	public final class PotentialGroupWrapper extends PotentialCalculation {
		private PotentialGroupWrapper next;
		public final IteratorDirective localDirective = new IteratorDirective();
		private PotentialGroup potential;
		
		public PotentialGroupWrapper set(PotentialGroup pg) {
			potential = pg;
			return this;
		}
		
		public void release() {
			next = firstWrapper;
			firstWrapper = this;
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
				link.potential.calculate(atomSet, localDirective, PotentialCalculation.this);
			}//end for	    	
		}

	}
	    
    public interface Summable {
        public double sum();
        public Summable reset();
    }
}//end of PotentialCalculation