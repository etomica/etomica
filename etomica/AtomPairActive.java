package etomica;

/**
 * @author David Kofke
 *
 * Interface for a class that can perform an action on an AtomPair.
 */

/* History
 * 
 * 01/25/03 (DAK) new
 */
public interface AtomPairActive extends AtomSetActive {
	
	/**
	 * Performs on the given pair the action defined by this interface.
	 * @param pair atom pair subject to the action
	 */
	public void actionPerformed(AtomPair pair);
	
	
	/**
	 * Provides instance of a AtomPairActive.Wrapper class devoted to this
	 * action.  It is expected that the same instance is returned with every
	 * call.
	 * @return Wrapper
	 * @see etomica.AtomPairActive.InnerWrapper
	 */
	public InnerWrapper innerWrapper();
	
	
	/**
	 * Provides instance of a AtomPairActive.OuterWrapper class devoted to this
	 * action.  It is expected that the same instance is returned with every
	 * call.
	 * @return InnerWrapper
	 * @see etomica.AtomPairActive.OuterWrapper
	 */
	public OuterWrapper outerWrapper();
	
	/**
	 * Wrapper class that makes an AtomPairActive suitable for input to an
	 * AtomIterator. This is needed to enable an AtomPair iterator to perform an
	 * AtomPairActive on all its pairs. Class is AtomActive and contains the
	 * desired AtomPairActive. The first Atom is set externally  before being
	 * fed into an AtomIterator, which iterates values for the second Atom of
	 * the pair; the pair is sent to the wrapped AtomPairActive on each
	 * iteration.
	 */
		public static final class InnerWrapper implements AtomActive {
			public final AtomPair pair;
			private final AtomPairActive action;
			public InnerWrapper(AtomPairActive action, AtomPair pair) {
				this.action = action;
				this.pair = pair;
			}
			public void actionPerformed(AtomSet atomSet) {
				actionPerformed((Atom)atomSet);
			}
			public void actionPerformed(Atom a) {
				pair.atom2 = a;
				pair.reset();
				action.actionPerformed(pair);
			}
			public void actionPerformed() {
				pair.reset();
				action.actionPerformed(pair);
			}
		}//end of InnerWrapper
    
		public static final class OuterWrapper implements AtomActive {
			public final IteratorDirective directive = new IteratorDirective(IteratorDirective.UP);
			public Atom innerBasis;
			public AtomIterator aiInner;
			private final AtomPairActive action;
			private final InnerWrapper wrapper;
			public OuterWrapper(AtomPairActive action) {
				this.action = action;
				wrapper = action.innerWrapper();
			}
			public void actionPerformed(AtomSet atomSet) {
				actionPerformed((Atom)atomSet);
			}
			public void actionPerformed(Atom a) {
				wrapper.pair.atom1 = a;
				aiInner.all(innerBasis, directive.set(a), wrapper);
			}
		}//end of OuterWrapper
            
    

}
