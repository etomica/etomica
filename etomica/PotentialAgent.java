package etomica;

/**
 * @author David Kofke
 *
 * Agent of a potential instance that is used to perform 
 * potential calculations.
 */
public abstract class PotentialAgent {
	
	public PotentialAgent(Potential potential) {
		this.potential = potential;
	}
        
	public abstract void calculate(IteratorDirective id, PotentialCalculation pc);
    
	//Sets the basis for iteration
	public abstract Potential set(Atom[] atoms);
	public abstract Potential set(SpeciesMaster s);    

	public final Potential potential;
	
	/**
	 * Linked list of PotentialAgent instances.
	 */
	public static final class List {
    	private final Linker header = new Linker();
    	public void add(PotentialAgent agent) {
    		addFirst(agent);  
    	}
		public synchronized void addLast(PotentialAgent agent) {
			Linker linker = new Linker(agent, header, header.previous);//put at end of list
			header.previous.next = linker;
			header.previous = linker;    		    
		}
		public synchronized void addFirst(PotentialAgent agent) {
			Linker linker = new Linker(agent, header.next, header);//put at end of list
			header.next.previous = linker;
			header.next = linker;    		    
		}
    	public synchronized void remove(PotentialAgent agent) {
    		for(Linker linker=header.next; linker!=header; linker=linker.next) {
    			if(linker.agent == agent) {
    				remove(linker);
    				break;
    			}
    		}
    	}
    	private void remove(Linker linker) {
    		linker.previous = linker.next;
    		linker.next = linker.previous;
    	}
    	public synchronized void clear() {
    		header.next = header.previous = header;
    	}
    	public Iterator iterator() {return new Iterator(this.header);}
	}//end of List
	
	public static class Iterator {
		private final Linker header;
		private Linker current;
		private Iterator(Linker header) {
			this.header = header;
			current = header;
		}
		public boolean hasNext() {return current != header;}
		public void reset() {current = header;}
		public PotentialAgent next() {
			current = current.next;
			return current.previous.agent;
		}		
	}//end of Iterator
	
	private static class Linker {
		private Linker next, previous;
		private final PotentialAgent agent;
		private Linker() {
			this.agent = null;
			this.next = this;
			this.previous = this;
		}
		private Linker(PotentialAgent agent, Linker next, Linker previous) {
			this.agent = agent;
			this.next = next;
			this.previous = previous;
		}
	}//end of Linker
    


}//end of Potential
