package etomica.utility;

/**
 * Defines a function as chain of other functions.
 * For example, defines F(x) = f3(f2(f1(x))).
 * Defaults to identity function if no subfunctions are specified.
 */
public final class FunctionChain implements Function {
    
    private Link head, tail;
    
    /**
     * Calculates function by applying subfunctions, in the ordered they were added.
     */
    public double f(double x) {
        if(head == null) return x;
        for(Link link=head; link!=null; link=link.next) {
            x = link.function.f(x);
        }
        return x;
    }
    
    /**
     * Inverts function by calling inverse methods of subfunctions, in reverse order of their addition.
     */
    public double inverse(double f) {
        if(tail == null) return f;
        for(Link link=tail; link!=null; link=link.previous) {
            f = link.function.inverse(f);
        }
        return f;
    }
    
    /**
     * Evaluates derivative via the chain rule.
     */
    public double dfdx(double x) {
        if(head == null) return 1.0;
        double df = 1.0;
        for(Link link=head; link!=null; link=link.next) {
            df = link.function.dfdx(x);
            x = link.function.f(x);
        }
        return df;
    }
    
    /**
     * Adds the given function to the chain of functions.
     * Functions are applied sequentially in the order they are added.
     * First function added is the first applied, and so on.
     */
    public void add(Function newF) {
        tail = new Link(newF, tail);//new function goes to tail of list
        if(head==null) head = tail;//make head also if first added to list
    }
    
    /**
     * Removes the given function from the chain.
     * If the function is not in the chain, returns with no error or warning message.
     * If function is present multiple times in chain, remove the first one added.
     */
    public void remove(Function f) {
        for(Link link=head; link!=null; link=link.next) {
            if(link.function == f) {//found it
                if(link.previous != null) link.previous.next = link.next;
                else head = link.next;
                
                if(link.next != null) link.next.previous = link.previous;
                else tail = link.previous;
                
                break;
            }//end if
        }//end for
    }
    
    /**
     * Class for making linked list of Functions
     */
    public static class Link {
        public Link next, previous;
        public final Function function;
        public Link(Function func, Link prev) {
            function = func;
            previous = prev;
            if(prev != null) prev.next = this;
        }
    }
}