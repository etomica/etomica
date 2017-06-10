/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.function;

/**
 * Defines a function as chain of other functions.
 * For example, defines F(x) = f3(f2(f1(x))).
 * Defaults to identity function if no subfunctions are specified.
 */
public final class FunctionChain implements Function, java.io.Serializable {
    
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
            f = ((FunctionInvertible)link.function).inverse(f);
        }
        return f;
    }
    
    /**
     * Evaluates derivative via the chain rule.
     */
    public double df(int n, double x) {
        if(head == null) return 1.0;
        if(n < 0) {
            throw new IllegalArgumentException("order of derivative must be non-negative");
        }
        if(n == 0) {
            return f(x);
        }
        if(n > 1) {
            throw new IllegalArgumentException("Only first-order derivative available for FunctionChain");
        }
        double df = 1.0;
        for(Link link=head; link!=null; link=link.next) {
            df = ((FunctionDifferentiable)link.function).df(1, x);
            x = link.function.f(x);
        }
        return df;
    }
    
    /**
     * Adds the given function to the chain of functions.
     * Functions are applied sequentially in the order they are added.
     * First function added is the first applied, and so on.
     */
    public void add(IFunction newF) {
        tail = new Link(newF, tail);//new function goes to tail of list
        if(head==null) head = tail;//make head also if first added to list
    }
    
    /**
     * Removes the given function from the chain.
     * If the function is not in the chain, returns with no error or warning message.
     * If function is present multiple times in chain, remove the first one added.
     */
    public void remove(IFunction f) {
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
    public static class Link implements java.io.Serializable {
        public Link next, previous;
        public final IFunction function;
        public Link(IFunction func, Link prev) {
            function = func;
            previous = prev;
            if(prev != null) prev.next = this;
        }
    }
}
