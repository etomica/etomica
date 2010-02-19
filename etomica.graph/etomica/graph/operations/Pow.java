package etomica.graph.operations;

import java.util.Set;

import etomica.graph.model.Graph;

public class Pow implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    assert (params instanceof PowParameters);
    Unary pcopy = new PCopy();
    Binary multiply = new Mul();
    Set<Graph> result = pcopy.apply(argument, params);
    for (byte pow = 2; pow < ((PowParameters) params).exponent(); pow++) {
      result = multiply.apply(result, argument, params);
    }
    return result;
  }
}