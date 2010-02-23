package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.CommandValue;
import etomica.graph.engine.Parser.Value;

public class CommandValueImpl extends CommandImpl implements CommandValue {

  private Value value;

  public CommandValueImpl(String command, Value value) {

    super(command);
    this.value = value;
  }

  public Value getValue() {

    return value;
  }
}