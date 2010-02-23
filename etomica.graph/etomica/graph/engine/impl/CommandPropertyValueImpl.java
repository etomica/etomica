package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.CommandPropertyValue;
import etomica.graph.engine.Parser.Property;
import etomica.graph.engine.Parser.Value;

public class CommandPropertyValueImpl extends CommandImpl implements CommandPropertyValue {

  private Property property;
  private Value value;

  public CommandPropertyValueImpl(String command, Property property, Value value) {

    super(command);
    this.property = property;
    this.value = value;
  }

  public Property getProperty() {

    return this.property;
  }

  public Value getValue() {

    return this.value;
  }
}
