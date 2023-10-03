package etomica.chem.elements;

    public class Nitrogen_1 extends ElementChemical {

        protected Nitrogen_1(String symbol) {
            this(symbol, 14.01);
        }

        protected Nitrogen_1(String symbol, double mass) {
            super(symbol, mass, 7);
        }

        public static final Nitrogen_1 INSTANCE = new Nitrogen_1("N_1");

        private static final long serialVersionUID = 1L;
}
