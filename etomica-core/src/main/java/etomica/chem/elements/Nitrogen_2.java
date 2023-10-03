package etomica.chem.elements;

    public class Nitrogen_2 extends ElementChemical {

        protected Nitrogen_2(String symbol) {
            this(symbol, 14.01);
        }

        protected Nitrogen_2(String symbol, double mass) {
            super(symbol, mass, 7);
        }

        public static final Nitrogen_2 INSTANCE = new Nitrogen_2("N_2");

        private static final long serialVersionUID = 1L;
}
