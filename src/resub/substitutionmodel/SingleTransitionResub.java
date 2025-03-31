package resub.substitutionmodel;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import resub.math.ResubMathUtil;

@Description("A special case of resub where there are just 2 epochs")
public class SingleTransitionResub extends SubstitutionModel.Base implements Loggable {
	
	final private static boolean EXPAND_FREQS_RENORM = true;
	final private String AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY";
	
	
	final public Input<GeneralSubstitutionModel> substModelInput = new Input<>("substModel", ""
			+ "substitution model for the branches in the tree ", null, Validate.REQUIRED);
	
	final public Input<BooleanParameter> refineInput = new Input<>("refine", "a vector of booleans (one per aaRS coalescent node)."
			+ "				Suppose there is an epoch node with two child states: A and B, where A and B are sets of states."
			+ "				If refine=0, then the state of this node will be A union B."
			+ "				If refine=1, this state will be either A or B.", Input.Validate.REQUIRED);
	
	final public Input<BooleanParameter> expandInput = new Input<>("expand", "a vector of booleans (one per aaRS coalescent node)."
			+ "				Suppose there is an epoch node with two child states: A and B, where A and B are sets of states."
			+ "				If refine=0 and expand=0, this state will be A."
			+ "				If refine=0 and expand=1, this state will be B."
			+ "				If refine=1, this value is ignored", Input.Validate.REQUIRED);
	
	final public Input<BooleanParameter> useResubInput = new Input<>("useResub", "should resub even be used?", Input.Validate.REQUIRED);
	final public Input<RealParameter> transitionHeightInput = new Input<>("h1", "the age (height) when it changes between 2 epochs", Input.Validate.REQUIRED);
	
	final public Input<RealParameter> transitionFreqInput = new Input<>("pi", "probability going to state1 at transition (defaults to normalised equilibrium frequencies)", Input.Validate.OPTIONAL);
	
	
	final public Input<String> state1Input = new Input<>("state1", "one of the two states to merge", Input.Validate.REQUIRED);
	final public Input<String> state2Input = new Input<>("state2", "one of the two states to merge", Input.Validate.REQUIRED);
	
	
	boolean needsUpdate;
	GeneralSubstitutionModel substModel;
	GeneralSubstitutionModel substModelSmall;
	
	double[] tmpMatrix3;
	double[] transformerMatrix;
	
	
	Frequencies freqsTopInput;
	RealParameter freqsTopRealParameterInput;
	RealParameter ratesTopInput;
	
	int state1, state2;
	
	double[] rootFrequencies;
	
	
	int stateBeingKept;
	int stateBeingDropped;
	
	// Matrices used for various rate matrix calculations. They are initialised once to save time
	double[] tmpMatrix;
	double[] tmpMatrix2;
	
	double[][] rateMatrix;
	
	@Override
    public void initAndValidate() {
        
        this.needsUpdate = true;
        
        frequencies = frequenciesInput.get();
        
        // Initialise merge
        this.nrOfStates = this.frequencies.getFreqs().length;
        this.substModel = substModelInput.get();
        expandInput.get().setDimension(1);
        refineInput.get().setDimension(1);
        
        
        if (this.substModel.getStateCount() != nrOfStates) {
        	throw new IllegalArgumentException("Number of states in the substModel must match number of frequencies " + this.substModel.getStateCount() + " != " + nrOfStates);
        }
        
        
        // Initialise rates for smaller matrix
        List<Double> rates = new ArrayList<>();
        for (int i = 0; i < this.nrOfStates*(this.nrOfStates-1);  i++) {
        	rates.add(1.0);
        }

        
        ratesTopInput = new RealParameter();
        ratesTopInput.initByName("value", rates);
        
        
        // Initialise freqs
        List<Double> f = new ArrayList<>();
        for (int i = 0; i < this.nrOfStates;  i++) {
        	f.add(1.0 / this.nrOfStates);
        }
        freqsTopRealParameterInput = new RealParameter();
        freqsTopRealParameterInput.initByName("value", f);
        freqsTopInput = new Frequencies();
        freqsTopInput.initByName("frequencies", freqsTopRealParameterInput);
        
        // Init small subst model
        substModelSmall = new GeneralSubstitutionModel();
        substModelSmall.initByName("rates", ratesTopInput, "frequencies", freqsTopInput);
        

        
        this.tmpMatrix = new double[(this.nrOfStates+1)*(this.nrOfStates+1)];
        this.tmpMatrix2 = new double[(this.nrOfStates+1)*(this.nrOfStates+1)];
        this.rateMatrix = new double[this.nrOfStates][this.nrOfStates];
        
        
        this.tmpMatrix3 = new double[(this.nrOfStates+1)*(this.nrOfStates+1)];
        this.transformerMatrix = new double[(this.nrOfStates+1)*(this.nrOfStates+1)];
        
        // State numbers
        String[] aa = AMINO_ACIDS.split("");
        for (int i = 0; i < aa.length; i++) {
        	String a = aa[i];
        	if (a.equals(state1Input.get())) {
        		state1 = i;
        	}
        	if (a.equals(state2Input.get())) {
        		state2 = i;
        	}
        }
        
        if (state1 == state2) {
        	throw new IllegalArgumentException("Error: please ensure that state1 != state2");
        }

        

        
        update();

        
	}
	
	
	@Override
    protected boolean requiresRecalculation() {
		boolean update = modelIsFilthy();
		
		//Log.warning("filthy: " + update);
		
        if (update || modelIsDirty()) {
        	needsUpdate = true;
        }
        return update;
    }
	
	public boolean modelIsDirty() {
		if (InputUtil.isDirty(refineInput)) return true;
		if (InputUtil.isDirty(expandInput)) return true;
		if (InputUtil.isDirty(useResubInput)) return true;
        if (modelIsFilthy()) return true;
        return false;
	}
	
	

	
	public boolean modelIsFilthy() {
		if (InputUtil.isDirty(refineInput)) return true;
		if (InputUtil.isDirty(expandInput)) return true;
		if (InputUtil.isDirty(useResubInput)) return true;
		if (InputUtil.isDirty(transitionHeightInput)) return true;
        if (InputUtil.isDirty(substModelInput)) return true;
        if (InputUtil.isDirty(frequenciesInput)) return true;
        if (InputUtil.isDirty(transitionFreqInput)) return true;
        return false;
	}
	
	
	@Override
    public void store() {
		substModelSmall.store();
		super.store();
	}
	
	
	
	@Override
	public void restore() {
		this.needsUpdate = true;
		substModelSmall.restore();
		super.restore();
	}

	


	@Override
	public double[] getFrequencies() {
		if (!useResubInput.get().getValue()) {
			return frequencies.getFreqs();
		}
		update();
		return rootFrequencies;
	}
	


	@Override
	public void init(PrintStream out) {
		String model1 = "resub.no";
		String model2 = "resub" + state1Input.get() + state2Input.get();
		String model3 = "resub"  + state1Input.get();
		String model4 = "resub"  + state2Input.get();
		out.print(model1 + "\t" + model2 + "\t" + model3 + "\t" + model4 + "\t");
	}



	@Override
	public void log(long sample, PrintStream out) {
		
		
		boolean useResub = useResubInput.get().getValue();
		boolean isRefining = this.refineInput.get().getValue(0);
		boolean expandLeft = this.expandInput.get().getValue(0);
		
		
		String value1 = useResub ? "0" : "1";
		String value2 = "0", value3 = "0", value4 = "0";
		if (useResub && isRefining) {
			value2 = "1";
		}else if (useResub && !isRefining) {
			if (expandLeft) {
				value3 = "1";
			}else {
				value4 = "1";
			}
		}
		
		
		
		out.print(value1 + "\t" + value2 + "\t" + value3 + "\t" + value4 + "\t");
		

		
	}




	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

	
	
	
	private double[] convertRatesToSymmetric(double[] ratesIn) {
		
		if (ratesIn.length == this.nrOfStates * (this.nrOfStates-1)) {
			return ratesIn;
		}
		
		double[] ratesOut = new double[this.nrOfStates * (this.nrOfStates - 1)];
        int k = 0;
        for (int i = 0; i < this.nrOfStates; i++) {
            for (int j = i+1; j < this.nrOfStates; j++) {
            	double rateVal = ratesIn[k];
            	this.setRate(ratesOut, i, j, rateVal);
            	this.setRate(ratesOut, j, i, rateVal);
            	k++;
            }
        }
        
        return ratesOut;
		
		
	}
	
	
	/**
	 * Update the model. Most of the core logic is here
	 */
	public void update() {
		
		if (!needsUpdate) return;
		needsUpdate = false;
		
		// Get the ground level rate matrix (no frequencies)
        substModelInput.get().setupRelativeRates();	
        
        
        // Prepare rates and frequencies for reduced matrix
		boolean isRefining = this.refineInput.get().getValue(0);
		boolean expandLeft = this.expandInput.get().getValue(0);
		
		double[] ratesBottom = substModelInput.get().getRelativeRates();
		double[] freqsBottom = substModelInput.get().getFrequencies();
		
		
		// Ensure the rates have the right dimension
		ratesBottom = this.convertRatesToSymmetric(ratesBottom);
		
		
		double[] ratesTop = new double[this.nrOfStates * (this.nrOfStates-1)];
		double[] frequenciesTop = new double[this.nrOfStates];
		for (int i = 0; i < ratesBottom.length; i ++) ratesTop[i] = ratesBottom[i];
		for (int i = 0; i < frequenciesTop.length; i ++) frequenciesTop[i] = freqsBottom[i];
		
		if (isRefining) {
			
			//Log.warning("refining");
			
			double f1 = freqsBottom[state1];
			double f2 = freqsBottom[state2];
			for (int state = 0; state < this.nrOfStates; state ++) {
				double r1 = (this.getRate(ratesBottom, state, state1) + this.getRate(ratesBottom, state1, state)) / 2;
				double r2 = (this.getRate(ratesBottom, state, state2) + this.getRate(ratesBottom, state2, state)) / 2;
				double r = (r1*f1 + r2*f2) / (f1 + f2);
				this.setRate(ratesTop, state, state1, r);
				this.setRate(ratesTop, state1, state, r);
			}
			
			
			for (int state = 0; state < this.nrOfStates; state ++) {
				this.setRate(ratesTop, state, state2, 0);
				this.setRate(ratesTop, state2, state, 0);
			}
			
			
			frequenciesTop[state1] = f1 + f2;
			frequenciesTop[state2] = 0;
			stateBeingDropped = state2;
			stateBeingKept = state1;
			
			
		}else {
			
			
			//Log.warning("expanding");
			
			// Take the left child
			if (expandLeft) {
				for (int state = 0; state < this.nrOfStates; state ++) {
					this.setRate(ratesTop, state, state2, 0);
					this.setRate(ratesTop, state2, state, 0);
				}
				
				
				if (!EXPAND_FREQS_RENORM) {
					double f1 = freqsBottom[state1];
					double f2 = freqsBottom[state2];
					frequenciesTop[state1] = f1 + f2;
				}else {
					double f1 = freqsBottom[state1];
					frequenciesTop[state1] = f1;
				}
				
				frequenciesTop[state2] = 0;
				stateBeingDropped = state2;
				stateBeingKept = state1;
				
			}
			
			
			
			// Take the right child
			else {
				for (int state = 0; state < this.nrOfStates; state ++) {
					this.setRate(ratesTop, state, state1, 0);
					this.setRate(ratesTop, state1, state, 0);
				}
				
				if (!EXPAND_FREQS_RENORM) {
					double f1 = freqsBottom[state1];
					double f2 = freqsBottom[state2];
					frequenciesTop[state1] = f1 + f2;
				}else {
					double f2 = freqsBottom[state2];
					frequenciesTop[state2] = f2;
				}
				
				frequenciesTop[state1] = 0;
				stateBeingDropped = state1;
				stateBeingKept = state2;
				
			}
			
			
		}
		
		
		// Normalise frequencies
		double freqSum = 0;
		for (int i = 0; i < frequenciesTop.length; i ++) {
			freqSum += frequenciesTop[i];
		}
		for (int i = 0; i < frequenciesTop.length; i ++) {
			frequenciesTop[i] = frequenciesTop[i] / freqSum;
			freqsTopRealParameterInput.setValue(i, frequenciesTop[i]);
		}
		

		
		// Set top rates 
		for (int i = 0; i < ratesTop.length; i ++) {
			ratesTopInput.setValue(i, ratesTop[i]);
		}
		
		// Transformer matrix
		double p1, p2;
		if (transitionFreqInput.get() != null) {
			p1 = transitionFreqInput.get().getArrayValue();
			p2 = 1-p1;
		}else {
			p1 = freqsBottom[stateBeingKept];
			p2 = freqsBottom[stateBeingDropped];
			double psum = p1 + p2;
			p1 = p1 / psum;
			p2 = p2 / psum;
		}
		
		for (int from = 0; from < this.nrOfStates; from++) {
			for (int to = 0; to < this.nrOfStates; to++) {
				
				// Diagonal matrix, except for the row being expanded/refined, which has p1 and p2
				int index = this.getIndex(from, to);
				
				
				if (from != stateBeingKept && from == to) {
					this.transformerMatrix[index] = 1;
				}
				
				else if (from == stateBeingKept && to == stateBeingKept) {
					this.transformerMatrix[index] = p1;
				}
				
				else if (from == stateBeingKept && to == stateBeingDropped) {
					this.transformerMatrix[index] = p2;
				}
				
				else {
					this.transformerMatrix[index] = 0;
				}
			}
		}
		

		
		// Set top frequencies
		rootFrequencies = frequenciesTop;
		
		
		// Update subst model
		freqsTopInput.doUpdate();
        this.substModelSmall.doUpdate();
		
		
	}
	

	private void setRate(double[] rates, int to, int from, double value) {
		
		if (from == to) return;
		
		
		int k = 0;
		for (int i = 0; i < this.nrOfStates; i ++) {
			for (int j = 0; j < this.nrOfStates; j ++) {
				if (i == j) continue;
				if (i == to && j == from) {
					rates[k] = value;
					return;
				}
				k ++;
			}
		}
	}
	
	
	
	private double getRate(double[] rates, int to, int from) {
		
		if (from == to) return -1;
		
		
		int k = 0;
		for (int i = 0; i < this.nrOfStates; i ++) {
			for (int j = 0; j < this.nrOfStates; j ++) {
				if (i == j) continue;
				if (i == to && j == from) {
					return rates[k];
				}
				k ++;
			}
		}
		
		return -2;
	}
	
	

	
	@Override
	public void getTransitionProbabilities(Node node, double startHeight, double endHeight, double rate, double[] outMatrix) {
		
		
		
		
		double distance = (startHeight - endHeight) * rate;
		if (distance < 0) {
			
			throw new IllegalArgumentException("Error 1: negative branch length detected! " + distance);
		}
		
		double transitionHeight = this.transitionHeightInput.get().getValue();
		
		
		// Make sure the aux arrays are the right length
		if (tmpMatrix.length != outMatrix.length) {
			tmpMatrix = new double[outMatrix.length];
		}
		if (tmpMatrix2.length != outMatrix.length) {
			tmpMatrix2 = new double[outMatrix.length];
		}
		if (tmpMatrix3.length != outMatrix.length) {
			tmpMatrix3 = new double[outMatrix.length];
		}
		
		
	
		// Case 1: use vanilla subst model
		if (!useResubInput.get().getValue() || startHeight < transitionHeight) {
			
			//Log.warning("ZZZ");
			
			this.substModel.getTransitionProbabilities(node, startHeight, endHeight, rate, outMatrix);
			tidyAndValidateProbs(outMatrix, "full matrix 1", false);
			return;
		}

		
		update();
		
		
		//if (Math.abs(endHeight - transitionHeight) < 1e-8) {
			//Log.warning("XXX");
		//}
		
		// Case 2: old sequence -> use 19 character alphabet
		if (endHeight > transitionHeight) {
			//Log.warning("QQQQ" + endHeight + " - " + transitionHeight);
			this.substModelSmall.getTransitionProbabilities(node, startHeight, endHeight, rate, outMatrix);
			tidyAndValidateProbs(outMatrix, "reduced matrix 1", true);
			return;
		}
		
		// Case 3: branch crosses the transition -> multiply the big and small matrices together
		this.substModelSmall.getTransitionProbabilities(node, startHeight, transitionHeight, rate, this.tmpMatrix);
		this.substModel.getTransitionProbabilities(node, transitionHeight, endHeight, rate, this.tmpMatrix2);
		
		tidyAndValidateProbs(tmpMatrix, "reduced matrix 2", true);
		tidyAndValidateProbs(tmpMatrix2, "full matrix 2", false);
		tidyAndValidateProbs(transformerMatrix, "transformer", false);
		

		if (Math.abs(endHeight - transitionHeight) < 1e-16) {
			
			//Log.warning("XXX");
			
			// Psmall * transformer
			ResubMathUtil.multiplyMatrices(this.tmpMatrix, this.transformerMatrix, outMatrix, this.nrOfStates);
			tidyAndValidateProbs(outMatrix, "multiplied matrix 1 threshold", false);
			
		}else {
			
			//Log.warning("YYY");
			
			// Psmall * transformer
			ResubMathUtil.multiplyMatrices(this.tmpMatrix, this.transformerMatrix, this.tmpMatrix3, this.nrOfStates);
			tidyAndValidateProbs(tmpMatrix3, "multiplied matrix 1", false);
			
			// (Psmall * transformer) * Plarge
			ResubMathUtil.multiplyMatrices(this.tmpMatrix3, this.tmpMatrix2, outMatrix, this.nrOfStates);
			tidyAndValidateProbs(outMatrix, "multiplied matrix 2", false);
			
		}
			
		

		
	}
	
	
	
	
	private int getIndex(int from, int to) {
		return from*this.nrOfStates + to;
	}
	
	

	/**
	 * Ensure that each row of the transition probability matrix sums to 1
	 * @param matrix
	 * @param name
	 */
	private void tidyAndValidateProbs(double[] matrix, String name, boolean reduced) {
		
		//if (true) return;
		
		// Correct the probs in the row being dropped, just in case of numerical instabilities
		if (reduced) {
			int k = 0;
			for (int i = 0; i < this.nrOfStates; i ++) {
				for (int j = 0; j < this.nrOfStates; j ++) {
					if (i == stateBeingDropped || j == stateBeingDropped) {
						if (i == j) {
							matrix[k] = 1;
						}else {
							matrix[k] = 0;
						}
					}
					k++;
				}
				
			}
		}
		
		
		int k = 0;
		for (int i = 0; i < this.nrOfStates; i ++) {
			double rowsum = 0;
			for (int j = 0; j < this.nrOfStates; j ++) {
				double prob = matrix[k];
				rowsum += prob;
				k++;
			}
			
			if (Math.abs(1.0-rowsum) > 1e-4) {
				Log.warning(name + " sums to " + rowsum + " on row " + i + " refine: " + this.refineInput.get().getValue(0));
				throw new IllegalArgumentException();
			}
			
			
		}
		
		
	}
	

	
	/**
	 * Get the string of this state number
	 * @param stateNr
	 * @return
	 */
	public String getState(int stateNr) {
		if (stateNr < 0) return "X";
		//DataType dt = dataTypeInput.get();
		//String s = dt.encodingToString(new int[] { stateNr });
		
		
		String[] aa = AMINO_ACIDS.split("");
		return aa[stateNr];

		
//		
//		//String traitNr = traitsInput.get().getTaxonValues()[stateNr];
//		for (String stateName : treeInput.get().getTaxaNames()) {
//			if (traitsInput.get().getStringValue(stateName).equals(""+stateNr)) {
//				if (!stateName.equals(s)) {
//					Log.warning("State " + stateNr + " has mismatch between " + stateName + " and " + s);
//				}
//				return stateName;
//			}
//		}
//		return "X";
		
	}
	
	

	@Override
    public boolean canReturnComplexDiagonalization() {
		
		// Set to true so that beagle will let us do the matrix exponentiation work
        return true;
    }
	
	@Override
	public EigenDecomposition getEigenDecomposition(Node arg0) {
		return null;
	}

	@Override
	public boolean canHandleDataType(DataType dataType) {
		return dataType instanceof Aminoacid;
	}

}
