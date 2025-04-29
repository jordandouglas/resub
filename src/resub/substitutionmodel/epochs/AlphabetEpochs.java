package resub.substitutionmodel.epochs;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Loggable;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;


@Description("Several epochs amino acid alphabets")
public class AlphabetEpochs extends CalculationNode implements Loggable, Function {
	
	
	// Epoch inputs
	final public Input<Integer> nepochsInput = new Input<>("nepochs", "Number of epochs, fixed throughout the analysis", Input.Validate.REQUIRED);
	final public Input<RealParameter> teInput = new Input<>("te", "The transition time of the epoch", Input.Validate.REQUIRED);
	final public Input<IntegerParameter> alphaInput = new Input<>("alpha", "State index 1 to consider for merging", Input.Validate.REQUIRED);
	final public Input<IntegerParameter> betaInput = new Input<>("beta", "State index 2 to consider for merging", Input.Validate.REQUIRED);
	final public Input<IntegerParameter> indicatorInput = new Input<>("indicator", "Model indicator, 0=no merging, 1=refinement, 2=expand left, 3=expand right", Input.Validate.REQUIRED);
	final public Input<Function> upperInitInput = new Input<>("upperInit", "Upper limit of initial values for te", Input.Validate.OPTIONAL);
	
	final public Input<List<AlphaBetaConstraint>> constraintsInput = new Input<>("constraint", "ensures that at least one epoch contains these two states", new ArrayList<>());
	
	
	
	// Substitution model (optional)
	final public Input<GeneralSubstitutionModel> substModelInput = new Input<>("substModel", "An amino acid substitution model that applies to the youngest epoch");
	final public Input<Frequencies> frequenciesInput = new Input<>("frequencies", "substitution model equilibrium state frequencies in the youngest epoch");
	//final public Input<Boolean> complexEigenInput = new Input<>("complex", "whether to use complex number eigen decomposition on reduced matrices (potentially more numerically stable)", false);
	
	final public Input<RealParameter> prefineInput = new Input<>("prefine", "for state initialisation only");
	final public Input<RealParameter> pactiveInput = new Input<>("pactive", "for state initialisation only");
	
	//final public Input<DataType> datatypeInput = new Input<>("datatype", "type of data", new Aminoacid() );
	
	
	boolean needsUpdate;
	int nepochs;
	int nOldEpochs;
	final DataType datatype = new Aminoacid();
	final static String states = "ACDEFGHIKLMNPQRSTVWY";
	List<String> stateList;
	List<AlphabetEpoch> epochs = new ArrayList<>();
	List<AlphaBetaConstraint> constraints;
	
	@Override
	public void initAndValidate() {
		
		
		
		this.nepochs = nepochsInput.get();
		this.nOldEpochs = this.nepochs - 1;
		
		int nstates = datatype.getStateCount();
		
		if (this.nepochs < 2 || this.nepochs > nstates - 1) {
			throw new IllegalArgumentException("nepochs must be between 2 and the number of states-1");
		}
		
		// List of amino acids
		stateList = Arrays.asList(states.split(""));
		this.constraints = constraintsInput.get();
		
		// Validate substitution model and frequencies, if provided
		if (substModelInput.get() != null) {
			if (frequenciesInput.get() == null) {
				throw new IllegalArgumentException("Please provide amino acid 'frequencies' to acocmpany the substModel");
			}
			
			if (!substModelInput.get().canHandleDataType(datatype)) {
				throw new IllegalArgumentException("Please ensure the 'substModel' supports amino acids");
			}
			
			if (frequenciesInput.get().getFreqs().length != stateList.size()) {
				throw new IllegalArgumentException("Please ensure there are " + stateList.size() + " elements in 'frequencues'");
			}
			
		}
		

		// Set dimensions
		boolean needsInit = teInput.get().getDimension() != this.nOldEpochs
							|| alphaInput.get().getDimension() != this.nOldEpochs
							|| betaInput.get().getDimension() != this.nOldEpochs
							|| indicatorInput.get().getDimension() != this.nOldEpochs
							|| teInput.get().getValue() <= 0 
							|| (upperInitInput.get() != null && teInput.get().getValue() >= upperInitInput.get().getArrayValue()
							|| alphaInput.get().getValue() == betaInput.get().getValue());
		
		teInput.get().setDimension(this.nOldEpochs);
		alphaInput.get().setDimension(this.nOldEpochs);
		betaInput.get().setDimension(this.nOldEpochs);
		indicatorInput.get().setDimension(this.nOldEpochs);
		
		
		
		// Set mins and maxes
		
		teInput.get().setLower(0.0);
		alphaInput.get().setLower(0);
		betaInput.get().setLower(0);
		indicatorInput.get().setLower(0);
		indicatorInput.get().setUpper(2);
		
		
		
		// The upper limit is variable but this is the maximum possible limit
		alphaInput.get().setUpper(nstates-1); 
		betaInput.get().setUpper(nstates-1); 
		
		
		//Log.warning("Getting upper " + indicatorInput.get().getUpper() + " " + this.nOldEpochs);
		
		
		// Initialise values
		if (needsInit){
			
			
			// Ensure that no 2 epochs have the same beta
			List<Integer> availableBetas = new ArrayList<>();
			for (int i = 0; i < this.stateList.size(); i ++) availableBetas.add(i);
			
			
			double pactive = pactiveInput.get() == null ? 0.5 : pactiveInput.get().getArrayValue();
			double prefine = prefineInput.get() == null ? 0.5 : prefineInput.get().getArrayValue();
			
			double dt = upperInitInput.get() == null ? 1.0 : upperInitInput.get().getArrayValue()/(this.nOldEpochs+1.0);
			int constraintNumber = 0;
			
			for(int i = 0; i < this.nOldEpochs; i ++) {
				
				// Stagger the epoch ages - going backwards in time
				teInput.get().setValue(i, (i+1)*dt);
				
				// Random resub state to start with
				//indicatorInput.get().setValue(i, Randomizer.nextInt(3)); 
				
				// Initialise randomly based on model parameters
				int indicator = Randomizer.nextFloat() >= pactive ? 0 : Randomizer.nextFloat() <= prefine ? 1 : 2;
				Log.warning("setting to " + indicator);
				indicatorInput.get().setValue(i, indicator);
				
				// Ensure this is a valid state
				if (constraintNumber < this.constraints.size()) {
					AlphaBetaConstraint constraint = this.constraints.get(constraintNumber);
					alphaInput.get().setValue(i, constraint.getState1());
					betaInput.get().setValue(i, constraint.getState2());
					availableBetas.remove(constraint.getState2());
					constraintNumber++;
				}else {
					
					// Set alpha and beta to random states such that alpha != beta
					int k = Randomizer.nextInt(availableBetas.size());
					int beta = availableBetas.remove(k);
					betaInput.get().setValue(i, beta); // Set each beta to a different index
					
					// Pick an alpha from the remaining ones
					k = Randomizer.nextInt(availableBetas.size());
					alphaInput.get().setValue(i, availableBetas.get(k)); 
				}
				
			}
		}
		
		// Create the epochs
		epochs.clear();
		for(int i = 0; i < this.nOldEpochs; i ++) {
			AlphabetEpoch epoch = new AlphabetEpoch(i, teInput.get(), alphaInput.get(), betaInput.get(), indicatorInput.get(), stateList.size(), false);
			epochs.add(epoch);
		}
		
		
		// Update the epochs
		this.needsUpdate = true;
		update();
		
	}
	

	
	public boolean update() {
		
		
		if (!needsUpdate) return true; 
		
		
		// Sort epochs by age, such that the inactive epochs are at the end
		Collections.sort(epochs);
		
		// Present day alphabet
		List<String> stateListLarge = new ArrayList<>(stateList);
		double[] ratesCurrentEpoch = getSubstRatesPresent();
		double[] freqsCurrentEpoch = getFrequenciesPresent();
		
		
		//Log.warning("____");
		
		
		// Set the alphabet of each epoch
		for(int i = 0; i < this.epochs.size(); i ++) {
			AlphabetEpoch epoch = this.epochs.get(i);
			
			
			
			if (epoch.getAlpha() == epoch.getBeta()) return false;
			
			//Log.warning("epoch " + i + " alpha=" + epoch.getAlpha() + "  beta=" + epoch.getBeta() + " I=" + epoch.getIndicator());
			
			
			double pStay = 1;
			List<String> stateListSmall = new ArrayList<>(stateListLarge);
			if (epoch.isActive()) {
				
				
				String alphaState = stateListLarge.get(epoch.getAlpha());
				String betaState = stateListLarge.get(epoch.getBeta());
				
				
				// Invalid state
				if (alphaState.equals(AlphabetEpoch.BLANK_STATE) || betaState.equals(AlphabetEpoch.BLANK_STATE)) {
					//System.out.println("invalid because of " + alphaState + " " + betaState + " epoch rank " + i);
					return false;
				}
				
				if (ratesCurrentEpoch != null && freqsCurrentEpoch != null) {
				
					// The proportion of states that remain on alpha after the transition
					pStay = freqsCurrentEpoch[epoch.getAlpha()] / (freqsCurrentEpoch[epoch.getAlpha()] + freqsCurrentEpoch[epoch.getBeta()]);
					
					
					// Refinement: replace both states with their merged state at position alpha
					if (epoch.isRefinement()) {
						
						//Log.warning("refinement");
						
						// Update subst model
						ratesCurrentEpoch = mergeRates(epoch.getAlpha(), epoch.getBeta(), ratesCurrentEpoch, freqsCurrentEpoch, stateListSmall.size());
						
						// Update frequencies
						freqsCurrentEpoch = mergeFrequencies(epoch.getAlpha(), epoch.getBeta(), freqsCurrentEpoch);
						
						// Merge the two state names
						String mergedStateName = alphaState + "/" + betaState;
						stateListSmall.set(stateListSmall.indexOf(alphaState), mergedStateName);
						stateListSmall.set(stateListSmall.indexOf(betaState), AlphabetEpoch.BLANK_STATE);
						
						
					}
					
					
					// Expansion: delete beta and keep alpha
					else  {
							
						//Log.warning("expansion");
						
						// Update subst model
						ratesCurrentEpoch = removeState(epoch.getBeta(), ratesCurrentEpoch, stateListSmall.size());
						
						// Update frequencies
						freqsCurrentEpoch = removeFrequency(epoch.getBeta(), freqsCurrentEpoch);
						//freqsCurrentEpoch = mergeFrequencies(epoch.getAlpha(), epoch.getBeta(), freqsCurrentEpoch);
						
						// Remove a state from the state list
						stateListSmall.set(stateListSmall.indexOf(betaState), AlphabetEpoch.BLANK_STATE);
							
					}
					
				
				}



			}
			
			
			
			epoch.setRank(i);
			epoch.setStatesAtEpoch(new ArrayList<>(stateListSmall), new ArrayList<>(stateListLarge));
			epoch.prepareSubstitionModel(ratesCurrentEpoch, freqsCurrentEpoch, pStay);
			stateListLarge = stateListSmall;
			
			//System.out.println(epoch);
			//if (freqsCurrentEpoch != null) {
				//System.out.println(epoch.getFrequencyString());
				//System.out.println(epoch.getRateString());
			//}
			
			
			
		}
		
		
		needsUpdate = false;
		return true;
		
	}
	
	

	/**
	 * Return amino acid subst rates in the youngest epoch, if available
	 * @return
	 */
	public double[] getSubstRatesPresent() {
		if (this.substModelInput.get() == null) return null;
		this.substModelInput.get().setupRelativeRates();
		double[] rates = substModelInput.get().getRelativeRates();
//		Log.warning("r A<=>C:" + rates[0]);
//		Log.warning("r A<=>R:" + rates[13]);
//		Log.warning("r W<=>Y:" + rates[rates.length-1]);
		return rates;
	}

	
	/**
	 * Return amino acid frequencies in the youngest epoch, if available
	 * @return
	 */
	public double[] getFrequenciesPresent() {
		if (this.frequenciesInput.get() == null) return null;
		return this.frequenciesInput.get().getFreqs();
	}

	public List<AlphabetEpoch> getSortedEpochs() {
		if (!update()) return null;
		return epochs;
	}
	
	
	public int getNActiveEpochs() {
		
		update();
		int count = 0;
		for(int i = 0; i < this.epochs.size(); i ++) {
			if (this.epochs.get(i).isActive()) count ++;
		}
		return count;
	}
	
	public DataType getDataType() {
		return this.datatype;
	}
	
	public RealParameter getTransitionAges() {
		return this.teInput.get();
	}
	
	public IntegerParameter getAlphas() {
		return this.alphaInput.get();
	}
	
	public IntegerParameter getBetas() {
		return this.betaInput.get();
	}
	
	public IntegerParameter getIndicators() {
		return this.indicatorInput.get();
	}

	public int getNOldEpochs() {
		return this.nOldEpochs;
	}

	@Override
	public void init(PrintStream out) {
		out.print(this.getID() + ".nactive\t");
		for(int i = 0; i < this.epochs.size(); i ++) {
			out.print("transition." + i + "\t");
			out.print("age." + i + "\t");
		}
	}

	
	@Override
	public void log(long sample, PrintStream out) {
		
		
		// Update
		update();
		
		// Print epochs by index number, not sorted. Then they can be compared with the logged times
		out.print(this.getNActiveEpochs() + "\t");
		for(int i = 0; i < this.epochs.size(); i ++) {
			AlphabetEpoch epoch = this.getEpochByRank(i);
			out.print(epoch.getTransitionString() + "\t");
			out.print((epoch.isActive() ? epoch.getAge() : 0) + "\t"); // Print 0 if inactive
		}
		
	}
	
	@Override
    public void store() {
		for(int i = 0; i < this.epochs.size(); i ++) {
			this.getEpochSubstitionModelByRank(i).store();
		}
		super.store();
	}
	
	
	
	@Override
	public void restore() {
		this.needsUpdate = true;
		for(int i = 0; i < this.epochs.size(); i ++) {
			this.getEpochSubstitionModelByRank(i).restore();
		}
		super.restore();
	}

	

	@Override
	public void close(PrintStream out) {
		
	}

	
	public AlphabetEpoch getEpochByIndex(int epochIndex) {
		this.update();
		for(int i = 0; i < this.epochs.size(); i ++) {
			if (this.epochs.get(i).getIndex() == epochIndex) return this.epochs.get(i);
		}
		return null;
	}
	
	
	public AlphabetEpoch getEpochByRank(int rank) {
		this.update();
		return this.epochs.get(rank);
	}
	
	public GeneralSubstitutionModel getEpochSubstitionModelByRank(int rank) {
		this.update();
		return this.epochs.get(rank).getSubstitionModel();
	}

	public double[] getEpochTransportMatrixByRank(int rank) {
		this.update();
		return this.epochs.get(rank).getTransportMatrix();
	}
	

	/**
	 * Normalise vector of frequencies so that they sum to 1
	 * @param freqs
	 */
	private void normaliseFrequencies(double[] freqs) {
		
		if (freqs == null) return;
		
		double sum = 0;
		for (int i = 0; i < freqs.length; i ++) {
			sum += freqs[i];
		}
		for (int i = 0; i < freqs.length; i ++) {
			freqs[i] = freqs[i] / sum;
		}
		
	}
	
	
	/**
	 * Remove the index, renormalise, and return a shorter array
	 * @param indexToRemove
	 * @param freqs
	 * @return
	 */
	private double[] removeFrequency(int indexToRemove, double[] freqs) {
		
		if (freqs == null) return null;
		
		double[] newFreqs = new double[freqs.length];
		for (int i = 0; i < freqs.length; i ++) {
			if (i == indexToRemove) {
				newFreqs[i] = 0;
			}else {
				newFreqs[i] = freqs[i];
			}
			
		}
		
		normaliseFrequencies(newFreqs); // Essential step
		return newFreqs;
		
	}

	
	/**
	 * Combine the two at position alpha
	 * @param alpha
	 * @param beta
	 * @param freqs
	 * @return
	 */
	private double[] mergeFrequencies(int alpha, int beta, double[] freqs) {
		
		if (freqs == null) return null;
		
		double[] newFreqs = new double[freqs.length];
		for (int i = 0; i < freqs.length; i ++) {
			if (i == beta) {
				newFreqs[i] = 0;
			}else if (i == alpha) {
				newFreqs[i] = freqs[alpha] + freqs[beta];
			}else {
				newFreqs[i] = freqs[i];
			}
			
		}
		
		normaliseFrequencies(newFreqs); // Not necessary but just being safe in case of numerical issues
		return newFreqs;
	}
	
	

	/**
	 * Remove this state from the matrix
	 * @param indexToRemove
	 * @param relativeRates
	 * @return
	 */
	private double[] removeState(int indexToRemove, double[] relativeRates, int nstates) {
		
		if (relativeRates == null) return null;
		
		// Copy everything over
		double[] reducedMatrix = new double[nstates * (nstates-1)];
		for (int i = 0; i < relativeRates.length; i ++) reducedMatrix[i] = relativeRates[i];
				
		for (int state = 0; state < nstates; state ++) {
			this.setRate(reducedMatrix, state, indexToRemove, 0);
			this.setRate(reducedMatrix, indexToRemove, state, 0);
		}
		
		
		return reducedMatrix;
	}

	
	/**
	 * Merge the two rates - requires frequencies
	 * @param alpha
	 * @param beta
	 * @param relativeRates
	 * @return
	 */
	private double[] mergeRates(int alpha, int beta, double[] relativeRates, double[] freqs, int nstates) {
		
		if (relativeRates == null) return null;
		
		// Copy everything over
		double[] reducedMatrix = new double[nstates * (nstates-1)];
		for (int i = 0; i < relativeRates.length; i ++) reducedMatrix[i] = relativeRates[i];
		
		
		// Merge
		double f1 = freqs[alpha];
		double f2 = freqs[beta];
		for (int state = 0; state < nstates; state ++) {
			
			if (state == alpha || state == beta) continue;
			
			
			double r1 = (this.getRate(relativeRates, state, alpha) + this.getRate(relativeRates, alpha, state)) / 2;
			double r2 = (this.getRate(relativeRates, state, beta) + this.getRate(relativeRates, beta, state)) / 2;
			double r = (r1*f1 + r2*f2) / (f1 + f2);
			this.setRate(reducedMatrix, state, alpha, r);
			this.setRate(reducedMatrix, alpha, state, r);
			
			
//			if (r1 <= 0 || r2 <= 0 || f1 <= 0 || f2 <= 0) {
//				Log.warning("state=" + state);
//				Log.warning("alpha=" + alpha);
//				Log.warning("beta=" + beta);
//				Log.warning("r1=" + this.getRate(relativeRates, state, alpha) + " = " + this.getRate(relativeRates, alpha, state));
//				Log.warning("r2=" + this.getRate(relativeRates, state, beta) + " = " + this.getRate(relativeRates, beta, state));
//				Log.warning("f1=" + f1);
//				Log.warning("f2=" + f2);
//				Log.warning("r=" + r);
//				Log.warning("-----------------");
//			}
			

			
		}
		
		
		for (int state = 0; state < nstates; state ++) {
			this.setRate(reducedMatrix, state, beta, 0);
			this.setRate(reducedMatrix, beta, state, 0);
		}
		
		return reducedMatrix;
		
	}
	
	
	private void setRate(double[] rates, int to, int from, double value) {
		
		if (from == to) return;
		
		int k = 0;
		for (int i = 0; i < this.getNStatesTotal(); i ++) {
			for (int j = 0; j < this.getNStatesTotal(); j ++) {
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
		for (int i = 0; i < this.getNStatesTotal(); i ++) {
			for (int j = 0; j < this.getNStatesTotal(); j ++) {
				if (i == j) continue;
				if (i == to && j == from) {
					return rates[k];
				}
				k ++;
			}
		}
		
		
		Log.warning("Dev error 26297 cannot find rate at i=" + to + " j=" + from);
		return -2;
	}
	


	public GeneralSubstitutionModel getSubstModel() {
		return substModelInput.get();
	}

	
	@Override
	protected boolean requiresRecalculation() {
		
		if (InputUtil.isDirty(alphaInput)) needsUpdate = true;
		if (InputUtil.isDirty(betaInput)) needsUpdate = true;
		if (InputUtil.isDirty(indicatorInput)) needsUpdate = true;
		if (InputUtil.isDirty(teInput)) needsUpdate = true;
		
		if (InputUtil.isDirty(substModelInput)) needsUpdate = true;
		if (InputUtil.isDirty(frequenciesInput)) needsUpdate = true;
		
		return needsUpdate;
		
	}

	public int getNStatesTotal() {
		return stateList.size();
	}
	
	public List<String> getStateList() {
		return new ArrayList<>(stateList); // Clone
	}

	public double[] getRootFrequencies() {
		this.update();
		
		int nepochs = this.getNActiveEpochs();
		if (nepochs == 0) {
			return frequenciesInput.get().getFreqs();
		}
		return this.getEpochByRank(nepochs-1).getFrequencies();
	}

	
	@Override
	public int getDimension() {
		return 1; // Only return 1 number: the youngest epoch age
	}

	@Override
	public double getArrayValue(int dim) {
		this.update();
		return this.getEpochByRank(dim).getEffectiveAge();
	}

	public void requestUpdate() {
		this.needsUpdate = true;
	}
	
	
	/**
	 * Are the alpha/beta constraints satisfied?
	 * @return
	 */
	public boolean constraintsAreSatisfied(){
		
		
		for (AlphaBetaConstraint constraint : this.constraints) {
			if (!constraint.isSatisfied(this)) return false;
		}
		
		return true;
		
	}
	
	

}



