/*
 * TreeLikelihood.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package resub.likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beagle.Beagle;
import beagle.BeagleFactory;
import beagle.BeagleFlag;
import beagle.BeagleInfo;
import beagle.InstanceDetails;
import beagle.ResourceDetails;
import beast.base.core.Description;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;


/**
 * BeagleTreeLikelihoodModel - implements a Likelihood Function for sequences on a tree.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Marc Suchard
 * @version $Id$
 */

@Description("Uses Beagle library to calculate Tree likelihood")
public class EpochBeagleTreeLikelihood extends EpochTreeLikelihood {

    // This property is a comma-delimited list of resource numbers (0 == CPU) to
    // allocate each BEAGLE instance to. If less than the number of instances then
    // will wrap around.
    // note: to use a different device, say device 2, start beast with
    // java -Dbeagle.resource.order=2 beast.app.BeastMCMC
    private static final String RESOURCE_ORDER_PROPERTY = "beagle.resource.order";
    private static final String PREFERRED_FLAGS_PROPERTY = "beagle.preferred.flags";
    private static final String REQUIRED_FLAGS_PROPERTY = "beagle.required.flags";
    private static final String SCALING_PROPERTY = "beagle.scaling";
    private static final String RESCALE_FREQUENCY_PROPERTY = "beagle.rescale";
    // Which scheme to use if choice not specified (or 'default' is selected):
    private static final PartialsRescalingScheme DEFAULT_RESCALING_SCHEME = PartialsRescalingScheme.DYNAMIC;

    private static int instanceCount = 0;
    private static List<Integer> resourceOrder = null;
    private static List<Integer> preferredOrder = null;
    private static List<Integer> requiredOrder = null;
    private static List<String> scalingOrder = null;

    private static final int RESCALE_FREQUENCY = 10000;
    private static final int RESCALE_TIMES = 1;

    boolean m_bUseAmbiguities, m_bUseTipLikelihoods;
    int m_nStateCount;
    int m_nNodeCount;
    private double [] matrices;
    private int epochCount;
    private double substModelThrehold;

    
    private double [] currentCategoryRates;
//    private double [] storedCurrentCategoryRates;
    private double [] currentFreqs;
    private double [] currentCategoryWeights;

    private double [] currentThresholds;
    private double [] storedThresholds;

    private int invariantCategory = -1;

    @Override
    public void initAndValidate() {
        boolean forceJava = Boolean.valueOf(System.getProperty("java.only"));
        if (forceJava) {
        	return;
        }
        
    	models = modelsInput.get().toArray(new SubstitutionModel[] {});
    	epochDates = epochDatesInput.get();
		if (models.length != epochDates.getDimension() + 1) {
			throw new IllegalArgumentException("The number of epoch dates should be one less than the number of substitution models");
		}
		epochDatesArray = new double[epochDates.getDimension()];
        partialsMap = new int[treeInput.get().getNodeCount()];
        epochCount = models.length;
        
        currentThresholds = new double[epochDates.getDimension()];
        storedThresholds  = new double[epochDates.getDimension()];

        initialize();
    }

    private boolean initialize() {
        m_nNodeCount = treeInput.get().getNodeCount();
        m_bUseAmbiguities = m_useAmbiguities.get();
        m_bUseTipLikelihoods = m_useTipLikelihoods.get();
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
        	throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        // substitutionModel = m_siteModel.substModelInput.get();
        branchRateModel = branchRateModelInput.get();
        if (branchRateModel == null) {
        	branchRateModel = new StrictClockModel();
        }
        m_branchLengths = new double[m_nNodeCount];
        storedBranchLengths = new double[m_nNodeCount];

        m_nStateCount = dataInput.get().getMaxStateCount();
        patternCount = dataInput.get().getPatternCount();

        //System.err.println("Attempt to load BEAGLE TreeLikelihood");

        eigenCount = epochCount;//this.branchSubstitutionModel.getEigenCount();

        double[] categoryRates = m_siteModel.getCategoryRates(null);
        // check for invariant rates category
        if (m_siteModel.hasPropInvariantCategory) {
	        for (int i = 0; i < categoryRates.length; i++) {
	        	if (categoryRates[i] == 0) {
	        		proportionInvariant = m_siteModel.getRateForCategory(i, null);
	                int stateCount = dataInput.get().getMaxStateCount();
	                int patterns = dataInput.get().getPatternCount();
	                calcConstantPatternIndices(patterns, stateCount);
	                invariantCategory = i;
	                
	                double [] tmp = new double [categoryRates.length - 1];
	                for (int k = 0; k < invariantCategory; k++) {
	                	tmp[k] = categoryRates[k];
	                }
	                for (int k = invariantCategory + 1; k < categoryRates.length; k++) {
	                	tmp[k-1] = categoryRates[k];
	                }
	                categoryRates = tmp;
	        		break;
	        	}
	        }
	        if (getConstantPattern() != null && getConstantPattern().size() > dataInput.get().getPatternCount()) {
	        	// if there are many more constant patterns than patterns (each pattern can
	        	// have a number of constant patters, one for each state) it is less efficient
	        	// to just calculate the TreeLikelihood for constant sites than optimising
	        	Log.debug("switch off constant sites optimisiation: calculating through separate TreeLikelihood category (as in the olden days)");
	        	invariantCategory = -1;
	        	proportionInvariant = 0;
	        	setConstantPattern(null);
	        	categoryRates = m_siteModel.getCategoryRates(null);
	        }
        }        

        this.categoryCount = m_siteModel.getCategoryCount() - (invariantCategory >= 0 ? 1 : 0);
        
        probabilities = new double[(m_nStateCount + 1) * (m_nStateCount + 1)];
        Arrays.fill(probabilities, 1.0);
        matrices = new double[m_nStateCount * m_nStateCount * categoryCount];

        tipCount = treeInput.get().getLeafNodeCount();

        internalNodeCount = m_nNodeCount - tipCount;

        int compactPartialsCount = tipCount;
        if (m_bUseAmbiguities) {
            // if we are using ambiguities then we don't use tip partials
            compactPartialsCount = 0;
        }

        // one partials buffer for each tip and two for each internal node (for store restore)
        partialBufferHelper = new BufferIndexHelper(m_nNodeCount, tipCount);

        // two eigen buffers for each decomposition for store and restore.
        eigenBufferHelper = new BufferIndexHelper(eigenCount, 0);

        // two matrices for each node less the root
        matrixBufferHelper = new BufferIndexHelper(m_nNodeCount, 0);

        // one scaling buffer for each internal node plus an extra for the accumulation, then doubled for store/restore
        scaleBufferHelper = new BufferIndexHelper(getScaleBufferCount(), 0);

        // Attempt to get the resource order from the System Property
        if (resourceOrder == null) {
            resourceOrder = parseSystemPropertyIntegerArray(RESOURCE_ORDER_PROPERTY);
        }
        if (preferredOrder == null) {
            preferredOrder = parseSystemPropertyIntegerArray(PREFERRED_FLAGS_PROPERTY);
        }
        if (requiredOrder == null) {
            requiredOrder = parseSystemPropertyIntegerArray(REQUIRED_FLAGS_PROPERTY);
        }
        if (scalingOrder == null) {
            scalingOrder = parseSystemPropertyStringArray(SCALING_PROPERTY);
        }

        // first set the rescaling scheme to use from the parser
        rescalingScheme = PartialsRescalingScheme.DEFAULT;// = rescalingScheme;
        rescalingScheme = DEFAULT_RESCALING_SCHEME;
        int[] resourceList = null;
        long preferenceFlags = 0;
        long requirementFlags = 0;

        if (scalingOrder.size() > 0) {
            this.rescalingScheme = PartialsRescalingScheme.parseFromString(
                    scalingOrder.get(instanceCount % scalingOrder.size()));
        }

        if (resourceOrder.size() > 0) {
            // added the zero on the end so that a CPU is selected if requested resource fails
            resourceList = new int[]{resourceOrder.get(instanceCount % resourceOrder.size()), 0};
            if (resourceList[0] > 0) {
                preferenceFlags |= BeagleFlag.PROCESSOR_GPU.getMask(); // Add preference weight against CPU
            }
        }

        if (preferredOrder.size() > 0) {
            preferenceFlags = preferredOrder.get(instanceCount % preferredOrder.size());
        }

        if (requiredOrder.size() > 0) {
            requirementFlags = requiredOrder.get(instanceCount % requiredOrder.size());
        }

        if (scaling.get().equals(Scaling.always)) {
        	this.rescalingScheme = PartialsRescalingScheme.ALWAYS;
        }
        if (scaling.get().equals(Scaling.none)) {
        	this.rescalingScheme = PartialsRescalingScheme.NONE;
        }
        
        // Define default behaviour here
        if (this.rescalingScheme == PartialsRescalingScheme.DEFAULT) {
            //if GPU: the default is^H^Hwas dynamic scaling in BEAST, now NONE
            if (resourceList != null && resourceList[0] > 1) {
                //this.rescalingScheme = PartialsRescalingScheme.DYNAMIC;
                this.rescalingScheme = PartialsRescalingScheme.NONE;
            } else { // if CPU: just run as fast as possible
                //this.rescalingScheme = PartialsRescalingScheme.NONE;
                // Dynamic should run as fast as none until first underflow
                this.rescalingScheme = PartialsRescalingScheme.DYNAMIC;
            }
        }

        if (this.rescalingScheme == PartialsRescalingScheme.AUTO) {
            preferenceFlags |= BeagleFlag.SCALING_AUTO.getMask();
            useAutoScaling = true;
        } else {
//                preferenceFlags |= BeagleFlag.SCALING_MANUAL.getMask();
        }
        String r = System.getProperty(RESCALE_FREQUENCY_PROPERTY);
        if (r != null) {
            rescalingFrequency = Integer.parseInt(r);
            if (rescalingFrequency < 1) {
                rescalingFrequency = RESCALE_FREQUENCY;
            }
        }

        if (preferenceFlags == 0 && resourceList == null) { // else determine dataset characteristics
            if (m_nStateCount == 4 && patternCount < 10000) // TODO determine good cut-off
                preferenceFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
        }

//        if (substitutionModel.canReturnComplexDiagonalization()) {
//            requirementFlags |= BeagleFlag.EIGEN_COMPLEX.getMask();
//        }

        instanceCount++;

        try {
	        beagle = BeagleFactory.loadBeagleInstance(
	                tipCount,
	                partialBufferHelper.getBufferCount(),
	                compactPartialsCount,
	                m_nStateCount,
	                patternCount,
	                eigenBufferHelper.getBufferCount(),            // eigenBufferCount
	                matrixBufferHelper.getBufferCount(),
	                categoryCount,
	                scaleBufferHelper.getBufferCount(), // Always allocate; they may become necessary
	                resourceList,
	                preferenceFlags,
	                requirementFlags
	        );
        } catch (Exception e) {
        	beagle = null;
        }
        if (beagle == null) {
            return false;
        }

//        beagle = new BeagleDebugger(beagle, true);
        
        
        InstanceDetails instanceDetails = beagle.getDetails();
        ResourceDetails resourceDetails = null;

        if (instanceDetails != null) {
            resourceDetails = BeagleFactory.getResourceDetails(instanceDetails.getResourceNumber());
            if (resourceDetails != null) {
                StringBuilder sb = new StringBuilder("  Using BEAGLE version: " + BeagleInfo.getVersion()
                		+ " resource ");
                sb.append(resourceDetails.getNumber()).append(": ");
                sb.append(resourceDetails.getName()).append("\n");
                if (resourceDetails.getDescription() != null) {
                    String[] description = resourceDetails.getDescription().split("\\|");
                    for (String desc : description) {
                        if (desc.trim().length() > 0) {
                            sb.append("    ").append(desc.trim()).append("\n");
                        }
                    }
                }
                sb.append("    with instance flags: ").append(instanceDetails.toString());
                Log.info.println(sb.toString());
            } else {
                Log.warning.println("  Error retrieving BEAGLE resource for instance: " + instanceDetails.toString());
                beagle = null;
                return false;
            }
        } else {
        	Log.warning.println("  No external BEAGLE resources available, or resource list/requirements not met, using Java implementation");
            beagle = null;
            return false;
        }
        Log.warning.println("  " + (m_bUseAmbiguities ? "Using" : "Ignoring") + " ambiguities in tree likelihood.");
        Log.warning.println("  " + (m_bUseTipLikelihoods ? "Using" : "Ignoring") + " character uncertainty in tree likelihood.");
        Log.warning.println("  With " + patternCount + " unique site patterns.");

        
        Node [] nodes = treeInput.get().getNodesAsArray();
        for (int i = 0; i < tipCount; i++) {
        	int taxon = getTaxonIndex(nodes[i].getID(), dataInput.get());  
            if (m_bUseAmbiguities || m_bUseTipLikelihoods) {
                setPartials(beagle, i, taxon);
            } else {
                setStates(beagle, i, taxon);
            }
        }

        if (dataInput.get().isAscertained) {
            ascertainedSitePatterns = true;
        }

        double[] patternWeights = new double[patternCount];
        for (int i = 0; i < patternCount; i++) {
            patternWeights[i] = dataInput.get().getPatternWeight(i);
        }
        beagle.setPatternWeights(patternWeights);

        if (this.rescalingScheme == PartialsRescalingScheme.AUTO &&
                resourceDetails != null &&
                (resourceDetails.getFlags() & BeagleFlag.SCALING_AUTO.getMask()) == 0) {
            // If auto scaling in BEAGLE is not supported then do it here
            this.rescalingScheme = PartialsRescalingScheme.DYNAMIC;
            Log.warning.println("  Auto rescaling not supported in BEAGLE, using : " + this.rescalingScheme.getText());
        } else {
        	Log.warning.println("  Using rescaling scheme : " + this.rescalingScheme.getText());
        }

        if (this.rescalingScheme == PartialsRescalingScheme.DYNAMIC) {
            everUnderflowed = false; // If false, BEAST does not rescale until first under-/over-flow.
        }

        updateSubstitutionModel = true;
        updateSiteModel = true;
        // some subst models (e.g. WAG) never become dirty, so set up subst models right now
        setUpSubstModel();
        // set up sitemodel
        
        beagle.setCategoryRates(categoryRates);
        currentCategoryRates = categoryRates;
        currentFreqs = new double[m_nStateCount];
        currentCategoryWeights = new double[categoryRates.length];
        
        return true;
    }

    
    private static List<Integer> parseSystemPropertyIntegerArray(String propertyName) {
        List<Integer> order = new ArrayList<>();
        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    int n = Integer.parseInt(part.trim());
                    order.add(n);
                } catch (NumberFormatException nfe) {
                	Log.warning.println("Invalid entry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
    }

    private static List<String> parseSystemPropertyStringArray(String propertyName) {

        List<String> order = new ArrayList<>();

        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    String s = part.trim();
                    order.add(s);
                } catch (NumberFormatException nfe) {
                	Log.warning.println("Invalid getEigenDecompositionentry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
    }
    
    
    protected int getScaleBufferCount() {
        return internalNodeCount + 1;
    }

    /**
     * Sets the partials from a sequence in an alignment.
     *
     * @param beagle        beagle
     * @param nodeIndex     nodeIndex
     * @param taxon the taxon
     */
    protected final void setPartials(Beagle beagle,
                                     int nodeIndex, int taxon) {
        Alignment data = dataInput.get();

        double[] partials = new double[patternCount * m_nStateCount * categoryCount];

        int v = 0;
        for (int i = 0; i < patternCount; i++) {

        	double[] tipProbabilities = data.getTipLikelihoods(taxon,i);
            if (tipProbabilities != null) {
            	for (int state = 0; state < m_nStateCount; state++) {
            		partials[v++] = tipProbabilities[state];
            	}
            }
            else {
            	int stateCount = data.getPattern(taxon, i);
                boolean[] stateSet = data.getStateSet(stateCount);
                for (int state = 0; state < m_nStateCount; state++) {
                	 partials[v++] = (stateSet[state] ? 1.0 : 0.0);                
                }
            }
        }

        // if there is more than one category then replicate the partials for each
        int n = patternCount * m_nStateCount;
        int k = n;
        for (int i = 1; i < categoryCount; i++) {
            System.arraycopy(partials, 0, partials, k, n);
            k += n;
        }

        beagle.setPartials(nodeIndex, partials);
    }

    public int getPatternCount() {
        return patternCount;
    }

    void setUpSubstModel() {
        // we are currently assuming a no-category model...
        // TODO More efficient to update only the substitution model that changed, instead of all
        for (int i = 0; i < eigenCount; i++) {
        	if (!models[i].canReturnComplexDiagonalization()) {
	            //EigenDecomposition ed = m_substitutionModel.getEigenDecomposition(i, 0);
	            EigenDecomposition ed = models[i].getEigenDecomposition(null);
	
	            eigenBufferHelper.flipOffset(i);
	
	            beagle.setEigenDecomposition(
	                    eigenBufferHelper.getOffsetIndex(i),
	                    ed.getEigenVectors(),
	                    ed.getInverseEigenVectors(),
	                    ed.getEigenValues());
	        }
    	}
    }

    /**
     * Sets the partials from a sequence in an alignment.
     *
     * @param beagle        beagle
     * @param nodeIndex     nodeIndex
     * @param taxon         the taxon
     */
    protected final void setStates(Beagle beagle,
                                   int nodeIndex, int taxon) {
        Alignment data = dataInput.get();
        int i;

        int[] states = new int[patternCount];

        for (i = 0; i < patternCount; i++) {
            int code = data.getPattern(taxon, i);
            int[] statesForCode = data.getDataType().getStatesForCode(code);
            if (statesForCode.length==1)
                states[i] = statesForCode[0];
            else
                states[i] = code; // Causes ambiguous states to be ignored.
        }

        beagle.setTipStates(nodeIndex, states);
    }

    /**
     *
     * @param taxon the taxon name as a string
     * @param data the alignment
     * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
     *         or -1 if the taxon is not in the alignment.
     */
    private int getTaxonIndex(String taxon, Alignment data) {    	
        int taxonIndex = data.getTaxonIndex(taxon);
        if (taxonIndex == -1) {
        	if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
            }
            if (taxonIndex == -1) {
            	throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
            }
        }
        return taxonIndex;
	}
    
    
//    public void setStates(int tipIndex, int[] states) {
//        System.err.println("BTL:setStates");
//        beagle.setTipStates(tipIndex, states);
//        makeDirty();
//    }
//
//    public void getStates(int tipIndex, int[] states) {
//        System.err.println("BTL:getStates");
//        beagle.getTipStates(tipIndex, states);
//    }


    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
	protected boolean requiresRecalculation() {
        hasDirt = Tree.IS_CLEAN;
        substModelThrehold = Double.POSITIVE_INFINITY;
        
        double[] categoryRates = m_siteModel.getCategoryRates(null);
        if (getConstantPattern() != null) {
            double [] tmp = new double [categoryRates.length - 1];
            for (int k = 0; k < invariantCategory; k++) {
            	tmp[k] = categoryRates[k];
            }
            for (int k = invariantCategory + 1; k < categoryRates.length; k++) {
            	tmp[k-1] = categoryRates[k];
            }
            categoryRates = tmp;
        }
        for (int i = 0; i < categoryRates.length; i++) {
        	if (categoryRates[i] != currentCategoryRates[i]) {
        		updateSiteModel = true;
        		break;
        	}
        }
        //updateSiteModel |= m_siteModel.isDirtyCalculation();

        // make sure nodes crossing thresholds that changed get updated
        Node [] nodes = treeInput.get().getNodesAsArray();
        boolean hasDirtyNodes = false;
        for (int i = 0; i < epochDates.getDimension(); i++) {
        	double t = epochDates.getArrayValue(i);
        	double t2 = currentThresholds[i]; 
        	double maxTime = Math.max(t, t2) + 1e-10;
        	if (t2 != t) {
        		for (Node node : nodes) {
        			if (!node.isRoot() && 
        				((node.getHeight() <= t  && node.getParent().getHeight() >= t ) ||
        				 (node.getHeight() <= t2 && node.getParent().getHeight() >= t2))) {
        				// mark all nodes to the root as requiring the transition matrix
        				// to be updated. This is done by setting its m_branchLengths to -1.
        				while (node != null && node.getHeight() <= maxTime) {
            				m_branchLengths[node.getNr()] = -1;
        					node = node.getParent();
        				}
        				hasDirtyNodes = true;
        			}
        		}
        		currentThresholds[i] = t;
        	}
        }
        
        for (int i = 0; i < models.length; i++) {
        	if (((CalculationNode) models[i]).isDirtyCalculation()) {
        		updateSubstitutionModel = true;
        		substModelThrehold = Math.min(i == 0 ? 0 : epochDatesArray[i-1], substModelThrehold);
// System.err.println("updateSubstitutionModel["+i+"] = true " + substModelThrehold);
        	}
        }
        

        
        
        if (dataInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (updateSubstitutionModel || m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            //m_nHasDirt = Tree.IS_FILTHY;
            return true;
        }

        return treeInput.get().somethingIsDirty() || hasDirtyNodes || updateSubstitutionModel;
    }

    /**
     * Stores the additional state other than model components
     */
    @Override
    public void store() {
        partialBufferHelper.storeState();
        eigenBufferHelper.storeState();
        matrixBufferHelper.storeState();

        if (useScaleFactors || useAutoScaling) { // Only store when actually used
            scaleBufferHelper.storeState();
            System.arraycopy(scaleBufferIndices, 0, storedScaleBufferIndices, 0, scaleBufferIndices.length);
//            storedRescalingCount = rescalingCount;
        }
        super.store();
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
        System.arraycopy(currentThresholds, 0, storedThresholds, 0, currentThresholds.length);
    }

    @Override
    public void restore() {
  		updateSiteModel = true; // this is required to upload the categoryRates to BEAGLE after the restore
        
        partialBufferHelper.restoreState();
        eigenBufferHelper.restoreState();
        matrixBufferHelper.restoreState();

        if (useScaleFactors || useAutoScaling) {
            scaleBufferHelper.restoreState();
            int[] tmp2 = storedScaleBufferIndices;
            storedScaleBufferIndices = scaleBufferIndices;
            scaleBufferIndices = tmp2;
//            rescalingCount = storedRescalingCount;
        }

//        updateRestrictedNodePartials = true;
        super.restore();
        //double[] tmp = m_branchLengths;
        //m_branchLengths = storedBranchLengths;
        //storedBranchLengths = tmp;
        double[] tmp = currentThresholds;
        currentThresholds = storedThresholds;
        storedThresholds = tmp;
    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    @Override
    public double calculateLogP() {
    	
    	

        if (patternLogLikelihoods == null) {
            patternLogLikelihoods = new double[patternCount];
        }

        if (matrixUpdateIndices == null) {
            matrixUpdateIndices = new int[eigenCount][m_nNodeCount];
            branchLengths = new double[eigenCount][m_nNodeCount];
            branchUpdateCount = new int[eigenCount];
            scaleBufferIndices = new int[internalNodeCount];
            storedScaleBufferIndices = new int[internalNodeCount];
        }

        if (operations == null) {
            operations = new int[epochCount][internalNodeCount * Beagle.OPERATION_TUPLE_SIZE];
            operationCount = new int[epochCount];
        }

        recomputeScaleFactors = false;

        if (this.rescalingScheme == PartialsRescalingScheme.ALWAYS) {
            useScaleFactors = true;
            recomputeScaleFactors = true;
        } else if (this.rescalingScheme == PartialsRescalingScheme.DYNAMIC && everUnderflowed) {
            useScaleFactors = true;
            if (rescalingCountInner < RESCALE_TIMES) {
                recomputeScaleFactors = true;
                hasDirt = Tree.IS_FILTHY;// makeDirty();
//                System.err.println("Recomputing scale factors");
            }

            rescalingCountInner++;
            rescalingCount++;
            if (rescalingCount > RESCALE_FREQUENCY) {
                rescalingCount = 0;
                rescalingCountInner = 0;
            }
        } else if (this.rescalingScheme == PartialsRescalingScheme.DELAYED && everUnderflowed) {
            useScaleFactors = true;
            recomputeScaleFactors = true;
            hasDirt = Tree.IS_FILTHY;
            rescalingCount++;
        }

        //for (int i = 0; i < eigenCount; i++) {
        //   branchUpdateCount[i] = 0;
        //}
        operationListCount = 0;

        // operationCount[0] = 0;
		reset();

        final Node root = treeInput.get().getRoot();
        traverse(root, null, true);

        if (updateSubstitutionModel) {
            setUpSubstModel();
        }

        if (updateSiteModel) {
            double[] categoryRates = m_siteModel.getCategoryRates(null);
            if (getConstantPattern() != null) {
	            double [] tmp = new double [categoryRates.length - 1];
	            for (int k = 0; k < invariantCategory; k++) {
	            	tmp[k] = categoryRates[k];
	            }
	            for (int k = invariantCategory + 1; k < categoryRates.length; k++) {
	            	tmp[k-1] = categoryRates[k];
	            }
	            categoryRates = tmp;
	        }
            for (int i = 0; i < categoryRates.length; i++) {
            	if (categoryRates[i] != currentCategoryRates[i]) {
                    beagle.setCategoryRates(categoryRates);
                    i = categoryRates.length;
            	}
            }
            currentCategoryRates = categoryRates;
        }

        for (int i = 0; i < eigenCount; i++) {
            if (!models[i].canReturnComplexDiagonalization()) {
                if (branchUpdateCount[i] > 0) {
                    beagle.updateTransitionMatrices(
                            eigenBufferHelper.getOffsetIndex(i),
                            matrixUpdateIndices[i],
                            null,
                            null,
                            branchLengths[i],
                            branchUpdateCount[i]);
                }
            }
        }

//        if (COUNT_TOTAL_OPERATIONS) {
//            for (int i = 0; i < eigenCount; i++) {
//                totalMatrixUpdateCount += branchUpdateCount[i];
//            }
//            
//            for (int i = 0; i <= numRestrictedPartials; i++) {
//                totalOperationCount += operationCount[i];
//            }
//        }

        double logL;
        boolean done;
        boolean firstRescaleAttempt = true;

        do {
            for (int i = 0; i < eigenCount; i++) {
            	if (operationCount[i] > 0) {
            		beagle.updatePartials(operations[i], operationCount[i], Beagle.NONE);
            	}
            }

            int rootIndex = partialBufferHelper.getOffsetIndex(root.getNr());

            double[] categoryWeights = m_siteModel.getCategoryProportions(null);
            if (getConstantPattern() != null) {
	            double [] tmp = new double [categoryWeights.length - 1];
	            for (int k = 0; k < invariantCategory; k++) {
	            	tmp[k] = categoryWeights[k];
	            }
	            for (int k = invariantCategory + 1; k < categoryWeights.length; k++) {
	            	tmp[k-1] = categoryWeights[k];
	            }
	            categoryWeights = tmp;
            }
            
            double[] frequencies = null;
            	
        	// Get the frequencies from the oldest epoch that is younger than the root
        	double rootHeight = this.treeInput.get().getRoot().getHeight();
        	for (int i = epochCount-1; i >= 0; i --) {
        		if (i == 0 || epochDatesArray[i-1] <= rootHeight) {
        			
        			//Log.warning("top epoch is " + i);
        			
        			if (i == 0 && rootFrequenciesInput.get() != null) {
        				frequencies = rootFrequenciesInput.get().getFreqs();
        			}else {
        				frequencies = models[i].getFrequencies();
        			}
            		break;
        		}
        	}
        	
        	
        	
			//System.out.println(this.getClass().getSimpleName() + " " + epochDatesArray[0]);
//          System.out.print("epoch freqs: ");
//          for (int i = 0; i < frequencies.length; i ++) {
//          	System.out.print(frequencies[i] + ",");
//          }
//          System.out.println();
        	


            int cumulateScaleBufferIndex = Beagle.NONE;
            if (useScaleFactors) {

                if (recomputeScaleFactors) {
                    scaleBufferHelper.flipOffset(internalNodeCount);
                    cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
                    beagle.resetScaleFactors(cumulateScaleBufferIndex);
                    beagle.accumulateScaleFactors(scaleBufferIndices, internalNodeCount, cumulateScaleBufferIndex);
                } else {
                    cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
                }
            } else if (useAutoScaling) {
                beagle.accumulateScaleFactors(scaleBufferIndices, internalNodeCount, Beagle.NONE);
            }

            // these could be set only when they change but store/restore would need to be considered
            
            for (int i = 0; i < categoryWeights.length; i++) {
            	if (categoryWeights[i] != currentCategoryWeights[i]) {
                    beagle.setCategoryWeights(0, categoryWeights);
            		i = categoryWeights.length;
            	}
            }
            currentCategoryWeights = categoryWeights;
            for (int i = 0; i < frequencies.length; i++) {
            	if (frequencies[i] != currentFreqs[i]) {
                    beagle.setStateFrequencies(0, frequencies);
            		i = frequencies.length;
            	}
            }
            currentFreqs = frequencies;

            double[] sumLogLikelihoods = new double[1];

            beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0},
                    new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);

            
            //System.out.println(this.getClass().getSimpleName() + " " + sumLogLikelihoods[0]);
            logL = sumLogLikelihoods[0];

            if (ascertainedSitePatterns) {
                // Need to correct for ascertainedSitePatterns
                beagle.getSiteLogLikelihoods(patternLogLikelihoods);
                logL = getAscertainmentCorrectedLogLikelihood(dataInput.get(),
                        patternLogLikelihoods, dataInput.get().getWeights(), frequencies);
            } else if (invariantCategory >= 0) {
                beagle.getSiteLogLikelihoods(patternLogLikelihoods);
                int [] patternWeights = dataInput.get().getWeights();
                proportionInvariant = m_siteModel.getProportionInvariant();
                
                
    	        for (int k : getConstantPattern()) {
    	        	int i = k / m_nStateCount;
    	        	int j = k % m_nStateCount;
    	        	patternLogLikelihoods[i] = (Math.log(Math.exp(patternLogLikelihoods[i]) + proportionInvariant * frequencies[j]));
    	        }
        	
	            logL = 0.0;
	            for (int i = 0; i < patternCount; i++) {
	                logL += patternLogLikelihoods[i] * patternWeights[i];
	            }
            }

            if (Double.isNaN(logL) || Double.isInfinite(logL)) {
                everUnderflowed = true;
                logL = Double.NEGATIVE_INFINITY;

                if (firstRescaleAttempt && (rescalingScheme == PartialsRescalingScheme.DYNAMIC || rescalingScheme == PartialsRescalingScheme.DELAYED)) {
                    // we have had a potential under/over flow so attempt a rescaling                	
                	useScaleFactors = true;
                    recomputeScaleFactors = true;

                    // for (int i = 0; i < eigenCount; i++) {
                    //   branchUpdateCount[i] = 0;
                    //}
                    //operationCount[0] = 0;
                    
                    reset();

                    // traverse again but without flipping partials indices as we
                    // just want to overwrite the last attempt. We will flip the
                    // scale buffer indices though as we are recomputing them.
                    traverse(root, null, false);

                    done = false; // Run through do-while loop again
                    firstRescaleAttempt = false; // Only try to rescale once
                } else {
                    // we have already tried a rescale, not rescaling or always rescaling
                    // so just return the likelihood...
                    done = true;
                }
            } else {
                done = true; // No under-/over-flow, then done
            }

        } while (!done);

        // If these are needed...
        //beagle.getSiteLogLikelihoods(patternLogLikelihoods);

        //********************************************************************
        // after traverse all nodes and patterns have been updated --
        //so change flags to reflect this.
//        for (int i = 0; i < m_nNodeCount; i++) {
//            updateNode[i] = false;
//        }

        updateSubstitutionModel = false;
        updateSiteModel = false;
        //********************************************************************

        logP = logL;
        return logL;
    }

	private void reset() {
	    Arrays.fill(branchUpdateCount, 0);
	    Arrays.fill(operationCount, 0);
	
	    for (int i = 0;i < epochDatesArray.length; i++) {
	    	epochDatesArray[i] = epochDates.getArrayValue(i);
	    }
	}
	
    protected void setPartials(int number, double[] partials) {
        beagle.setPartials(partialBufferHelper.getOffsetIndex(number), partials);
    }

    private double getAscertainmentCorrectedLogLikelihood(Alignment patternList,
                                                          double[] patternLogLikelihoods,
                                                          int[] patternWeights,
                                                          double [] frequencies) {
    	if (getConstantPattern() != null) {
	        proportionInvariant = m_siteModel.getProportionInvariant();
	        for (int k : getConstantPattern()) {
	        	int i = k / m_nStateCount;
	        	int j = k % m_nStateCount;
	        	patternLogLikelihoods[i] = (Math.log(Math.exp(patternLogLikelihoods[i]) + proportionInvariant * frequencies[j]));
	        }
    	}
    	
        double logL = 0.0;
        double ascertainmentCorrection = patternList.getAscertainmentCorrection(patternLogLikelihoods);
        for (int i = 0; i < patternCount; i++) {
            logL += (patternLogLikelihoods[i] - ascertainmentCorrection) * patternWeights[i];
        }
        return logL;
    }

    /**
     * Traverse the tree calculating partial likelihoods.
     *
     * @param node           node
     * @param operatorNumber operatorNumber
     * @param flip           flip
     * @return boolean
     */
    private int traverse(Node node, int[] operatorNumber, boolean flip) {

        int nodeNum = node.getNr();

        //Node parent = node.getParent();

        if (operatorNumber != null) {
            operatorNumber[0] = -1;
        }

        // First update the transition probability matrix(ices) for this branch
        int update = (node.isDirty() | hasDirt); 
        if (!node.isRoot() && node.getParent().getHeight() >= substModelThrehold) {
        	update |= Tree.IS_DIRTY;
        }
//        if (parent!=null) {
//        	update |= parent.isDirty();
//        }
        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;
		double startTime = node.getHeight();
		int epoch = Arrays.binarySearch(epochDatesArray, startTime);
		if (epoch < 0) {
			epoch = -epoch-1;
		}

		if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeNum])) {
            m_branchLengths[nodeNum] = branchTime;
            if (branchTime < 0.0) {
                throw new RuntimeException("Negative branch length: " + branchTime);
            }

    		double endTime = node.getParent().getHeight();
    		int endEpoch = Arrays.binarySearch(epochDatesArray, endTime);
    		if (endEpoch < 0) {
    			endEpoch = -endEpoch-1;
    		}

            
            if (flip) {
                // first flip the matrixBufferHelper
                matrixBufferHelper.flipOffset(nodeNum);
            }

            // then set which matrix to update
            // final int eigenIndex = 0;// = m_substitutionModel.getBranchIndex(node);
            final int updateCount = branchUpdateCount[epoch];
            matrixUpdateIndices[epoch][updateCount] = matrixBufferHelper.getOffsetIndex(nodeNum);

            if (epoch != endEpoch || models[epoch].canReturnComplexDiagonalization()) {
                for (int i = 0; i < this.categoryCount; i++) {
                    final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                    getTransitionProbabilities(endEpoch, epoch, node, node.getParent().getHeight(), node.getHeight(), jointBranchRate, probabilities);
                    //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
                    System.arraycopy(probabilities, 0, matrices,  m_nStateCount * m_nStateCount * i, m_nStateCount * m_nStateCount);
                }
            	int matrixIndex = matrixBufferHelper.getOffsetIndex(nodeNum);
            	beagle.setTransitionMatrix(matrixIndex, matrices, 1);

            	// do not BEAGLE exponentiate this matrix
            	branchUpdateCount[epoch]--;
            }

            branchLengths[epoch][updateCount] = branchTime;
            branchUpdateCount[epoch]++;

            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            Node child1 = node.getLeft();
            final int[] op1 = {-1};
            final int update1 = traverse(child1, op1, flip);

            Node child2 = node.getRight();
            final int[] op2 = {-1};
            final int update2 = traverse(child2, op2, flip);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                int x = operationCount[epoch] * Beagle.OPERATION_TUPLE_SIZE;

                if (flip) {
                    // first flip the partialBufferHelper
                    partialBufferHelper.flipOffset(nodeNum);
                }

                final int[] operations = this.operations[epoch];

                // System.out.println(x + " " + nodeNum + " " + partialBufferHelper.getOffsetIndex(nodeNum));
                operations[x] = partialBufferHelper.getOffsetIndex(nodeNum);

                if (useScaleFactors) {
                    // get the index of this scaling buffer
                    int n = nodeNum - tipCount;

                    if (recomputeScaleFactors) {
                        // flip the indicator: can take either n or (internalNodeCount + 1) - n
                        scaleBufferHelper.flipOffset(n);

                        // store the index
                        scaleBufferIndices[n] = scaleBufferHelper.getOffsetIndex(n);

                        operations[x + 1] = scaleBufferIndices[n]; // Write new scaleFactor
                        operations[x + 2] = Beagle.NONE;

                    } else {
                        operations[x + 1] = Beagle.NONE;
                        operations[x + 2] = scaleBufferIndices[n]; // Read existing scaleFactor
                    }

                } else {

                    if (useAutoScaling) {
                        scaleBufferIndices[nodeNum - tipCount] = partialBufferHelper.getOffsetIndex(nodeNum);
                    }
                    operations[x + 1] = Beagle.NONE; // Not using scaleFactors
                    operations[x + 2] = Beagle.NONE;
                }

                operations[x + 3] = partialBufferHelper.getOffsetIndex(child1.getNr()); // source node 1
                operations[x + 4] = matrixBufferHelper.getOffsetIndex(child1.getNr()); // source matrix 1
                operations[x + 5] = partialBufferHelper.getOffsetIndex(child2.getNr()); // source node 2
                operations[x + 6] = matrixBufferHelper.getOffsetIndex(child2.getNr()); // source matrix 2

                operationCount[epoch]++;

                update |= (update1 | update2);

            }
        }

        return update;

    }

    
    private void getTransitionProbabilities(int startEpoch, int endEpoch, Node node, double startTime, double endTime,
			double jointBranchRate, double[] probabilities) {
    	
    	
    	//Log.warning(startTime + "-" + endTime + " at transition " + epochDatesArray[0]);

		double [] p2 = new double[probabilities.length];
		double [] p3 = new double[probabilities.length];
		for (int k = startEpoch; k > endEpoch; k--) {
			//Log.warning(models[k].getClass().getSimpleName() + " -> " + startTime + "-" + epochDatesArray[k-1]);
			models[k].getTransitionProbabilities(node, startTime, epochDatesArray[k-1], jointBranchRate, p2);
			if (k < startEpoch) {
				System.arraycopy(probabilities, 0, p3, 0, probabilities.length);
				multiplyMatrices(p3, p2, probabilities);
			} else {
				System.arraycopy(p2, 0, probabilities, 0, probabilities.length);
			}
			startTime = epochDatesArray[k-1];
		}

		//Log.warning(models[endEpoch].getClass().getSimpleName() + " -> " + startTime + "-" + endTime);
		models[endEpoch].getTransitionProbabilities(node, startTime, endTime, jointBranchRate, p2);
		if (startEpoch != endEpoch) {
			System.arraycopy(probabilities, 0, p3, 0, probabilities.length);
			multiplyMatrices(p3, p2, probabilities);
		} else {
			System.arraycopy(p2, 0, probabilities, 0, probabilities.length);
		}
    }
    
    

	public void multiplyMatrices(double[] left, double[] right, double[] result) {
		
		int index3 = 0;
		for (int rowNum = 0; rowNum < m_nStateCount; rowNum ++) {
			
			for (int colNum = 0; colNum < m_nStateCount; colNum ++) {
				
				double sum = 0;
				
				// Move along the columns of left and rows of right 
				int index1 = rowNum * m_nStateCount;
				int index2 = colNum;
				for (int pos = 0; pos < m_nStateCount; pos++) {
					double term1 = left[index1];
					double term2 = right[index2]; 
					//double term2 = right[index1]; // Transposed
					sum += term1*term2;

					index1 ++;
					index2 += m_nStateCount;
				}

				result[index3] = sum;
				index3 ++;
			}
		}
	}
	

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    private int eigenCount;
    private int[][] matrixUpdateIndices;
    private double[][] branchLengths;
    private int[] branchUpdateCount;
    private int[] scaleBufferIndices;
    private int[] storedScaleBufferIndices;

    private int[][] operations;
    private int operationListCount;
    private int[] operationCount;

    protected BufferIndexHelper partialBufferHelper;
    public BufferIndexHelper getPartialBufferHelper() {return partialBufferHelper;}
    
    private /*final*/ BufferIndexHelper eigenBufferHelper;
    protected BufferIndexHelper matrixBufferHelper;
    public BufferIndexHelper getMatrixBufferHelper() {return matrixBufferHelper;}
    protected BufferIndexHelper scaleBufferHelper;

    protected /*final*/ int tipCount;
    protected /*final*/ int internalNodeCount;
    protected /*final*/ int patternCount;

    private PartialsRescalingScheme rescalingScheme = DEFAULT_RESCALING_SCHEME;
    private int rescalingFrequency = RESCALE_FREQUENCY;
    protected boolean useScaleFactors = false;
    private boolean useAutoScaling = false;
    private boolean recomputeScaleFactors = false;
    private boolean everUnderflowed = false;
    private int rescalingCount = 0;
    private int rescalingCountInner = 0;

    
    /**
     * the number of rate categories
     */
    protected int categoryCount;

    /**
     * an array used to transfer tip partials
     */
    protected double[] tipPartials;

    /**
     * the BEAGLE library instance
     */
    protected Beagle beagle;
    
    public Beagle getBeagle() {return beagle;}

    /**
     * Flag to specify that the substitution model has changed
     */
    protected boolean updateSubstitutionModel;
    protected boolean storedUpdateSubstitutionModel;

    /**
     * Flag to specify that the site model has changed
     */
    protected boolean updateSiteModel;
    protected boolean storedUpdateSiteModel;

    /**
     * Flag to specify if site patterns are acertained
     */

    private boolean ascertainedSitePatterns = false;

    public class BufferIndexHelper {
        /**
         * @param maxIndexValue the number of possible input values for the index
         * @param minIndexValue the minimum index value to have the mirrored buffers
         */
        BufferIndexHelper(int maxIndexValue, int minIndexValue) {
            this.maxIndexValue = maxIndexValue;
            this.minIndexValue = minIndexValue;

            offsetCount = maxIndexValue - minIndexValue;
            indexOffsets = new int[offsetCount];
            storedIndexOffsets = new int[offsetCount];
        }

        public int getBufferCount() {
            return 2 * offsetCount + minIndexValue;
        }

        void flipOffset(int i) {
            if (i >= minIndexValue) {
                indexOffsets[i - minIndexValue] = offsetCount - indexOffsets[i - minIndexValue];
            } // else do nothing
        }

        public int getOffsetIndex(int i) {
            if (i < minIndexValue) {
                return i;
            }
            return indexOffsets[i - minIndexValue] + i;
        }

        void getIndices(int[] outIndices) {
            for (int i = 0; i < maxIndexValue; i++) {
                outIndices[i] = getOffsetIndex(i);
            }
        }

        void storeState() {
            System.arraycopy(indexOffsets, 0, storedIndexOffsets, 0, indexOffsets.length);

        }

        void restoreState() {
            int[] tmp = storedIndexOffsets;
            storedIndexOffsets = indexOffsets;
            indexOffsets = tmp;
        }

        private final int maxIndexValue;
        private final int minIndexValue;
        private final int offsetCount;

        private int[] indexOffsets;
        private int[] storedIndexOffsets;

    } // class BufferIndexHelper

    public enum PartialsRescalingScheme {
        DEFAULT("default"), // whatever our current favourite default is
        NONE("none"),       // no scaling
        DYNAMIC("dynamic"), // rescale when needed and reuse scaling factors
        ALWAYS("always"),   // rescale every node, every site, every time - slow but safe
        DELAYED("delayed"), // postpone until first underflow then switch to 'always'
        AUTO("auto");       // BEAGLE automatic scaling - currently playing it safe with 'always'
//        KICK_ASS("kickAss"),// should be good, probably still to be discovered

        PartialsRescalingScheme(String text) {
            this.text = text;
        }

        public String getText() {
            return text;
        }

        private final String text;

        public static PartialsRescalingScheme parseFromString(String text) {
            for (PartialsRescalingScheme scheme : PartialsRescalingScheme.values()) {
                if (scheme.getText().compareToIgnoreCase(text) == 0)
                    return scheme;
            }
            return DEFAULT;
        }
    }
    
    @Override
    public double [] getPatternLogLikelihoods() {
        beagle.getSiteLogLikelihoods(patternLogLikelihoods);
		return patternLogLikelihoods.clone();
	}
    

}