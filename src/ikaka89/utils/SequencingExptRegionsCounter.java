package ikaka89.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.math.stats.StatUtil;
import org.seqcode.projects.seed.stats.FeatureStatistics;

import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;



/**
 * Simple class to count the number of tags mapping to a given set of regions for a given list of sequencing experiments
 * @author akshaykakumanu
 *
 */
public class SequencingExptRegionsCounter {
	
	protected ExptConfig econf;
	protected GenomeConfig gconf;
	
	protected double minSigCtrlFoldDifference;
	protected List<Region> regions;
	protected String label;
	
	public SequencingExptRegionsCounter(ExptConfig econfig, GenomeConfig gconfig) {
		econf = econfig;
		gconf = gconfig;
	}
	
	public void printPooledCounts(){
		if(!econf.getCacheAllData()){
			Collections.sort(regions);
		}
		ExperimentManager manager = new ExperimentManager(econf);
		double[][] counts = new double[regions.size()][manager.getConditions().size()];
		StringBuilder header = new StringBuilder();
		header.append("#Region");
		header.append("\t");
		int ConditionCount = 0;
		for(ExperimentCondition ec :  manager.getConditions()){
			header.append(ec.getName());
			header.append("\t");
			int regionCount = 0;
			for(Region r : regions){
				for(Sample s : ec.getSignalSamples()){
					counts[regionCount][ConditionCount] += s.countHits(r);
				}
				regionCount++;
			}
			ConditionCount++;
		}
		header.deleteCharAt(header.length()-1);
		System.out.println(header.toString());
		StringBuilder CountsSb = new StringBuilder();
		for(int r=0; r<regions.size(); r++){
			CountsSb.append(regions.get(r).getLocationString());
			CountsSb.append("\t");
			for(int c=0; c<manager.getConditions().size(); c++){
				CountsSb.append(counts[r][c]);
				CountsSb.append("\t");
			}
			CountsSb.deleteCharAt(CountsSb.length()-1);
			CountsSb.append("\n");
		}
		System.out.println(CountsSb.toString());
		
		manager.close();
		
	}
	
	public void printCounts(){
		if(!econf.getCacheAllData()){
			Collections.sort(regions);
		}
		ExperimentManager manager = new ExperimentManager(econf);
		float[][] Counts = new float[regions.size()][manager.getSamples().size()];
		StringBuilder header = new StringBuilder();
		header.append("#Region");
		header.append("\t");
		int SampleCount = 0;
		for(Sample sample : manager.getSamples()){
			header.append(sample.getName());
			header.append("\t");
			int RegionCount = 0;
			for(Region r : regions){
				Counts[RegionCount][SampleCount] = sample.countHits(r);
				RegionCount++;
			}
			SampleCount++;
		}
		
		
		header.deleteCharAt(header.length()-1);
		System.out.println(header.toString());
		
		StringBuilder CountsSb = new StringBuilder();
		for(int r=0; r<regions.size(); r++){
			CountsSb.append(regions.get(r).getLocationString());
			CountsSb.append("\t");
			for(int c=0; c<manager.getSamples().size(); c++){
				CountsSb.append(Counts[r][c]);
				CountsSb.append("\t");
			}
			CountsSb.deleteCharAt(CountsSb.length()-1);
			CountsSb.append("\n");
		}
		System.out.println(CountsSb.toString());
		
		manager.close();
	}
	
	public void printBinomialPvals(){
		if(!econf.getCacheAllData()){
			Collections.sort(regions);
		}
		ExperimentManager manager = new ExperimentManager(econf);
		double[][] pvals = new double[regions.size()][manager.getConditions().size()];
		StringBuilder header = new StringBuilder();
		header.append("#Region");
		header.append("\t");
		int ConditionCount = 0;
		for(ExperimentCondition ec :  manager.getConditions()){
			header.append(ec.getName());
			header.append("\t");
			int regionCount = 0;
			for(Region r : regions){
				double pooledSigCount = 0;
				double pooledCtrlCount =0;
				for(Sample s : ec.getSignalSamples()){
					pooledSigCount += s.countHits(r);
				}
				
				if(ec.getControlSamples().size() != 0 ){
					for(Sample ctrls : ec.getControlSamples()){
						pooledCtrlCount += ctrls.countHits(r);
					}
				}else{
					pooledCtrlCount += ec.getTotalSignalCount()*r.getWidth()/econf.getMappableGenomeLength();
				}
				pooledCtrlCount = pooledCtrlCount*ec.getPooledSampleControlScaling();
				pvals[regionCount][ConditionCount] = StatUtil.binomialPValue(pooledCtrlCount, (pooledCtrlCount+pooledSigCount), 1/(1+this.minSigCtrlFoldDifference));
				
				regionCount++;
			}
			ConditionCount++;
		}
		header.deleteCharAt(header.length()-1);
		System.out.println(header.toString());
		StringBuilder pvalsSb = new StringBuilder();
		for(int r=0; r<regions.size(); r++){
			pvalsSb.append(regions.get(r).getLocationString());
			pvalsSb.append("\t");
			for(int c=0; c<manager.getConditions().size(); c++){
				pvalsSb.append(pvals[r][c]);
				pvalsSb.append("\t");
			}
			pvalsSb.deleteCharAt(pvalsSb.length()-1);
			pvalsSb.append("\n");
		}
		System.out.println(pvalsSb.toString());
		
	}
	
	public void printInputFoldChanges(){
		if(!econf.getCacheAllData()){
			Collections.sort(regions);
		}
		ExperimentManager manager = new ExperimentManager(econf);
		
		ArrayList<String> header = new ArrayList<String>();
		
		HashMap<String,double[]> pooledSigCounts = new HashMap<String,double[]>();
		HashMap<String,double[]> pooledCtrlCounts = new HashMap<String,double[]>();
		
		//double[] pooledSigCounts= new double[regions.size()];
		//double[] pooledCtrlCounts = new double[regions.size()];
		
		for(ExperimentCondition ec: manager.getConditions()){
			
			header.add(ec.getName());
			if(!pooledSigCounts.containsKey(ec.getName())){
				pooledSigCounts.put(ec.getName(), new double[regions.size()]);
				pooledCtrlCounts.put(ec.getName(), new double[regions.size()]);
			}
			
			for(Sample sigs : ec.getSignalSamples()){
				int RegionCount=0;
				for(Region r : regions){
					pooledSigCounts.get(ec.getName())[RegionCount] = pooledSigCounts.get(ec.getName())[RegionCount] + sigs.countHits(r);
					RegionCount++;
				}
			}
			
			if(ec.getControlSamples().size() != 0 ){
				for(Sample ctrls : ec.getControlSamples()){
					int RegionCount=0;
					for(Region r : regions){
						pooledCtrlCounts.get(ec.getName())[RegionCount] = pooledCtrlCounts.get(ec.getName())[RegionCount] + ctrls.countHits(r);
						RegionCount++;
					}
				}
			}else{
				int RegionCount=0;
				for(Region r : regions){
					pooledCtrlCounts.get(ec.getName())[RegionCount] = ec.getTotalSignalCount()*r.getWidth()/econf.getMappableGenomeLength();
					RegionCount++;
				}
			}
			// Scale the control counts
			int RegionCount=0;
			for(Region r : regions){
				pooledCtrlCounts.get(ec.getName())[RegionCount] = pooledCtrlCounts.get(ec.getName())[RegionCount]*ec.getPooledSampleControlScaling();
				RegionCount++;
			}
			
			
		}
		
		StringBuilder fcSb = new StringBuilder();
		fcSb.append("#Region");fcSb.append("\t");
		
		for(String h : header){
			fcSb.append(h);fcSb.append("\t");
		}
		fcSb.deleteCharAt(fcSb.length()-1);
		fcSb.append("\n");
		
		for(int r=0; r<regions.size(); r++){
			fcSb.append(regions.get(r).getLocationString());
			fcSb.append("\t");
			for(String ec : header){
				fcSb.append(pooledSigCounts.get(ec)[r]/pooledCtrlCounts.get(ec)[r]);fcSb.append("\t");
			}
			fcSb.deleteCharAt(fcSb.length()-1);
			fcSb.append("\n");
		}
		System.out.println(fcSb.toString());
		
		manager.close();
	}
	
	
	public void printInputEnrichment(){
		if(!econf.getCacheAllData()){
			Collections.sort(regions);
		}
		ExperimentManager manager = new ExperimentManager(econf);
		double[] enrichment = new double[regions.size()];
		StringBuilder header = new StringBuilder();
		header.append("#Region");
		header.append("\t");
		if(label != null){
			header.append(label);
		}
		
		for(ExperimentCondition ec: manager.getConditions()){
			if(label == null){
				header.append(ec.getName());
				header.append("\t");
			}
			double[] pooledSigCounts= new double[regions.size()];
			double[] pooledCtrlCounts = new double[regions.size()];
			
			for(Sample sigs : ec.getSignalSamples()){
				int RegionCount=0;
				for(Region r : regions){
					pooledSigCounts[RegionCount] = pooledSigCounts[RegionCount] + sigs.countHits(r);
					RegionCount++;
				}
			}
			
			if(ec.getControlSamples().size() != 0 ){
				for(Sample ctrls : ec.getControlSamples()){
					int RegionCount=0;
					for(Region r : regions){
						pooledCtrlCounts[RegionCount] = pooledCtrlCounts[RegionCount] + ctrls.countHits(r);
						RegionCount++;
					}
				}
			}else{
				int RegionCount=0;
				for(Region r : regions){
					pooledCtrlCounts[RegionCount] = ec.getTotalSignalCount()*r.getWidth()/econf.getMappableGenomeLength();
					RegionCount++;
				}
			}
			// Scale the control counts
			int RegionCount=0;
			for(Region r : regions){
				pooledCtrlCounts[RegionCount] = pooledCtrlCounts[RegionCount]*ec.getPooledSampleControlScaling();
				RegionCount++;
			}
			
			
			// Do bionomial testing 
			RegionCount=0;
			for(Region r : regions){
				double pval = StatUtil.binomialPValue(pooledCtrlCounts[RegionCount], (pooledCtrlCounts[RegionCount]+pooledSigCounts[RegionCount]), 1/(1+this.minSigCtrlFoldDifference));
				enrichment[RegionCount] = -1*Math.log10(pval);
				RegionCount++;
			}
		}
		if(label == null){header.deleteCharAt(header.length()-1);}
		
		System.out.println(header.toString());
		
		StringBuilder enrichmentSb = new StringBuilder();
		
		for(int r=0; r<regions.size(); r++){
			enrichmentSb.append(regions.get(r).getLocationString());
			enrichmentSb.append("\t");
			enrichmentSb.append(enrichment[r]);
			enrichmentSb.append("\n");
		}
		System.out.println(enrichmentSb.toString());
		
		manager.close();
	}
	
	
	
	public void setRegions(List<Region> regs){regions = regs;}
	public void setMinSigCntrlFoldDiff(double min){minSigCtrlFoldDifference = min;}
	public void setExptLable(String lab){label = lab;}
	
	public static void main(String[] args){
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		
		ArgParser ap = new ArgParser(args);
		
		SequencingExptRegionsCounter counter = new SequencingExptRegionsCounter(econ, gcon);
		
		int win = Args.parseInteger(args, "win", 200);
		
		String label = ap.hasKey("label") ? ap.getKeyValue("label") : "IP";
		
		double minSigCtrlFoldDifference = Args.parseDouble(args,"minfolddiff",1);
		
		if (!ap.hasKey("refTSSs") && !ap.hasKey("peaks") && !ap.hasKey("regions")){
			System.err.println("Provied either a peaks file or a refTSSs file or a regions !!!");
			System.exit(1);
		}
		
		List<Region> reg = new ArrayList<Region>();
		
		if(ap.hasKey("peaks") || ap.hasKey("refTSSs")){
			@SuppressWarnings("unchecked")
			List<Point> points = (List<Point>) (ap.hasKey("refTSSs") ? RegionFileUtilities.loadStrandedPointFromRefTssFile(gcon.getGenome(), ap.getKeyValue("refTSSs"))
					: RegionFileUtilities.loadPeaksFromPeakFile(gcon.getGenome(), ap.getKeyValue("peaks"), win));
			for(Point p: points){
				reg.add(p.expand(win/2));
			}
		}else{
			reg = RegionFileUtilities.loadRegionsFromPeakFile(gcon.getGenome(), ap.getKeyValue("regions"), -1);
		}
		
		
		counter.setRegions(reg);
		counter.setMinSigCntrlFoldDiff(minSigCtrlFoldDifference);
		counter.setExptLable(label);
		if(ap.hasKey("counts")){
			counter.printCounts();
		}else if(ap.hasKey("enrichment")){
			counter.printInputEnrichment();
		}else if(ap.hasKey("pooledCounts")){
			counter.printPooledCounts();
		}else if(ap.hasKey("enrichmentPvals")){
			counter.printBinomialPvals();
		}else if(ap.hasKey("foldchangesOverinput")){
			counter.printInputFoldChanges();
		}
		
	}
	

}
