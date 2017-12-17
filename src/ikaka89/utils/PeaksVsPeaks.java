package ikaka89.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;


public class PeaksVsPeaks {
	
	protected Genome gen;
	protected List<Point> peaksA;
	protected List<Region> regionsA;
	protected List<Point> peaksB;
	protected List<Region> regionsB;
	
	protected List<List<Point>> peaksBSet;
	protected List<List<Region>> regionsBSet;
	
	protected int overlapD; // min distance between the closest peaks to call them overlapping 
	protected int win; // for expanding a peak into a peak region	
	
	public PeaksVsPeaks(Genome g) {
		gen = g;
	}
	
	
	
	// Settors
	public void setPeaksA(List<Point> a){peaksA = a;}
	public void setPeaksB(List<Point> b){peaksB = b;}
	public void setRegsA(List<Region> a){regionsA = a;}
	public void setRegsB(List<Region>  b){regionsB = b;}
	public void setOverlapD(int o){overlapD =o;}
	public void setWin(int w){win = w;}
	public void setPeaksBSet(List<List<Point>> pset){peaksBSet = pset;}
	public void setRegionsBSet(List<List<Region>> rset){regionsBSet = rset;}
	
	/**
	 * Given a list of peaks; prints the number of other peaks in that list
	 * that are within a distance of rad from each peak
	 * @param rad
	 */
	public void printClosePeaksInAList(int rad){
		HashMap<String, List<Point>> peaksAbyChrom = hashbychrom(peaksA);
		HashMap<Point, Integer> closebyPeaks = new HashMap<Point,Integer>();
		for(String chrom : peaksAbyChrom.keySet()){
			List<Point> currChrPeaks = peaksAbyChrom.get(chrom);
			Collections.sort(currChrPeaks);
			for(Point pa : currChrPeaks){
				int count = 0;
				for(Point pAprime : currChrPeaks){
					if(!pa.getLocationString().equals(pAprime.getLocationString())){
						if(pa.distance(pAprime) < rad)
							count++;
						else
							break;
					}
				}
				closebyPeaks.put(pa, count);
			}
		}

		StringBuilder sb = new StringBuilder();
		for(Point p : closebyPeaks.keySet()){
			sb.append(p.getLocationString());sb.append("\t");sb.append(closebyPeaks.get(p));sb.append("\n");
		}
		sb.deleteCharAt(sb.length()-1);
		System.out.println(sb.toString());
	}

	
	public void printUnique(int rad){
		
		HashMap<String, List<Point>> peaksAbyChrom = hashbychrom(peaksA);
		for(String chrom : peaksAbyChrom.keySet()){
			List<Point> outPeaks = new ArrayList<Point>();
			for(Point pa : peaksAbyChrom.get(chrom)){
				boolean add = true;
				for(Point addedA : outPeaks){
					if(addedA.expand(rad).contains(pa)){
						add = false;
						break;
					}
				}
				if(add)
					outPeaks.add(pa);
			}

			for(Point outP : outPeaks){
				System.out.println(outP.getLocationString());
			}

		}

	}
	
	public void printClosestOffsetPeak(){
		HashMap<String, List<Point>> peaksAbyChrom = hashbychrom(peaksA);
		HashMap<String, List<Point>> peaksBbyChrom = hashbychrom(peaksB);
		for(String chrom : peaksAbyChrom.keySet()){
			for(Point pa : peaksAbyChrom.get(chrom)){
				int mind=Integer.MAX_VALUE;
				Point nearestPeak=null;
				for(Point pb : peaksBbyChrom.get(chrom)){
					if(pb.getChrom().equals(pa.getChrom())){
						if(pa.distance(pb) < Math.abs(mind)){
							mind = pa.offset(pb);
						nearestPeak = pb;
						}
					}
				}
				if(Math.abs(mind) <= overlapD){
					System.out.println(pa.getLocationString()+"\t"+nearestPeak.getLocationString()+"\t"+Integer.toString(mind));
				}else{
					System.out.println(pa.getLocationString()+"\tNA");
				}
			}

		}

	}
	
	
	public void printClosestPeaks(){
		HashMap<String, List<Point>> peaksAbyChrom = hashbychrom(peaksA);
		HashMap<String, List<Point>> peaksBbyChrom = hashbychrom(peaksB);
		for(String chrom : peaksAbyChrom.keySet()){
			if(!peaksBbyChrom.containsKey(chrom)){
				for(Point pa : peaksAbyChrom.get(chrom)){
					System.out.println(pa.getLocationString()+"\t-");
				}
				continue;
			}
				
			for(Point pa : peaksAbyChrom.get(chrom)){
				int mind=Integer.MAX_VALUE;
				Point nearestPeak=null;
				for(Point pb : peaksBbyChrom.get(chrom)){
					if(pb.getChrom().equals(pa.getChrom())){
						if(pa.distance(pb) < mind){
							mind = pa.distance(pb);
							nearestPeak = pb;
						}
					}
				}
				if(mind <= overlapD){
					System.out.println(pa.getLocationString()+"\t"+nearestPeak.getLocationString()+"\t"+Integer.toString(mind));
				}else{
					System.out.println(pa.getLocationString()+"\t-");
				}
			
			}
		}
		
	}
	
	public void printPeakSets(){
		HashMap<String, List<Point>> peaksAbyChrom = hashbychrom(peaksA);
		List<HashMap<String,List<Point>>> peaksBSetbyChrom  = new ArrayList<HashMap<String,List<Point>>>();
		for(List<Point> set : peaksBSet){
			peaksBSetbyChrom.add(hashbychrom(set));
		}
		
		for(String chrom : peaksAbyChrom.keySet()){
			for(Point pa: peaksAbyChrom.get(chrom)){
				StringBuilder sb= new StringBuilder();
				sb.append(pa.getLocationString());sb.append("\t");
				Point[] nearestPeaks = new Point[peaksBSet.size()];
				int[] nearestDis = new int[peaksBSet.size()];
				Arrays.fill(nearestDis, Integer.MAX_VALUE);
				int setId=0;
				for(HashMap<String,List<Point>> bset: peaksBSetbyChrom){
					if(bset.containsKey(chrom)){
						for(Point pb : bset.get(chrom)){
							if(pb.getChrom().equals(pa.getChrom())){
								if(pa.distance(pb) < nearestDis[setId]){
									nearestDis[setId] = pa.distance(pb);
									nearestPeaks[setId] =pb;
								}
							}
						}
					}
					setId++;
				}
				for(int i=0; i<nearestDis.length;i++){
					if(nearestDis[i]<=overlapD){
						sb.append(nearestPeaks[i].getLocationString());sb.append("\t");
					}else{
						sb.append("-");sb.append("\t");
					}
				}
				sb.deleteCharAt(sb.length()-1);
				System.out.println(sb.toString());
			}
		}
	}

	public HashMap<String, List<Point>> hashbychrom(List<Point> pts){
		HashMap<String, List<Point>> byChr = new HashMap<String, List<Point>>();
		for(Point p : pts){
			if(!byChr.containsKey(p.getChrom()))
				byChr.put(p.getChrom(), new ArrayList<Point>());
			byChr.get(p.getChrom()).add(p);
		}
		return byChr;
	}
	
	
	public static void main(String[] args){
		GenomeConfig gconfig = new GenomeConfig(args);
		ArgParser ap = new ArgParser(args);
		
		int win = Args.parseInteger(args, "win", 150);
		
		PeaksVsPeaks analyzer = new PeaksVsPeaks(gconfig.getGenome());
		
		String peaksAfile = Args.parseString(args, "peaksA", null);
		if(peaksAfile == null){
			System.err.println("Provide ChIP-Seq peaks file A!!");
			return;
		}
		List<Point> peaksA = RegionFileUtilities.loadPeaksFromPeakFile(gconfig.getGenome(), peaksAfile, win);
		List<Region> regionsA = RegionFileUtilities.loadRegionsFromPeakFile(gconfig.getGenome(), peaksAfile, win);
		analyzer.setPeaksA(peaksA);
		analyzer.setRegsA(regionsA);
		
		if(!ap.hasKey("peaksB") && !ap.hasKey("peaksBset") && !ap.hasKey("uniq") && !ap.hasKey("closebyCount")){
			System.err.println("Provide ChIP-Seq peaks file to compare with the A list!!");
			return;
		}
		
		String peaksBfile = null;
		if(ap.hasKey("peaksB")){
			peaksBfile = Args.parseString(args, "peaksB", null);
			List<Point> peaksB = RegionFileUtilities.loadPeaksFromPeakFile(gconfig.getGenome(), peaksBfile, win);
			List<Region> regionsB = RegionFileUtilities.loadRegionsFromPeakFile(gconfig.getGenome(), peaksBfile, win);
			analyzer.setPeaksB(peaksB);
			analyzer.setRegsB(regionsB);
		}
		
		Collection<String> peaksBfileList = null;
		if(ap.hasKey("peaksBset")){
			peaksBfileList =  Args.parseStrings(args, "peaksBset");
			List<List<Point>> bPeaksSet = new ArrayList<List<Point>>();
			List<List<Region>> bRegionSet = new ArrayList<List<Region>>();
			for(String s : peaksBfileList){
				List<Point> currPoints = RegionFileUtilities.loadPeaksFromPeakFile(gconfig.getGenome(), s, win);
				bPeaksSet.add(currPoints);
				List<Region> currRegions = RegionFileUtilities.loadRegionsFromPeakFile(gconfig.getGenome(), s, win);
				bRegionSet.add(currRegions);
			}
			analyzer.setPeaksBSet(bPeaksSet);
			analyzer.setRegionsBSet(bRegionSet);
		}
		
		
		int overlapD = Args.parseInteger(args, "overlapD", 100);
		analyzer.setOverlapD(overlapD);
		
		if(ap.hasKey("peaksB") && !ap.hasKey("printoffset"))
			analyzer.printClosestPeaks();
		if(ap.hasKey("peaksBset"))
			analyzer.printPeakSets();
		if(ap.hasKey("peaksB") && ap.hasKey("printoffset"))
			analyzer.printClosestOffsetPeak();
		if(!ap.hasKey("peaksB") && !ap.hasKey("peaksBset") && ap.hasKey("uniq")){
			int rad = Args.parseInteger(args, "rad", 50);
			analyzer.printUnique(rad);
		}
		if(!ap.hasKey("peaksB") && !ap.hasKey("peaksBset") && ap.hasKey("closebyCount")){
			int rad = Args.parseInteger(args, "rad", 50);
			analyzer.printClosePeaksInAList(rad);
		}
	}
	
	
	
	

}
