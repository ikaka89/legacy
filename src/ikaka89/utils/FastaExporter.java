package ikaka89.utils;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.gseutils.NotFoundException;
import org.seqcode.gseutils.Pair;
import org.seqcode.viz.metaprofile.EventMetaMaker;


public class FastaExporter {
	
	
	public List<Region> regs;
	public Genome gen;
	public SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();;
	public List<String> fasta = new ArrayList<String>();
	public List<String> names = new ArrayList<String>();
	public List<String> strand = new ArrayList<String>();
	
	
	public FastaExporter(int win, List<Region> rs, Genome gen) {
		regs = new ArrayList<Region>();
		regs.addAll(rs);
		this.gen = gen;
		
	}
	
	public void execute(boolean cache, String seqPathName){
		if(cache){
			seqgen.useCache(cache);
			seqgen.setGenomePath(seqPathName);
		}
		
		for(int rind=0; rind<regs.size(); rind++){
			String fa = seqgen.execute(regs.get(rind)).toUpperCase();
			if(strand.size()>0 && strand.get(rind).equals("-")){
				fasta.add(SequenceUtils.reverseComplement(fa));
			}else{
				fasta.add(fa);
			}
		}
		int count=0;
		for(int i=0; i< regs.size(); i++){
			if(names.size()>0){
				System.out.println(">"+names.get(i));
			}else{
				System.out.println(">seq"+Integer.toString(count));
			}
			System.out.println(fasta.get(i));
			count++;
		}
	}
	
	
	// Mutators
	public void setRegions(List<Region> regions){regs = regions;}
	public void setNames(List<String> ns){names.addAll(ns);}
	public void setStrands(List<String>sts){strand.addAll(sts);}
	
	public static void main(String[] args) throws IOException{
		try{
			ArgParser ap =  new ArgParser(args);
			int window = Args.parseInteger(args, "win", -1);
			boolean cache = ap.hasKey("cache");
			String seqPathName = "";
			if(cache){
				seqPathName = Args.parseString(args, "seq", "");
			}
			
			if(!ap.hasKey("locations") && !ap.hasKey("regions")){
				System.err.println("Provide either a location or a regions file !!");
				System.exit(1);
			}
			
			Genome gen;
			if(ap.hasKey("species") || ap.hasKey("genome") || ap.hasKey("gen")){
				Pair<Species, Genome> orggen = Args.parseGenome(args);
				gen = orggen.cdr();
			} 
			else{
				if(ap.hasKey("geneinfo") || ap.hasKey("g")){
					String infofilename = ap.hasKey("geneinfo") ? ap.getKeyValue("geneinfo") : ap.getKeyValue("g"); 
					gen = new Genome("Genome",new File(infofilename),true);
				}
				else{
					gen = null;
				}
			}
			
			boolean hasnames = ap.hasKey("hasnames");
			boolean isStranded = ap.hasKey("stranded");
			List<String> ns = new ArrayList<String>();
			List<String> strs = new ArrayList<String>();
			if(hasnames){
				BufferedReader br = new BufferedReader(new FileReader(ap.getKeyValue("locations")));
				String line = null;
				while((line=br.readLine()) != null){
					String[] pieces = line.split("\t");
					ns.add(pieces[1]);
					if(isStranded){
						String[] subpieces = pieces[0].split(":");
						strs.add(subpieces[2]);
					}
				}
				br.close();
			}
			List<Region> search_locs = new ArrayList<Region>();
			if(isStranded){
				List<StrandedRegion> srs  = new ArrayList<StrandedRegion>();
				srs.addAll(RegionFileUtilities.loadStrandedRegionsFromMotifFile(gen, ap.getKeyValue("locations"), window));
				for(StrandedRegion sr:  srs){
					search_locs.add(new Region(gen, sr.getChrom(), sr.getStart(), sr.getEnd()));
				}
			}else{
				search_locs.addAll(RegionFileUtilities.loadRegionsFromPeakFile(gen, ap.getKeyValue("locations"), window));
			}
			
			FastaExporter exporter = new FastaExporter(window, search_locs, gen);
			if(hasnames)
				exporter.setNames(ns);
			if(isStranded)
				exporter.setStrands(strs);
			exporter.execute(cache, seqPathName);
		}catch(NotFoundException e){
			e.printStackTrace();
		} 
	}

}
