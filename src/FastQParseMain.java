import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Date;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.StreamSupport;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class FastQParseMain {
	private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("###,###.#######"); //standard number format
	private static final DateFormat DATE_FORMAT = new SimpleDateFormat("EEE MMM dd HH:mm:ss zzz yyyy"); //standard date format
	private static final int BUFFER_SIZE = 1048576; //buffer size for buffered reader/writer
	private static final int BUFFER_SIZE_GZIP = 1048576; //buffer size for gzip stream
	private static final int BUFFER_SIZE_LOG = 16384; //buffer size for log/stats files
	private static final String description2 = "+"; //the second description (3rd line in each read)
	
	private static PrintWriter logWriter; //prints to log file
	
	private static File inputFile; //.fastq file with all data
	private static File inputFile2; //.fastq file with data for paired-ends sequencing
	private static File sampleFile; //.txt file with sample DNA sequence
	private static File indexFile; //index file for UMI
	private static File indexFile2; //index file for paired-ends sequencing
	private static String outputDir; //directory where output files go
	
	private static int randUMILength; //length of random UMI sequence (could be 0)
	
	private static boolean removeFirstDup = false; //remove duplicates that are not the first
	private static boolean removeBestDup = false; //remove duplicates that are not the best
	
	private static boolean replaceOriginal = false; //replace original file
	
	private static int maxOffsetB = 0; //tries offsets from 0 to maxOffset when matching barcode/enzyme
	private static int maxOffsetA = 0; //tries offsets from 0 to maxOffset when matching adapter
	
	private static double removeDNAWithNPercent = 2.0; //if the percentage of DNA that is N is greater than this value then it will be undetermined
	
	private static double qualityFilter = 0.0; //quality filter threshold
	private static boolean filterAlgorithm = false; //false = average, true = error sum
	
	private static boolean inputGZIP = false; //is input file gzipped (sample file is not gzipped)
	private static boolean outputGZIP = false; //is output file gzipped
	
	private static boolean saveTemp = false; //save to temp file or not
	
	private static boolean removeBarRand = true; //remove barcode and random from beginning of DNA (removes quality too)
	private static boolean removeEnzyme = true; //remove enzyme from beginning of DNA (removes quality too)
	
	private static double editMaxB = 1.0; //how many edits to allow (0 is exact match) for barcodes
	private static double editMaxA = 1.0; //how many edits to allow (0 is exact match) for adapters
	private static double editMaxM = 1.0; //how many edits to allow (0 is exact match) for merging paired-ends
	
	private static boolean saveDup = false; //save duplicates in separate files
	
	private static long printProcessedInterval = 80000000L; //interval to print the number of DNA processed while running
	private static long printDuplicateInterval = 5000000L; //interval to print number of duplicates removed
	
	private static boolean mergePairedEnds = false; //whether to merge corresponding paired ends reads
	
	private static int minLength = 0; //every DNA sequence with total length < this will be undetermined
	private static int maxLength = Integer.MAX_VALUE; //every DNA sequence with total length > this will be undetermined
	
	private static ArrayList<Adapter> adaptersF = new ArrayList<Adapter>(); //adapters to remove from forwards reads (file 1)
	private static ArrayList<Adapter> adaptersR = new ArrayList<Adapter>(); //adapters to remove from reversed reads (file 2)
	
	private static int minOverlapA = Integer.MAX_VALUE; //minimum overlap when matching adapters
	private static int minOverlapB = Integer.MAX_VALUE; //minimum overlap when matching barcodes
	
	private static boolean trimAlgorithm = false; //false = local average, true = sum
	private static int qualityTrimQScore1 = 0; //quality under this will be trimmed in quality trimming (from start)
	private static int qualityTrimQScore2 = 0; //quality under this will be trimmed in quality trimming (from end)
	private static int qualityTrimLength = 1; //length that is needed for quality trimming
	
	private static double trimNPercent = 2.0; //trim leading and trailing 'N' percentage
	
	private static boolean allowIndelsB = false; //allow insertions and deletions for barcode/enzyme
	private static boolean allowIndelsA = false; //allow insertions and deletions for adapters
	
	private static boolean checkReversedReads = false; //whether to check for barcodes/enzymes in reversed reads
	private static boolean wildcard = false; //whether to check for undetermined bp 'N'
	private static boolean singleBarcodeMatchOnly = false; //whether to only allow one barcode match
	
	private static boolean removeUntrimmedReads = false;
	private static boolean removeNoAdapterReads = false;
	private static boolean removeNoMergeReads = false;
	
	private static boolean parallel = false; //use parallel streams or not
	private static int splitBatchSize = 2048; //batch size when using ReadSpliterator
	
	private static boolean simReversed = false; //generate simulated data
	private static boolean simUMI = false;
	private static boolean simMerging = false;
	private static int simReadLength = 100;
	private static long simIter = 1000;
	
	private static double probB = -1.0; //prior probabilities for probability based matching
	private static double probA = -1.0;
	private static double probM = -1.0;
	
	private enum Mode{ //features
		DEMULTIPLEX("Demultiplex", "demultiplex"), DEDUP("Deduplicate Reads", "dedup"),
		PAIRMERGE("Merge Paired-End Reads", "pairmerge"), FILTER("Filter Reads", "filter"),
		SIM_READS("Generate Simulated Reads", "sim");
		
		String description1;
		String description2;
		
		private Mode(String description1, String description2){
			this.description1 = description1;
			this.description2 = description2;
		}
	};
	private enum Stat{ //stats (description used as column titles in .stats file)
		SEQUENCE_FORWARDS("Forwards Barcodes"), SEQUENCE_REVERSED("Reversed Barcodes"), SEQUENCE_COUNT("Reads"),
		SEQUENCE_PERCENT("% of Total Reads"), MERGED_COUNT("Merged Count"), MERGED_PERCENT("% Merged"),
		DEDUP_COUNT("Deduplicated Count"), DEDUP_PERCENT("% Deduplicated"), REMOVEADAPTER_COUNT("Removed Adapter Count"),
		REMOVEADAPTER_PERCENT("% Removed Adapter"), QUALITYTRIM_COUNT("Quality Trimmed Count"), QUALITYTRIM_PERCENT("% Quality Trimmed"),
		BASE_COUNT("BP Count"), BASE_PERCENT("% of Total BP Count");
		
		String description;
		
		private Stat(String description){
			this.description = description;
		}
	};
	
	private HashMap<String, String> sampleMapF = new HashMap<String, String>();
	private ArrayList<String> sampleDNAF = new ArrayList<String>(); //names saved twice to keep them in order
	private ArrayList<String> sampleDNAR = new ArrayList<String>();
	
	private ArrayList<String> constEnzymesF; //constant enzyme sequence
	private ArrayList<String> constEnzymesR;
	
	private ArrayList<HashMap<Character, BitVector>> barcodePatternF;
	private ArrayList<HashMap<Character, BitVector>> barcodePatternR;
	private ArrayList<HashMap<Character, BitVector>> enzymePatternF;
	private ArrayList<HashMap<Character, BitVector>> enzymePatternR;
	private ArrayList<HashMap<Character, BitVector>> adapterPatternF;
	private ArrayList<HashMap<Character, BitVector>> adapterPatternR;
	
	private LongAdder totalDNAProcessed = new LongAdder(); //total DNA reads
	private long duplicatesRemoved = 0L; //total duplicates removed
	private LongAdder undeterminedDNA = new LongAdder(); //total undetermined DNA
	private LongAdder totalRemovedAdapters = new LongAdder(); //total DNA with adapters removed
	private LongAdder totalQualityTrimmed = new LongAdder(); //total DNA that has been quality trimmed
	private LongAdder baseCount = new LongAdder(); //total base count
	private LongAdder totalReadsMerged = new LongAdder(); //total number of reads that are actually merged
	
	private String enzymeF; //sample file enzyme name
	private String enzymeR;
	
	private boolean hasReversedBarcode = false;
	
	private int generatedDupFiles = 0; //total generated duplicate files
	private boolean generatedUndeterminedFile = false;
	
	private long startTime;
	
	private HashMap<String, EnumMap<Stat, String>> stats = new HashMap<String, EnumMap<Stat, String>>(); //stats for each sample
	
	//constructor handles printing to log and calling the methods to process the files
	public FastQParseMain(Mode mode) throws Exception{
		logWriter.println("Mode: " + mode.description1);
		System.out.println("Mode: " + mode.description1);
		logWriter.println();
		
		if(mode == Mode.DEMULTIPLEX){
			logWriter.println("Sample File: " + sampleFile.getAbsolutePath());
			logWriter.println("Processing Sample File...");
			logWriter.println();
			logWriter.flush();
			readSample(); //read sample file
			logWriter.println("Sample File Processing Finished.");
			logWriter.println("Number of Samples: " + sampleDNAF.size());
			logWriter.println("Has Reversed Barcodes: " + hasReversedBarcode);
			logWriter.println("Forwards Enzyme Name: " + enzymeF);
			String constEnzymeFString = constEnzymesF.toString();
			logWriter.println("Forwards Enzyme Sequence: " + constEnzymeFString.substring(1, constEnzymeFString.length() - 1));
			if(inputFile2 != null){
				logWriter.println("Reversed Enzyme Name: " + enzymeR);
				String constEnzymeRString = constEnzymesR.toString();
				logWriter.println("Reversed Enzyme Sequence: " + constEnzymeRString.substring(1, constEnzymeRString.length() - 1));
			}
			logWriter.println();
			
			logWriter.println("Forwards Read File: " + inputFile.getAbsolutePath());
			if(indexFile != null)
				logWriter.println("Forwards Index File: " + indexFile.getAbsolutePath());
			if(inputFile2 != null){
				logWriter.println("Reversed Read File: " + inputFile2.getAbsolutePath());
				if(indexFile2 != null)
					logWriter.println("Reversed Index File: " + indexFile2.getAbsolutePath());
			}
			logWriter.println("Output Directory: " + outputDir);
			logWriter.println("Length of Random UMI: " + randUMILength);
			logWriter.println("Maximum Right Offset For Barcodes/Enzymes: " + maxOffsetB);
			logWriter.println("Maximum Right Offset For Adapters: " + maxOffsetA);
			if(removeDNAWithNPercent >= 0 && removeDNAWithNPercent <= 1.0)
				logWriter.println("Remove Reads With % of N Greater Than: " + removeDNAWithNPercent);
			else if(removeDNAWithNPercent < 0)
				logWriter.println("Remove Reads With at Least 1 N");
			else
				logWriter.println("Remove Reads With N: false");
			logWriter.println("Quality Filter Threshold: " + qualityFilter);
			logWriter.println("Quality Filter Algorithm: " + (filterAlgorithm ? "Error Sum" : "Average"));
			logWriter.println("Only Keep First Duplicate: " + removeFirstDup);
			logWriter.println("Only Keep Best Duplicate: " + removeBestDup);
			logWriter.println("Is Input GZIP Format: " + inputGZIP);
			logWriter.println("Is Output GZIP Format: " + outputGZIP);
			logWriter.println("Save Temporary Files: " + saveTemp);
			logWriter.println("Remove Barcode and Random: " + removeBarRand);
			logWriter.println("Remove Enzyme: " + removeEnzyme);
			if(editMaxB < 0.0)
				logWriter.println("Edit Percentage to Count as Match For Barcodes/Enzymes: " + -editMaxB);
			else
				logWriter.println("Max Edits to Count as Match For Barcodes/Enzymes: " + editMaxB);
			if(editMaxA < 0.0)
				logWriter.println("Edit Percentage to Count as Match For Adapters: " + -editMaxA);
			else
				logWriter.println("Max Edits to Count as Match For Adapters: " + editMaxA);
			if(editMaxM < 0.0)
				logWriter.println("Edit Percentage to Count as Match For Merging Paired-Ends: " + -editMaxM);
			else
				logWriter.println("Max Edits to Count as Match For Merging Paired-Ends: " + editMaxM);
			logWriter.println("Allow Insertions and Deletions for Barcodes/Enzymes: " + allowIndelsB);
			logWriter.println("Allow Insertions and Deletions for Adapters: " + allowIndelsA);
			logWriter.println("Save Duplicates in Separate Files: " + saveDup);
			logWriter.println("Interval to Print Processed Reads: " + DECIMAL_FORMAT.format(printProcessedInterval));
			logWriter.println("Interval to Print Duplicate Reads: " + DECIMAL_FORMAT.format(printDuplicateInterval));
			logWriter.println("Minimum Read Length: " + minLength);
			logWriter.println("Maximum Read Length: " + maxLength);
			logWriter.println("Merge Paired-End Reads: " + mergePairedEnds);
			String adaptersFString = adaptersF.toString();
			logWriter.println("Forwards Read Adapters: " + adaptersFString.substring(1, adaptersFString.length() - 1));
			if(inputFile2 != null){
				String adaptersRString = adaptersR.toString();
				logWriter.println("Reversed Read Adapters: " + adaptersRString.substring(1, adaptersRString.length() - 1));
			}
			logWriter.println("Minimum Adapter Overlap: " + minOverlapA);
			logWriter.println("Minimum Barcode Overlap: " + minOverlapB);
			logWriter.println("Quality Trim Algorithm: " + (trimAlgorithm ? "Sum" : "Local Average"));
			logWriter.println("5' Quality Trim Score Threshold: " + qualityTrimQScore1);
			logWriter.println("3' Quality Trim Score Threshold: " + qualityTrimQScore2);
			logWriter.println("Quality Trim Length: " + qualityTrimLength);
			if(trimNPercent > 1.0)
				logWriter.println("Trim Leading and Trailing N: false");
			else
				logWriter.println("Percentage to Trim Leading and Trailing N By: " + trimNPercent);
			logWriter.println("Check Reversed Reads for Enzyme or Barcode: " + checkReversedReads);
			logWriter.println("Use Wildcard Characters: " + wildcard);
			logWriter.println("Remove Reads With Multiple Barcode Matches: " + singleBarcodeMatchOnly);
			logWriter.println("Remove Reads That Are Not Quality Trimmed: " + removeUntrimmedReads);
			logWriter.println("Remove Reads That Do Not Contain Any Adapters: " + removeNoAdapterReads);
			logWriter.println("Remove Reads That Are Not Merged: " + removeNoMergeReads);
			logWriter.println("Parallel: " + parallel);
			logWriter.println("Parallel Batch Size: " + splitBatchSize);
			logWriter.println("Useable Processors: " + Runtime.getRuntime().availableProcessors());
			if(probB < 0.0)
				logWriter.println("Probability Based Matching for Barcodes: false");
			else
				logWriter.println("Prior Probability for Barcode Matching: " + probB);
			if(probA < 0.0)
				logWriter.println("Probability Based Matching for Adapters: false");
			else
				logWriter.println("Prior Probability for Adapter Matching: " + probA);
			if(probM < 0.0)
				logWriter.println("Probability Based Matching for Merging Paired End Reads: false");
			else
				logWriter.println("Prior Probability for Paired End Merging: " + probM);
			logWriter.println();
			logWriter.println("Demultiplexing...");
			logWriter.println();
			logWriter.flush();
			
			//initialize the stats map
			EnumMap<Stat, String> map;
			for(int i = 0; i < sampleDNAF.size(); i++){
				map = new EnumMap<Stat, String>(Stat.class);
				map.put(Stat.SEQUENCE_FORWARDS, sampleDNAF.get(i));
				if(hasReversedBarcode)
					map.put(Stat.SEQUENCE_REVERSED, sampleDNAR.get(i));
				map.put(Stat.SEQUENCE_COUNT, "");
				map.put(Stat.SEQUENCE_PERCENT, "");
				map.put(Stat.MERGED_COUNT, "");
				map.put(Stat.MERGED_PERCENT, "");
				map.put(Stat.DEDUP_COUNT, DECIMAL_FORMAT.format(0));
				map.put(Stat.DEDUP_PERCENT, DECIMAL_FORMAT.format(0));
				map.put(Stat.REMOVEADAPTER_COUNT, "");
				map.put(Stat.REMOVEADAPTER_PERCENT, "");
				map.put(Stat.QUALITYTRIM_COUNT, "");
				map.put(Stat.QUALITYTRIM_PERCENT, "");
				map.put(Stat.BASE_COUNT, "");
				map.put(Stat.BASE_PERCENT, "");
				stats.put(sampleMapF.get(sampleDNAF.get(i)), map);
			}
			map = new EnumMap<Stat, String>(Stat.class);
			map.put(Stat.SEQUENCE_FORWARDS, "");
			if(hasReversedBarcode)
				map.put(Stat.SEQUENCE_REVERSED, "");
			map.put(Stat.SEQUENCE_COUNT, "");
			map.put(Stat.SEQUENCE_PERCENT, "");
			map.put(Stat.MERGED_COUNT, DECIMAL_FORMAT.format(0));
			map.put(Stat.MERGED_PERCENT, DECIMAL_FORMAT.format(0));
			map.put(Stat.DEDUP_COUNT, DECIMAL_FORMAT.format(0));
			map.put(Stat.DEDUP_PERCENT, DECIMAL_FORMAT.format(0));
			map.put(Stat.REMOVEADAPTER_COUNT, DECIMAL_FORMAT.format(0));
			map.put(Stat.REMOVEADAPTER_PERCENT, DECIMAL_FORMAT.format(0));
			map.put(Stat.QUALITYTRIM_COUNT, DECIMAL_FORMAT.format(0));
			map.put(Stat.QUALITYTRIM_PERCENT, DECIMAL_FORMAT.format(0));
			map.put(Stat.BASE_COUNT, "");
			map.put(Stat.BASE_PERCENT, "");
			stats.put("Undetermined", map);
			
			if(probB < 0.0){
				barcodePatternF = new ArrayList<HashMap<Character, BitVector>>();
				for(int i = 0; i < sampleDNAF.size(); i++){
					barcodePatternF.add(UtilMethods.genPatternMasks(sampleDNAF.get(i), allowIndelsB, wildcard));
				}
				
				enzymePatternF = new ArrayList<HashMap<Character, BitVector>>();
				for(int i = 0; i < constEnzymesF.size(); i++){
					enzymePatternF.add(UtilMethods.genPatternMasks(constEnzymesF.get(i), allowIndelsB, wildcard));
				}
				
				if(inputFile2 != null && checkReversedReads){
					barcodePatternR = new ArrayList<HashMap<Character, BitVector>>();
					for(int i = 0; i < sampleDNAR.size(); i++){
						barcodePatternR.add(UtilMethods.genPatternMasks(sampleDNAR.get(i), allowIndelsB, wildcard));
					}
					
					enzymePatternR = new ArrayList<HashMap<Character, BitVector>>();
					for(int i = 0; i < constEnzymesR.size(); i++){
						enzymePatternR.add(UtilMethods.genPatternMasks(constEnzymesR.get(i), allowIndelsB, wildcard));
					}
				}
			}
			
			if(probA < 0.0){
				adapterPatternF = new ArrayList<HashMap<Character, BitVector>>();
				for(int i = 0; i < adaptersF.size(); i++){
					adapterPatternF.add(UtilMethods.genPatternMasks(adaptersF.get(i).isStart ? adaptersF.get(i).str : UtilMethods.reverse(adaptersF.get(i).str), allowIndelsA, wildcard));
				}
				
				if(inputFile2 != null){
					adapterPatternR = new ArrayList<HashMap<Character, BitVector>>();
					for(int i = 0; i < adaptersR.size(); i++){
						adapterPatternR.add(UtilMethods.genPatternMasks(adaptersR.get(i).isStart ? adaptersR.get(i).str : UtilMethods.reverse(adaptersR.get(i).str), allowIndelsA, wildcard));
					}
				}
			}
			
			startTime = System.currentTimeMillis();
			int generatedSampleFiles = 0;
			
			//deduplicate after demultiplex
			if(removeFirstDup || removeBestDup){ //assume that index files are provided for deduplicating for simpler logic
				ConcurrentLinkedQueue<Strings> files = demultiplexFile(); //demultiplex
				logWriter.println();
				logWriter.println("Demultiplex Completed.");
				logWriter.println();
				logWriter.println("Run Time So Far: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
				logWriter.println();
				logWriter.println("Deduplicating...");
				logWriter.flush();
				for(Strings str : files){
					File file1 = new File(str.get(0));
					int pos = file1.getName().indexOf(".");
					String justName1 = pos > 0 ? file1.getName().substring(0, pos) : file1.getName();
					
					File file2 = new File(str.get(1));
					pos = file2.getName().indexOf(".");
					String justName2 = pos > 0 ? file2.getName().substring(0, pos) : file2.getName();
					
					logWriter.println();
					logWriter.println("Deduplicating File: " + justName1);
					logWriter.flush();
					//deduplicate
					long removed = deduplicate(str.get(0), inputFile2 == null ? null : str.get(2), str.get(1), inputFile2 == null ? null : str.get(3),
							outputDir, false, false, outputDir + "dup" + File.separatorChar + justName1 + "_dup.fastq.gz", outputDir + "dup" + File.separatorChar + justName2 + "_dup.fastq.gz");
					for(int j = 0; j < sampleDNAF.size(); j++){ //calculate duduplicated reads count
						if(justName1.contains(sampleMapF.get(sampleDNAF.get(j)))){
							map = stats.get(sampleMapF.get(sampleDNAF.get(j)));
							map.put(Stat.DEDUP_COUNT, DECIMAL_FORMAT.format(removed));
							map.put(Stat.DEDUP_PERCENT, DECIMAL_FORMAT.format((double)removed / Double.parseDouble(map.get(Stat.SEQUENCE_COUNT).replace(",", ""))));
							break;
						}
					}
					generatedSampleFiles += str.size();
					System.gc();
				}
				if(generatedUndeterminedFile)
					generatedSampleFiles += (inputFile2 == null ? 1 : 2) + (indexFile == null ? 0 : (inputFile2 == null ? 1 : 2));
				logWriter.println();
				logWriter.println("Deduplicating Completed.");
				if(!saveTemp){
					File tempFile = new File(outputDir + "temp" + File.separatorChar);
					if(tempFile.exists())
						UtilMethods.deleteFolder(tempFile);
				}
			}else{ //no deduplicating after demultiplex
				ConcurrentLinkedQueue<Strings> files = demultiplexFile(); //demultiplex
				for(Strings str : files)
					generatedSampleFiles += str.size();
				if(generatedUndeterminedFile)
					generatedSampleFiles += (inputFile2 == null ? 1 : 2) + (indexFile == null ? 0 : (inputFile2 == null ? 1 : 2));
				logWriter.println();
				logWriter.println("Demultiplex Completed.");
				File tempFile = new File(outputDir + "temp" + File.separatorChar);
				if(tempFile.exists())
					UtilMethods.deleteFolder(tempFile);
			}
			
			logWriter.println();
			logWriter.println("Total Run Time: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
			logWriter.println();
			logWriter.println("Number of Lines Processed: " + DECIMAL_FORMAT.format(totalDNAProcessed.sum() * 4));
			logWriter.println("Number of Reads Processed: " + DECIMAL_FORMAT.format(totalDNAProcessed.sum()));
			logWriter.println("Number of Generated Sample Files: " + generatedSampleFiles);
			logWriter.println("Number of Generated Duplicate Files: " + generatedDupFiles);
			
			//fill in info for total stats
			map = new EnumMap<Stat, String>(Stat.class);
			map.put(Stat.SEQUENCE_FORWARDS, "");
			if(hasReversedBarcode)
				map.put(Stat.SEQUENCE_REVERSED, "");
			map.put(Stat.SEQUENCE_COUNT, DECIMAL_FORMAT.format(totalDNAProcessed.sum()));
			map.put(Stat.SEQUENCE_PERCENT, DECIMAL_FORMAT.format(1));
			map.put(Stat.MERGED_COUNT, DECIMAL_FORMAT.format(totalReadsMerged.sum()));
			map.put(Stat.MERGED_PERCENT, DECIMAL_FORMAT.format((double)totalReadsMerged.sum() / (double)totalDNAProcessed.sum()));
			map.put(Stat.DEDUP_COUNT, DECIMAL_FORMAT.format(duplicatesRemoved));
			map.put(Stat.DEDUP_PERCENT, DECIMAL_FORMAT.format((double)duplicatesRemoved / (double)totalDNAProcessed.sum()));
			map.put(Stat.REMOVEADAPTER_COUNT, DECIMAL_FORMAT.format(totalRemovedAdapters.sum()));
			map.put(Stat.REMOVEADAPTER_PERCENT, DECIMAL_FORMAT.format((double)totalRemovedAdapters.sum() / (double)totalDNAProcessed.sum()));
			map.put(Stat.QUALITYTRIM_COUNT, DECIMAL_FORMAT.format(totalQualityTrimmed.sum()));
			map.put(Stat.QUALITYTRIM_PERCENT, DECIMAL_FORMAT.format((double)totalQualityTrimmed.sum() / (double)totalDNAProcessed.sum()));
			map.put(Stat.BASE_COUNT, DECIMAL_FORMAT.format(baseCount.sum()));
			map.put(Stat.BASE_PERCENT, DECIMAL_FORMAT.format(1));
			stats.put("Total", map);
			
			//print stats file with cool padding
			int width = 25;
			int space = 1;
			
			PrintWriter statsWriter = new PrintWriter(new BufferedWriter(new FileWriter(outputDir + "FastQParse_" + mode.description2 + ".stats"), BUFFER_SIZE_LOG));
			StringBuilder builder = new StringBuilder();
			builder.append("Sample");
			for(int i = 0; i < width - "Sample".length(); i++){
				builder.append(' ');
			}
			for(Stat column : map.keySet()){
				for(int i = 0; i < space; i++){
					builder.append(' ');
				}
				builder.append(column.description);
				for(int i = 0; i < width - column.description.length(); i++){
					builder.append(' ');
				}
			}
			statsWriter.println(builder.toString());
			for(int i = 0; i < sampleDNAF.size() + 2; i++){
				String row = i == sampleDNAF.size() ? "Undetermined" : (i == sampleDNAF.size() + 1 ? "Total" : sampleMapF.get(sampleDNAF.get(i)));
				builder = new StringBuilder();
				builder.append(row);
				for(int j = 0; j < width - row.length(); j++){
					builder.append(' ');
				}
				for(String item : stats.get(row).values()){
					for(int j = 0; j < space; j++){
						builder.append(' ');
					}
					builder.append(item);
					for(int j = 0; j < width - item.length(); j++){
						builder.append(' ');
					}
				}
				statsWriter.println(builder.toString());
			}
			statsWriter.close();
		}else if(mode == Mode.DEDUP){
			logWriter.println("Input File: " + inputFile.getAbsolutePath());
			logWriter.println("Index File: " + indexFile.getAbsolutePath());
			logWriter.println("Replace Original File: " + replaceOriginal);
			logWriter.println("Output Directory: " + outputDir);
			logWriter.println("Only Keep First Duplicate: " + removeFirstDup);
			logWriter.println("Only Keep Best Duplicate: " + removeBestDup);
			logWriter.println("Length of Random UMI: " + randUMILength);
			logWriter.println("Is Input GZIP Format: " + inputGZIP);
			logWriter.println("Is Output GZIP Format: " + outputGZIP);
			logWriter.println("Save Duplicates in a Separate File: " + saveDup);
			logWriter.println("Interval to Print Processed Reads: " + DECIMAL_FORMAT.format(printProcessedInterval));
			logWriter.println("Interval to Print Duplicate Reads: " + DECIMAL_FORMAT.format(printDuplicateInterval));
			logWriter.println();
			logWriter.println("Deduplicating...");
			logWriter.println();
			logWriter.flush();
			
			startTime = System.currentTimeMillis();
			
			int pos = inputFile.getName().indexOf(".");
			String justName1 = pos > 0 ? inputFile.getName().substring(0, pos) : inputFile.getName();
			
			pos = indexFile.getName().indexOf(".");
			String justName2 = pos > 0 ? indexFile.getName().substring(0, pos) : indexFile.getName();
			
			//deduplicate
			deduplicate(inputFile.getAbsolutePath(), null, indexFile.getAbsolutePath(), null, outputDir, replaceOriginal, true,
					replaceOriginal ? (inputFile.getParentFile().getAbsolutePath() + File.separatorChar + justName1 + "_dup.fastq.gz") : (outputDir + justName1 + "_dup.fastq.gz"),
							replaceOriginal ? (indexFile.getParentFile().getAbsolutePath() + File.separatorChar + justName2 + "_dup.fastq.gz") : (outputDir + justName2 + "_dup.fastq.gz"));
			
			logWriter.println();
			logWriter.println("Deduplication Completed.");
			logWriter.println();
			logWriter.println("Total Run Time: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
			logWriter.println();
			logWriter.println("Number of Lines Processed: " + DECIMAL_FORMAT.format(totalDNAProcessed.sum() * 4));
			logWriter.println("Number of Reads Processed: " + DECIMAL_FORMAT.format(totalDNAProcessed.sum()));
			logWriter.println("Deduplicated Count: " + DECIMAL_FORMAT.format(duplicatesRemoved));
			logWriter.println("% Deduplicated: " + DECIMAL_FORMAT.format((double)duplicatesRemoved / (double)totalDNAProcessed.sum()));
			logWriter.println("Total Number of BP: " + DECIMAL_FORMAT.format(baseCount.sum()));
		}else if(mode == Mode.PAIRMERGE){
			logWriter.println("Forwards Read File: " + inputFile.getAbsolutePath());
			logWriter.println("Reversed Read File: " + inputFile2.getAbsolutePath());
			logWriter.println("Output Directory: " + outputDir);
			logWriter.println("Is Input GZIP Format: " + inputGZIP);
			logWriter.println("Is Output GZIP Format: " + outputGZIP);
			if(editMaxM < 0.0)
				logWriter.println("Edit Percentage to Count as Match For Merging Paired-Ends: " + -editMaxM);
			else
				logWriter.println("Max Edits to Count as Match For Merging Paired-Ends: " + editMaxM);
			logWriter.println("Interval to Print Processed Reads: " + DECIMAL_FORMAT.format(printProcessedInterval));
			logWriter.println("Interval to Print Duplicate Reads: " + DECIMAL_FORMAT.format(printDuplicateInterval));
			logWriter.println("Use Wildcard Characters: " + wildcard);
			logWriter.println("Remove Reads That Are Not Merged: " + removeNoMergeReads);
			logWriter.println("Parallel: " + parallel);
			logWriter.println("Parallel Batch Size: " + splitBatchSize);
			logWriter.println("Useable Processors: " + Runtime.getRuntime().availableProcessors());
			if(probM < 0.0)
				logWriter.println("Probability Based Matching for Merging Paired End Reads: false");
			else
				logWriter.println("Prior Probability for Paired End Merging: " + probM);
			logWriter.println();
			logWriter.println("Merging Pairs...");
			logWriter.println();
			logWriter.flush();
			
			startTime = System.currentTimeMillis();
			pairMergeFiles(); //merge paired-end reads
			
			logWriter.println();
			logWriter.println("Pair Merging Completed.");
			logWriter.println();
			logWriter.println("Total Run Time: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
			logWriter.println();
			logWriter.println("Number of Lines Processed: " + DECIMAL_FORMAT.format(totalDNAProcessed.sum() * 4));
			logWriter.println("Number of Reads Processed: " + DECIMAL_FORMAT.format(totalDNAProcessed.sum()));
			logWriter.println("Number of Reads Merged: " + DECIMAL_FORMAT.format(totalReadsMerged.sum()));
			logWriter.println("% Merged: " + DECIMAL_FORMAT.format((double)totalReadsMerged.sum() / (double)totalDNAProcessed.sum()));
			logWriter.println("Number of Reads Removed: " + DECIMAL_FORMAT.format(undeterminedDNA.sum()));
			logWriter.println("% Reads Removed: " + DECIMAL_FORMAT.format((double)undeterminedDNA.sum() / (double)totalDNAProcessed.sum()));
			logWriter.println("Total Number of BP: " + DECIMAL_FORMAT.format(baseCount.sum()));
		}else if(mode == Mode.FILTER){
			logWriter.println("Input File: " + inputFile.getAbsolutePath());
			logWriter.println("Replace Original File: " + replaceOriginal);
			logWriter.println("Output Directory: " + outputDir);
			logWriter.println("Is Input GZIP Format: " + inputGZIP);
			logWriter.println("Is Output GZIP Format: " + outputGZIP);
			if(editMaxA < 0.0)
				logWriter.println("Edit Percentage to Count as Match For Adapters: " + -editMaxA);
			else
				logWriter.println("Max Edits to Count as Match For Adapters: " + editMaxA);
			logWriter.println("Allow Insertions and Deletions for Adapters: " + allowIndelsA);
			logWriter.println("Interval to Print Processed Reads: " + DECIMAL_FORMAT.format(printProcessedInterval));
			logWriter.println("Interval to Print Duplicate Reads: " + DECIMAL_FORMAT.format(printDuplicateInterval));
			logWriter.println("Minimum Read Length: " + minLength);
			logWriter.println("Maximum Read Length: " + maxLength);
			if(removeDNAWithNPercent >= 0 && removeDNAWithNPercent <= 1.0)
				logWriter.println("Remove Reads With % of N Greater Than: " + removeDNAWithNPercent);
			else if(removeDNAWithNPercent < 0)
				logWriter.println("Remove Reads With at Least 1 N");
			else
				logWriter.println("Remove Reads With N: false");
			logWriter.println("Quality Filter Threshold: " + qualityFilter);
			logWriter.println("Quality Filter Algorithm: " + (filterAlgorithm ? "Error Sum" : "Average"));
			String adaptersFString = adaptersF.toString();
			logWriter.println("Read Adapters: " + adaptersFString.substring(1, adaptersFString.length() - 1));
			logWriter.println("Maximum Right Offset For Adapters: " + maxOffsetA);
			logWriter.println("Minimum Adapter Overlap: " + minOverlapA);
			logWriter.println("Quality Trim Algorithm: " + (trimAlgorithm ? "Sum" : "Local Average"));
			logWriter.println("5' Quality Trim Score Threshold: " + qualityTrimQScore1);
			logWriter.println("3' Quality Trim Score Threshold: " + qualityTrimQScore2);
			logWriter.println("Quality Trim Length: " + qualityTrimLength);
			if(trimNPercent > 1.0)
				logWriter.println("Trim Leading and Trailing N: false");
			else
				logWriter.println("Percentage to Trim Leading and Trailing N By: " + trimNPercent);
			logWriter.println("Use Wildcard Characters: " + wildcard);
			logWriter.println("Remove Reads That Are Not Quality Trimmed: " + removeUntrimmedReads);
			logWriter.println("Remove Reads That Do Not Contain Any Adapters: " + removeNoAdapterReads);
			logWriter.println("Parallel: " + parallel);
			logWriter.println("Parallel Batch Size: " + splitBatchSize);
			logWriter.println("Useable Processors: " + Runtime.getRuntime().availableProcessors());
			if(probA < 0.0)
				logWriter.println("Probability Based Matching for Adapters: false");
			else
				logWriter.println("Prior Probability for Adapter Matching: " + probA);
			logWriter.println();
			logWriter.println("Filtering...");
			logWriter.println();
			logWriter.flush();
			
			if(probA < 0.0){
				adapterPatternF = new ArrayList<HashMap<Character, BitVector>>();
				for(int i = 0; i < adaptersF.size(); i++){
					adapterPatternF.add(UtilMethods.genPatternMasks(adaptersF.get(i).isStart ? adaptersF.get(i).str : UtilMethods.reverse(adaptersF.get(i).str), allowIndelsA, wildcard));
				}
				
				if(inputFile2 != null){
					adapterPatternR = new ArrayList<HashMap<Character, BitVector>>();
					for(int i = 0; i < adaptersR.size(); i++){
						adapterPatternR.add(UtilMethods.genPatternMasks(adaptersR.get(i).isStart ? adaptersR.get(i).str : UtilMethods.reverse(adaptersR.get(i).str), allowIndelsA, wildcard));
					}
				}
			}
			
			startTime = System.currentTimeMillis();
			filterFile(); //filter the file
			
			logWriter.println();
			logWriter.println("Filtering Completed.");
			logWriter.println();
			logWriter.println("Total Run Time: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
			logWriter.println();
			logWriter.println("Number of Lines Processed: " + DECIMAL_FORMAT.format(totalDNAProcessed.sum() * 4));
			logWriter.println("Number of Reads Processed: " + DECIMAL_FORMAT.format(totalDNAProcessed.sum()));
			logWriter.println("Number of Reads With Removed Adapter: " + DECIMAL_FORMAT.format(totalRemovedAdapters.sum()));
			logWriter.println("% Removed Adapter: " + DECIMAL_FORMAT.format((double)totalRemovedAdapters.sum() / (double)totalDNAProcessed.sum()));
			logWriter.println("Number of Reads Quality Trimmed: " + DECIMAL_FORMAT.format(totalQualityTrimmed.sum()));
			logWriter.println("% Quality Trimmed: " + DECIMAL_FORMAT.format((double)totalQualityTrimmed.sum() / (double)totalDNAProcessed.sum()));
			logWriter.println("Number of Reads Removed: " + DECIMAL_FORMAT.format(undeterminedDNA.sum()));
			logWriter.println("% Reads Removed: " + DECIMAL_FORMAT.format((double)undeterminedDNA.sum() / (double)totalDNAProcessed.sum()));
			logWriter.println("Total Number of BP: " + DECIMAL_FORMAT.format(baseCount.sum()));
		}else if(mode == Mode.SIM_READS){
			if(simMerging){
				simMergingReads();
			}else{
				readSample();
				if(simUMI)
					simUMIReads();
				else
					simReads();
			}
		}
	}
	
	//read the sample file for demultiplex
	private void readSample() throws Exception{
		EnzymeList.init();
		BufferedReader reader = new BufferedReader(new FileReader(sampleFile), BUFFER_SIZE_LOG);
		
		String tempLine;
		boolean flag = false;
		while((tempLine = reader.readLine()) != null){
			String[] line = tempLine.split("\\s+");
			sampleMapF.put(line[1].toUpperCase(), line[0]);
			sampleDNAF.add(line[1].toUpperCase());
			if(line.length > 4){ //if there are reversed barcodes
				sampleDNAR.add(line[4].toUpperCase());
				hasReversedBarcode = true;
			}
			if(!flag){
				enzymeF = line[2];
				constEnzymesF = EnzymeList.enzymes.get(enzymeF.toUpperCase());
				if(line.length > 3){ //if there are reversed enzymes
					enzymeR = line[3];
					constEnzymesR = EnzymeList.enzymes.get(enzymeR.toUpperCase());
				}else{ //reversed enzymes will be forwards enzymes
					enzymeR = enzymeF;
					constEnzymesR = constEnzymesF;
				}
				flag = true;
			}
		}
		
		reader.close();
	}
	
	private void simMergingReads() throws Exception{
		BufferedWriter inWriterF = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "simulated_reads_R1.fastq.gz"))));
		BufferedWriter inWriterR = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "simulated_reads_R2.fastq.gz"))));
		BufferedWriter outWriterF = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_merged_R1.fastq.gz"))));
		
		long simMerged = 0;
		
		for(long i = 0; i < simIter; i++){
			Random r = new Random();
			
			String mergeF;
			String mergeR;
			
			if(r.nextFloat() >= 0.3f){ //have merge location
				if(r.nextFloat() >= 0.5f){ //have merge location within edit range
					simMerged += 2;
					mergeF = UtilMethods.randSeq(r, (int)editMaxM + r.nextInt(simReadLength - (int)editMaxM) + 1);
					mergeR = UtilMethods.randEdit(r, UtilMethods.reverseComplement(mergeF), r.nextInt((int)editMaxM + 1), false);
				}else{ //have merge location outside edit range
					mergeF = UtilMethods.randSeq(r, (int)editMaxM + r.nextInt(simReadLength - (int)editMaxM) + 1);
					mergeR = UtilMethods.randEdit(r, UtilMethods.reverseComplement(mergeF), (int)editMaxM + r.nextInt(Math.min((int)editMaxM * 2 + 1, mergeF.length()) - (int)editMaxM), false);
				}
			}else{ //no merge location
				mergeF = "";
				mergeR = "";
			}
			
			String seqF = UtilMethods.randSeq(r, simReadLength - mergeF.length());
			String readF = seqF + mergeF;
			String qualF = UtilMethods.randQuality(r, simReadLength);
			
			String seqR = UtilMethods.randSeq(r, simReadLength - mergeR.length());
			String readR = seqR + mergeR;
			String qualR = UtilMethods.randQuality(r, simReadLength);
			
			String mergedReadF = seqF + mergeF + UtilMethods.reverseComplement(seqR);
			String mergedQualF = UtilMethods.makeStr('A', mergedReadF.length()); //the real quality values don't really matter
			
			inWriterF.write("SIMULATED READ");
			inWriterF.newLine();
			inWriterF.write(readF);
			inWriterF.newLine();
			inWriterF.write(description2);
			inWriterF.newLine();
			inWriterF.write(qualF);
			inWriterF.newLine();
			
			inWriterR.write("SIMULATED READ");
			inWriterR.newLine();
			inWriterR.write(readR);
			inWriterR.newLine();
			inWriterR.write(description2);
			inWriterR.newLine();
			inWriterR.write(qualR);
			inWriterR.newLine();
			
			outWriterF.write("SIMULATED READ");
			outWriterF.newLine();
			outWriterF.write(mergedReadF);
			outWriterF.newLine();
			outWriterF.write(description2);
			outWriterF.newLine();
			outWriterF.write(mergedQualF);
			outWriterF.newLine();
		}
		
		inWriterF.close();
		inWriterR.close();
		outWriterF.close();
		
		logWriter.println("Number of Reads: " + DECIMAL_FORMAT.format(simIter * 2));
		logWriter.println("Number of Reads Merged: " + DECIMAL_FORMAT.format(simMerged));
	}
	
	private void simReads() throws Exception{
		BufferedWriter inWriterF = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "simulated_reads_R1.fastq.gz"))));
		BufferedWriter inWriterR = null;
		BufferedWriter[] outWriterF = new BufferedWriter[sampleMapF.size()];
		BufferedWriter[] outWriterR = null;
		BufferedWriter undeterminedWriterF = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_sample_undetermined_R1.fastq.gz"))));
		BufferedWriter undeterminedWriterR = null;
		
		for(int i = 0; i < sampleMapF.size(); i++){
			outWriterF[i] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_sample_" + sampleMapF.get(sampleDNAF.get(i)) + "_R1.fastq.gz"))));
		}
		
		if(simReversed){
			inWriterR = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "simulated_reads_R2.fastq.gz"))));
			outWriterR = new BufferedWriter[sampleMapF.size()];
			undeterminedWriterR = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_sample_undetermined_R2.fastq.gz"))));
			
			for(int i = 0; i < sampleMapF.size(); i++){
				outWriterR[i] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_sample_" + sampleMapF.get(sampleDNAF.get(i)) + "_R2.fastq.gz"))));
			}
		}
		
		long[][] simCounts = new long[sampleMapF.size()][2];
		long simUndetermined = 0;
		long simTotalAdapter = 0;
		long simUndeterminedAdapter = 0;
		
		for(long i = 0; i < simIter; i++){
			Random r = new Random();
			
			int barcode = r.nextInt(sampleMapF.size());
			boolean correctF = r.nextFloat() >= 0.1f;
			boolean correctR = r.nextFloat() >= 0.1f;
			boolean adapterF = !adaptersF.isEmpty() && r.nextFloat() >= 0.3f;
			boolean adapterR = simReversed && !adaptersR.isEmpty() && r.nextFloat() >= 0.3f;
			boolean undetermined = false;
			
			String seqF = null;
			String readF = null;
			String seqR = null;
			String readR = null;
			
			if(correctF && (!simReversed || correctR)){ //both forwards and reversed are correct, or forwards is correct for single end
				simCounts[barcode][0] += simReversed ? 2 : 1;
				
				seqF = UtilMethods.randSeq(r, simReadLength);
				readF = UtilMethods.randEdit(r, sampleDNAF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB) +
						UtilMethods.randEdit(r, constEnzymesF.get(r.nextInt(constEnzymesF.size())), r.nextInt((int)editMaxB + 1), allowIndelsB) + seqF;
				
				if(simReversed){
					seqR = UtilMethods.randSeq(r, simReadLength);
					readR = UtilMethods.randEdit(r, hasReversedBarcode ? sampleDNAR.get(barcode) : sampleDNAF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB) +
							UtilMethods.randEdit(r, constEnzymesR.get(r.nextInt(constEnzymesR.size())), r.nextInt((int)editMaxB + 1), allowIndelsB) + seqR;
				}
			}else if(correctF){ //forwards is correct, reversed is not
				simUndetermined += 2;
				undetermined = true;
				
				seqF = UtilMethods.randSeq(r, simReadLength);
				readF = UtilMethods.randEdit(r, sampleDNAF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB) +
						UtilMethods.randEdit(r, constEnzymesF.get(r.nextInt(constEnzymesF.size())), r.nextInt((int)editMaxB + 1), allowIndelsB) + seqF;
				
				seqR = UtilMethods.randSeq(r, simReadLength);
				String randBarcode;
				String randEnzyme;
				if(r.nextFloat() < 0.5f){ //error in barcode
					randBarcode = UtilMethods.randEdit(r, hasReversedBarcode ? sampleDNAR.get(barcode) : sampleDNAF.get(barcode),
							(int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, (hasReversedBarcode ? sampleDNAR.get(barcode) : sampleDNAF.get(barcode)).length()) - (int)editMaxB), allowIndelsB);
					randEnzyme = UtilMethods.randEdit(r, constEnzymesR.get(r.nextInt(constEnzymesR.size())), r.nextInt((int)editMaxB + 1), allowIndelsB);
				}else{ //error in enzyme
					randBarcode = UtilMethods.randEdit(r, hasReversedBarcode ? sampleDNAR.get(barcode) : sampleDNAF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB);
					int enzyme = r.nextInt(constEnzymesR.size());
					randEnzyme = UtilMethods.randEdit(r, constEnzymesR.get(enzyme), (int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, constEnzymesR.get(enzyme).length()) - (int)editMaxB), allowIndelsB);
				}
				readR = randBarcode + randEnzyme + seqR;
			}else if(simReversed && correctR){ //reversed is correct, forwards is not
				simUndetermined += 2;
				undetermined = true;
				
				seqF = UtilMethods.randSeq(r, simReadLength);
				String randBarcode;
				String randEnzyme;
				if(r.nextFloat() < 0.5f){ //error in barcode
					randBarcode = UtilMethods.randEdit(r, sampleDNAF.get(barcode), (int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, sampleDNAF.get(barcode).length()) - (int)editMaxB), allowIndelsB);
					randEnzyme = UtilMethods.randEdit(r, constEnzymesF.get(r.nextInt(constEnzymesF.size())), r.nextInt((int)editMaxB + 1), allowIndelsB);
				}else{ //error in enzyme
					randBarcode = UtilMethods.randEdit(r, sampleDNAF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB);
					int enzyme = r.nextInt(constEnzymesF.size());
					randEnzyme = UtilMethods.randEdit(r, constEnzymesF.get(enzyme), (int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, constEnzymesF.get(enzyme).length()) - (int)editMaxB), allowIndelsB);
				}
				readF = randBarcode + randEnzyme + seqF;
				
				seqR = UtilMethods.randSeq(r, simReadLength);
				readR = UtilMethods.randEdit(r, hasReversedBarcode ? sampleDNAR.get(barcode) : sampleDNAF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB) +
						UtilMethods.randEdit(r, constEnzymesR.get(r.nextInt(constEnzymesR.size())), r.nextInt((int)editMaxB + 1), allowIndelsB) + seqR;
			}else{ //forwards is not correct, reversed could be correct
				simUndetermined += simReversed ? 2 : 1;
				undetermined = true;
				
				seqF = UtilMethods.randSeq(r, simReadLength);
				String randBarcode;
				String randEnzyme;
				if(r.nextFloat() < 0.5f){ //error in barcode
					randBarcode = UtilMethods.randEdit(r, sampleDNAF.get(barcode), (int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, sampleDNAF.get(barcode).length()) - (int)editMaxB), allowIndelsB);
					randEnzyme = UtilMethods.randEdit(r, constEnzymesF.get(r.nextInt(constEnzymesF.size())), r.nextInt((int)editMaxB + 1), allowIndelsB);
				}else{ //error in enzyme
					randBarcode = UtilMethods.randEdit(r, sampleDNAF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB);
					int enzyme = r.nextInt(constEnzymesF.size());
					randEnzyme = UtilMethods.randEdit(r, constEnzymesF.get(enzyme), (int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, constEnzymesF.get(enzyme).length()) - (int)editMaxB), allowIndelsB);
				}
				readF = randBarcode + randEnzyme + seqF;
				
				if(simReversed){
					seqR = UtilMethods.randSeq(r, simReadLength);
					if(r.nextFloat() < 0.5f){ //error in barcode
						randBarcode = UtilMethods.randEdit(r, hasReversedBarcode ? sampleDNAR.get(barcode) : sampleDNAF.get(barcode),
								(int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, (hasReversedBarcode ? sampleDNAR.get(barcode) : sampleDNAF.get(barcode)).length()) - (int)editMaxB), allowIndelsB);
						randEnzyme = UtilMethods.randEdit(r, constEnzymesR.get(r.nextInt(constEnzymesR.size())), r.nextInt((int)editMaxB + 1), allowIndelsB);
					}else{ //error in enzyme
						randBarcode = UtilMethods.randEdit(r, hasReversedBarcode ? sampleDNAR.get(barcode) : sampleDNAF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB);
						int enzyme = r.nextInt(constEnzymesR.size());
						randEnzyme = UtilMethods.randEdit(r, constEnzymesR.get(enzyme), (int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, constEnzymesR.get(enzyme).length()) - (int)editMaxB), allowIndelsB);
					}
					readR = randBarcode + randEnzyme + seqR;
				}
			}
			
			if(adapterF){
				if(r.nextFloat() < 0.5f){ //correct
					simTotalAdapter++;
					if(undetermined)
						simUndeterminedAdapter++;
					else
						simCounts[barcode][1]++;
					
					int adapter = r.nextInt(adaptersF.size());
					readF += UtilMethods.randEdit(r, adaptersF.get(adapter).str.substring(0, Math.min(minOverlapA, adaptersF.get(adapter).str.length()) + r.nextInt(adaptersF.get(adapter).str.length() - Math.min(minOverlapA, adaptersF.get(adapter).str.length()) + 1)), r.nextInt((int)editMaxA + 1), allowIndelsA);
				}else{ //wrong
					int adapter = r.nextInt(adaptersF.size());
					if(r.nextFloat() < 0.5f){ //length less that minimum overlap
						readF += UtilMethods.randEdit(r, adaptersF.get(adapter).str.substring(0, Math.min(minOverlapA, adaptersF.get(adapter).str.length()) - 1), r.nextInt((int)editMaxA + 1), allowIndelsA);
					}else{ //error in adapter
						readF += UtilMethods.randEdit(r, adaptersF.get(adapter).str.substring(0, Math.min(minOverlapA, adaptersF.get(adapter).str.length()) + r.nextInt(adaptersF.get(adapter).str.length() - Math.min(minOverlapA, adaptersF.get(adapter).str.length()) + 1)),
								r.nextInt((int)editMaxA + r.nextInt(Math.min((int)editMaxA * 2 + 1, adaptersF.get(adapter).str.length()) - (int)editMaxA)), allowIndelsA);
					}
				}
			}
			
			if(adapterR){
				if(r.nextFloat() < 0.5f){ //correct
					simTotalAdapter++;
					if(undetermined)
						simUndeterminedAdapter++;
					else
						simCounts[barcode][1]++;
					
					int adapter = r.nextInt(adaptersR.size());
					readR += UtilMethods.randEdit(r, adaptersR.get(adapter).str.substring(0, Math.min(minOverlapA, adaptersR.get(adapter).str.length()) + r.nextInt(adaptersR.get(adapter).str.length() - Math.min(minOverlapA, adaptersR.get(adapter).str.length()) + 1)), r.nextInt((int)editMaxA + 1), allowIndelsA);
				}else{ //wrong
					int adapter = r.nextInt(adaptersR.size());
					if(r.nextFloat() < 0.5f){ //length less that minimum overlap
						readR += UtilMethods.randEdit(r, adaptersR.get(adapter).str.substring(0, Math.min(minOverlapA, adaptersR.get(adapter).str.length()) - 1), r.nextInt((int)editMaxA + 1), allowIndelsA);
					}else{ //error in adapter
						readR += UtilMethods.randEdit(r, adaptersR.get(adapter).str.substring(0, Math.min(minOverlapA, adaptersR.get(adapter).str.length()) + r.nextInt(adaptersR.get(adapter).str.length() - Math.min(minOverlapA, adaptersR.get(adapter).str.length()) + 1)),
								r.nextInt((int)editMaxA + r.nextInt(Math.min((int)editMaxA * 2 + 1, adaptersR.get(adapter).str.length()) - (int)editMaxA)), allowIndelsA);
					}
				}
			}
			
			inWriterF.write("SIMULATED READ");
			inWriterF.newLine();
			inWriterF.write(readF);
			inWriterF.newLine();
			inWriterF.write(description2);
			inWriterF.newLine();
			inWriterF.write(UtilMethods.makeStr('A', readF.length()));
			inWriterF.newLine();
			
			if(undetermined){
				undeterminedWriterF.write("SIMULATED READ");
				undeterminedWriterF.newLine();
				undeterminedWriterF.write(readF);
				undeterminedWriterF.newLine();
				undeterminedWriterF.write(description2);
				undeterminedWriterF.newLine();
				undeterminedWriterF.write(UtilMethods.makeStr('A', readF.length()));
				undeterminedWriterF.newLine();
			}else{
				outWriterF[barcode].write("SIMULATED READ");
				outWriterF[barcode].newLine();
				outWriterF[barcode].write(seqF);
				outWriterF[barcode].newLine();
				outWriterF[barcode].write(description2);
				outWriterF[barcode].newLine();
				outWriterF[barcode].write(UtilMethods.makeStr('A', seqF.length()));
				outWriterF[barcode].newLine();
			}
			
			if(simReversed){
				inWriterR.write("SIMULATED READ");
				inWriterR.newLine();
				inWriterR.write(readR);
				inWriterR.newLine();
				inWriterR.write(description2);
				inWriterR.newLine();
				inWriterR.write(UtilMethods.makeStr('A', readR.length()));
				inWriterR.newLine();
				
				if(undetermined){
					undeterminedWriterR.write("SIMULATED READ");
					undeterminedWriterR.newLine();
					undeterminedWriterR.write(readR);
					undeterminedWriterR.newLine();
					undeterminedWriterR.write(description2);
					undeterminedWriterR.newLine();
					undeterminedWriterR.write(UtilMethods.makeStr('A', readR.length()));
					undeterminedWriterR.newLine();
				}else{
					outWriterR[barcode].write("SIMULATED READ");
					outWriterR[barcode].newLine();
					outWriterR[barcode].write(seqR);
					outWriterR[barcode].newLine();
					outWriterR[barcode].write(description2);
					outWriterR[barcode].newLine();
					outWriterR[barcode].write(UtilMethods.makeStr('A', seqR.length()));
					outWriterR[barcode].newLine();
				}
			}
		}
		
		inWriterF.close();
		undeterminedWriterF.close();
		for(int i = 0; i < sampleMapF.size(); i++){
			outWriterF[i].close();
		}
		if(simReversed){
			inWriterR.close();
			undeterminedWriterR.close();
			for(int i = 0; i < sampleMapF.size(); i++){
				outWriterR[i].close();
			}
		}
		
		logWriter.println("Sample\tReads\tAdapters");
		for(int i = 0; i < sampleMapF.size(); i++){
			logWriter.println(sampleMapF.get(sampleDNAF.get(i)) + "\t" + DECIMAL_FORMAT.format(simCounts[i][0]) + "\t" + DECIMAL_FORMAT.format(simCounts[i][1]));
		}
		logWriter.println("Undetermined\t" + DECIMAL_FORMAT.format(simUndetermined) + "\t" + DECIMAL_FORMAT.format(simUndeterminedAdapter));
		logWriter.println("Total\t" + DECIMAL_FORMAT.format(simReversed ? simIter * 2 : simIter) + "\t" + DECIMAL_FORMAT.format(simTotalAdapter));
	}
	
	private void simUMIReads() throws Exception{
		BufferedWriter inWriterF = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "simulated_reads_R1.fastq.gz"))));
		BufferedWriter inWriterF2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "simulated_indexes_R1.fastq.gz"))));
		BufferedWriter inWriterR = null;
		BufferedWriter inWriterR2 = null;
		BufferedWriter[] outWriterF = new BufferedWriter[sampleMapF.size()];
		BufferedWriter[] outWriterF2 = new BufferedWriter[sampleMapF.size()];
		BufferedWriter[] outWriterR = null;
		BufferedWriter[] outWriterR2 = null;
		BufferedWriter undeterminedWriterF = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_sample_undetermined_R1.fastq.gz"))));
		BufferedWriter undeterminedWriterF2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_index_undetermined_R1.fastq.gz"))));
		BufferedWriter undeterminedWriterR = null;
		BufferedWriter undeterminedWriterR2 = null;
		
		for(int i = 0; i < sampleMapF.size(); i++){
			outWriterF[i] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_sample_" + sampleMapF.get(sampleDNAF.get(i)) + "_R1.fastq.gz"))));
			outWriterF2[i] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_index_" + sampleMapF.get(sampleDNAF.get(i)) + "_R1.fastq.gz"))));
		}
		
		if(simReversed){
			inWriterR = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "simulated_reads_R2.fastq.gz"))));
			inWriterR2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "simulated_indexes_R2.fastq.gz"))));
			outWriterR = new BufferedWriter[sampleMapF.size()];
			outWriterR2 = new BufferedWriter[sampleMapF.size()];
			undeterminedWriterR = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_sample_undetermined_R2.fastq.gz"))));
			undeterminedWriterR2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_index_undetermined_R2.fastq.gz"))));
			
			for(int i = 0; i < sampleMapF.size(); i++){
				outWriterR[i] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_sample_" + sampleMapF.get(sampleDNAF.get(i)) + "_R2.fastq.gz"))));
				outWriterR2[i] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_index_" + sampleMapF.get(sampleDNAF.get(i)) + "_R2.fastq.gz"))));
			}
		}
		
		long[][] simCounts = new long[sampleMapF.size()][1];
		long simUndetermined = 0;
		
		for(long i = 0; i < simIter; i++){
			Random r = new Random();
			
			int barcode = r.nextInt(sampleMapF.size());
			boolean correctF = r.nextFloat() >= 0.1f;
			boolean correctR = r.nextFloat() >= 0.1f;
			boolean undetermined = false;
			
			String readF = null;
			String indexF = null;
			String readR = null;
			String indexR = null;
			
			if(correctF && (!simReversed || correctR)){
				simCounts[barcode][0] += simReversed ? 2 : 1;
				
				readF = UtilMethods.randSeq(r, simReadLength);
				indexF = UtilMethods.randEdit(r, sampleDNAF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB) + UtilMethods.randSeq(r, randUMILength);
				
				if(simReversed){
					readR = UtilMethods.randSeq(r, simReadLength);
					indexR = UtilMethods.randEdit(r, hasReversedBarcode ? sampleDNAR.get(barcode) : sampleDNAF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB) + UtilMethods.randSeq(r, randUMILength);
				}
			}else if(correctF){
				simUndetermined += 2;
				undetermined = true;
				
				readF = UtilMethods.randSeq(r, simReadLength);
				indexF = UtilMethods.randEdit(r, sampleDNAF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB) + UtilMethods.randSeq(r, randUMILength);
				
				readR = UtilMethods.randSeq(r, simReadLength);
				indexR = UtilMethods.randEdit(r, hasReversedBarcode ? sampleDNAR.get(barcode) : sampleDNAF.get(barcode),
						(int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, (hasReversedBarcode ? sampleDNAR.get(barcode) : sampleDNAF.get(barcode)).length()) - (int)editMaxB), allowIndelsB) + UtilMethods.randSeq(r, randUMILength);
			}else if(simReversed && correctR){
				simUndetermined += 2;
				undetermined = true;
				
				readF = UtilMethods.randSeq(r, simReadLength);
				indexF = UtilMethods.randEdit(r, sampleDNAF.get(barcode), (int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, sampleDNAF.get(barcode).length()) - (int)editMaxB), allowIndelsB) + UtilMethods.randSeq(r, randUMILength);
				
				readR = UtilMethods.randSeq(r, simReadLength);
				indexR = UtilMethods.randEdit(r, hasReversedBarcode ? sampleDNAR.get(barcode) : sampleDNAF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB) + UtilMethods.randSeq(r, randUMILength);
			}else{
				simUndetermined += simReversed ? 2 : 1;
				undetermined = true;
				
				readF = UtilMethods.randSeq(r, simReadLength);
				indexF = UtilMethods.randEdit(r, sampleDNAF.get(barcode), (int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, sampleDNAF.get(barcode).length()) - (int)editMaxB), allowIndelsB) + UtilMethods.randSeq(r, randUMILength);
				
				readR = UtilMethods.randSeq(r, simReadLength);
				indexR = UtilMethods.randEdit(r, hasReversedBarcode ? sampleDNAR.get(barcode) : sampleDNAF.get(barcode),
						(int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, (hasReversedBarcode ? sampleDNAR.get(barcode) : sampleDNAF.get(barcode)).length()) - (int)editMaxB), allowIndelsB) + UtilMethods.randSeq(r, randUMILength);
			}
			
			inWriterF.write("SIMULATED READ");
			inWriterF.newLine();
			inWriterF.write(readF);
			inWriterF.newLine();
			inWriterF.write(description2);
			inWriterF.newLine();
			inWriterF.write(UtilMethods.makeStr('A', readF.length()));
			inWriterF.newLine();
			
			inWriterF2.write("SIMULATED INDEX");
			inWriterF2.newLine();
			inWriterF2.write(indexF);
			inWriterF2.newLine();
			inWriterF2.write(description2);
			inWriterF2.newLine();
			inWriterF2.write(UtilMethods.makeStr('A', indexF.length()));
			inWriterF2.newLine();
			
			if(undetermined){
				undeterminedWriterF.write("SIMULATED READ");
				undeterminedWriterF.newLine();
				undeterminedWriterF.write(readF);
				undeterminedWriterF.newLine();
				undeterminedWriterF.write(description2);
				undeterminedWriterF.newLine();
				undeterminedWriterF.write(UtilMethods.makeStr('A', readF.length()));
				undeterminedWriterF.newLine();
				
				undeterminedWriterF2.write("SIMULATED INDEX");
				undeterminedWriterF2.newLine();
				undeterminedWriterF2.write(indexF);
				undeterminedWriterF2.newLine();
				undeterminedWriterF2.write(description2);
				undeterminedWriterF2.newLine();
				undeterminedWriterF2.write(UtilMethods.makeStr('A', indexF.length()));
				undeterminedWriterF2.newLine();
			}else{
				outWriterF[barcode].write("SIMULATED READ");
				outWriterF[barcode].newLine();
				outWriterF[barcode].write(readF);
				outWriterF[barcode].newLine();
				outWriterF[barcode].write(description2);
				outWriterF[barcode].newLine();
				outWriterF[barcode].write(UtilMethods.makeStr('A', readF.length()));
				outWriterF[barcode].newLine();
				
				outWriterF2[barcode].write("SIMULATED INDEX");
				outWriterF2[barcode].newLine();
				outWriterF2[barcode].write(indexF);
				outWriterF2[barcode].newLine();
				outWriterF2[barcode].write(description2);
				outWriterF2[barcode].newLine();
				outWriterF2[barcode].write(UtilMethods.makeStr('A', indexF.length()));
				outWriterF2[barcode].newLine();
			}
			
			if(simReversed){
				inWriterR.write("SIMULATED READ");
				inWriterR.newLine();
				inWriterR.write(readR);
				inWriterR.newLine();
				inWriterR.write(description2);
				inWriterR.newLine();
				inWriterR.write(UtilMethods.makeStr('A', readR.length()));
				inWriterR.newLine();
				
				inWriterR2.write("SIMULATED INDEX");
				inWriterR2.newLine();
				inWriterR2.write(indexR);
				inWriterR2.newLine();
				inWriterR2.write(description2);
				inWriterR2.newLine();
				inWriterR2.write(UtilMethods.makeStr('A', indexR.length()));
				inWriterR2.newLine();
				
				if(undetermined){
					undeterminedWriterR.write("SIMULATED READ");
					undeterminedWriterR.newLine();
					undeterminedWriterR.write(readR);
					undeterminedWriterR.newLine();
					undeterminedWriterR.write(description2);
					undeterminedWriterR.newLine();
					undeterminedWriterR.write(UtilMethods.makeStr('A', readR.length()));
					undeterminedWriterR.newLine();
					
					undeterminedWriterR2.write("SIMULATED INDEX");
					undeterminedWriterR2.newLine();
					undeterminedWriterR2.write(indexR);
					undeterminedWriterR2.newLine();
					undeterminedWriterR2.write(description2);
					undeterminedWriterR2.newLine();
					undeterminedWriterR2.write(UtilMethods.makeStr('A', indexR.length()));
					undeterminedWriterR2.newLine();
				}else{
					outWriterR[barcode].write("SIMULATED READ");
					outWriterR[barcode].newLine();
					outWriterR[barcode].write(readR);
					outWriterR[barcode].newLine();
					outWriterR[barcode].write(description2);
					outWriterR[barcode].newLine();
					outWriterR[barcode].write(UtilMethods.makeStr('A', readR.length()));
					outWriterR[barcode].newLine();
					
					outWriterR2[barcode].write("SIMULATED INDEX");
					outWriterR2[barcode].newLine();
					outWriterR2[barcode].write(indexR);
					outWriterR2[barcode].newLine();
					outWriterR2[barcode].write(description2);
					outWriterR2[barcode].newLine();
					outWriterR2[barcode].write(UtilMethods.makeStr('A', indexR.length()));
					outWriterR2[barcode].newLine();
				}
			}
		}
		
		inWriterF.close();
		inWriterF2.close();
		undeterminedWriterF.close();
		undeterminedWriterF2.close();
		for(int i = 0; i < sampleMapF.size(); i++){
			outWriterF[i].close();
			outWriterF2[i].close();
		}
		if(simReversed){
			inWriterR.close();
			inWriterR2.close();
			undeterminedWriterR.close();
			undeterminedWriterR2.close();
			for(int i = 0; i < sampleMapF.size(); i++){
				outWriterR[i].close();
				outWriterR2[i].close();
			}
		}
		
		logWriter.println("Sample\tReads");
		for(int i = 0; i < sampleMapF.size(); i++){
			logWriter.println(sampleMapF.get(sampleDNAF.get(i)) + "\t" + DECIMAL_FORMAT.format(simCounts[i][0]));
		}
		logWriter.println("Undetermined\t" + DECIMAL_FORMAT.format(simUndetermined));
		logWriter.println("Total\t" + DECIMAL_FORMAT.format(simReversed ? simIter * 2 : simIter));
	}
	
	private ConcurrentLinkedQueue<Strings> demultiplexFile() throws Exception{
		BufferedReader reader1; //forwards input
		BufferedReader reader2 = null; //reversed input
		BufferedReader reader3 = null; //forwards index input
		BufferedReader reader4 = null; //reversed index input
		if(inputGZIP){
			reader1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFile), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
			if(indexFile != null)
				reader3 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(indexFile), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
			if(inputFile2 != null){
				reader2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFile2), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
				if(indexFile2 != null)
					reader4 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(indexFile2), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
			}
		}else{
			reader1 = new BufferedReader(new FileReader(inputFile), BUFFER_SIZE);
			if(indexFile != null)
				reader3 = new BufferedReader(new FileReader(indexFile), BUFFER_SIZE);
			if(inputFile2 != null){
				reader2 = new BufferedReader(new FileReader(inputFile2), BUFFER_SIZE);
				if(indexFile2 != null)
					reader4 = new BufferedReader(new FileReader(indexFile2), BUFFER_SIZE);
			}
		}
		BufferedWriter[] writers1 = new BufferedWriter[sampleMapF.size() + 1]; //forwards output
		BufferedWriter[] writers2 = !mergePairedEnds && inputFile2 != null ? new BufferedWriter[sampleMapF.size()] : null; //reversed output
		BufferedWriter[] writers3 = indexFile != null ? new BufferedWriter[sampleMapF.size() + 1] : null; //forwards index output
		BufferedWriter[] writers4 = !mergePairedEnds && inputFile2 != null && indexFile2 != null ? writers4 = new BufferedWriter[sampleMapF.size()] : null; //reversed index output
		BufferedWriter[] undeterminedWriterR = !mergePairedEnds && inputFile2 != null ? undeterminedWriterR = new BufferedWriter[2] : null; //for reversed and reversed index undetermined reads
		
		ConcurrentLinkedQueue<Strings> files = new ConcurrentLinkedQueue<Strings>();
		LongAdder[][] DNACounts = new LongAdder[sampleMapF.size() + 1][5];
		for(int i = 0; i < DNACounts.length; i++){
			for(int j = 0; j < DNACounts[i].length; j++){
				DNACounts[i][j] = new LongAdder();
			}
		}
		Object[] locks = new Object[sampleMapF.size() + 1];
		for(int i = 0; i < locks.length; i++){
			locks[i] = new Object();
		}
		
		StreamSupport.stream(new ReadSpliterator<Read>(splitBatchSize, reader1, reader2, reader3, reader4), parallel).forEach((read) -> {
			totalDNAProcessed.add(inputFile2 == null ? 1 : 2);
			baseCount.add(read.readF[1].length() + (inputFile2 == null ? 0 : read.readR[1].length()));
			boolean flag = false;
			if(!parallel && totalDNAProcessed.sum() % printProcessedInterval == 0){
				logWriter.println("Reads Processed So Far: " + DECIMAL_FORMAT.format(totalDNAProcessed));
				logWriter.println("Run Time So Far: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
				logWriter.println("Current Memory Usage (GB): " + DECIMAL_FORMAT.format(UtilMethods.currentMemoryUsage()));
				logWriter.flush();
			}
			
			//save the quality and sequence lines so the undetermined reads are not trimmed or anything
			String tempF1 = read.readF[1];
			String tempF2 = read.readF[3];
			String tempR1 = null;
			String tempR2 = null;
			if(inputFile2 != null){
				tempR1 = read.readR[1];
				tempR2 = read.readR[3];
			}
			
			//check the percentage of N
			if(removeDNAWithNPercent > 1.0 || (removeDNAWithNPercent >= 0.0 && UtilMethods.percentN(read.readF[1]) <= removeDNAWithNPercent && (inputFile2 == null || UtilMethods.percentN(read.readR[1]) <= removeDNAWithNPercent)) ||
					(removeDNAWithNPercent < 0.0 && UtilMethods.countN(read.readF[1]) <= 0 && (inputFile2 == null || UtilMethods.countN(read.readR[1]) <= 0))){
				int barcodeIndex = -1;
				int barcodeEnd = -1;
				int barcodeLength = 0;
				int enzymeEnd = -1;
				int barcodeEnd2 = -1;
				int enzymeEnd2 = -1;
				int minEdit = Integer.MAX_VALUE;
				int barcodeMatchCount = 0;
				
				if(indexFile == null){ //check for enzyme and barcode in forwards reads
					for(int i = 0; i < sampleDNAF.size(); i++){
						ArrayList<Match> matches = null;
						if(probB < 0.0){
							matches = UtilMethods.searchWithN(read.readF[1], 0, Math.min(read.readF[1].length(), maxOffsetB + sampleDNAF.get(i).length() +
									(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * sampleDNAF.get(i).length()) : (int)editMaxB) : 0)), sampleDNAF.get(i), editMaxB, allowIndelsB, false, minOverlapB, wildcard, barcodePatternF.get(i));
						}else{
							matches = UtilMethods.searchWithProb(read.readF[1], 0, Math.min(read.readF[1].length(), maxOffsetB + sampleDNAF.get(i).length() +
									(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * sampleDNAF.get(i).length()) : (int)editMaxB) : 0)), read.readF[3], sampleDNAF.get(i), null, probB, minOverlapB, wildcard);
						}
						boolean isMatch = false;
						for(int j = 0; j < matches.size(); j++){
							for(int k = 0; k < constEnzymesF.size(); k++){
								ArrayList<Match> enzymeMatches = null;
								if(probB < 0.0){
									enzymeMatches = UtilMethods.searchWithN(read.readF[1], matches.get(j).end + 1/* + randUMILength*/, Math.min(read.readF[1].length(), maxOffsetB + matches.get(j).end + 1 + /*randUMILength + */constEnzymesF.get(k).length() +
											(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * constEnzymesF.get(k).length()) : (int)editMaxB) : 0)), constEnzymesF.get(k), editMaxB, allowIndelsB, true, Integer.MAX_VALUE, wildcard, enzymePatternF.get(k));
								}else{
									enzymeMatches = UtilMethods.searchWithProb(read.readF[1], matches.get(j).end + 1/* + randUMILength*/, Math.min(read.readF[1].length(), maxOffsetB + matches.get(j).end + 1 + /*randUMILength + */constEnzymesF.get(k).length() +
											(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * constEnzymesF.get(k).length()) : (int)editMaxB) : 0)), read.readF[3], constEnzymesF.get(k), null, probB, minOverlapB, wildcard);
								}
								if(!enzymeMatches.isEmpty()){
									int tempBarcodeEnd2 = -1;
									int tempBarcodeLength2 = 0;
									int tempEnzymeEnd2 = -1;
									if(inputFile2 != null && checkReversedReads){
										int minEdit2 = Integer.MAX_VALUE;
										ArrayList<Match> rMatches = null;
										if(probB < 0.0){
											rMatches = UtilMethods.searchWithN(read.readR[1], 0, Math.min(read.readR[1].length(), maxOffsetB + (hasReversedBarcode ? sampleDNAR.get(i).length() : sampleDNAF.get(i).length()) +
													(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * (hasReversedBarcode ? sampleDNAR.get(i).length() : sampleDNAF.get(i).length())) : (int)editMaxB) : 0)),
													hasReversedBarcode ? sampleDNAR.get(i) : UtilMethods.complement(sampleDNAF.get(i)), editMaxB, allowIndelsB, false, minOverlapB, wildcard, barcodePatternR.get(i));
										}else{
											rMatches = UtilMethods.searchWithProb(read.readR[1], 0, Math.min(read.readR[1].length(), maxOffsetB + (hasReversedBarcode ? sampleDNAR.get(i).length() : sampleDNAF.get(i).length()) +
													(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * (hasReversedBarcode ? sampleDNAR.get(i).length() : sampleDNAF.get(i).length())) : (int)editMaxB) : 0)), read.readR[3],
													hasReversedBarcode ? sampleDNAR.get(i) : UtilMethods.complement(sampleDNAF.get(i)), null, probB, minOverlapB, wildcard);
										}
										for(int ri = 0; ri < rMatches.size(); ri++){
											for(int rj = 0; rj < constEnzymesR.size(); rj++){
												ArrayList<Match> rEnzymeMatches = null;
												if(probB < 0.0){
													rEnzymeMatches = UtilMethods.searchWithN(read.readR[1], rMatches.get(ri).end + 1/* + randUMILength*/, Math.min(read.readR[1].length(), maxOffsetB + rMatches.get(ri).end + 1 + /*randUMILength + */constEnzymesR.get(rj).length() +
															(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * constEnzymesR.get(rj).length()) : (int)editMaxB) : 0)), constEnzymesR.get(rj), editMaxB, allowIndelsB, true, Integer.MAX_VALUE, wildcard, enzymePatternR.get(rj));
												}else{
													rEnzymeMatches = UtilMethods.searchWithProb(read.readR[1], rMatches.get(ri).end + 1/* + randUMILength*/, Math.min(read.readR[1].length(), maxOffsetB + rMatches.get(ri).end + 1 + /*randUMILength + */constEnzymesR.get(rj).length() +
															(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * constEnzymesR.get(rj).length()) : (int)editMaxB) : 0)), read.readR[3], constEnzymesR.get(rj), null, probB, minOverlapB, wildcard);
												}
												if(!rEnzymeMatches.isEmpty()){
													if(rMatches.get(ri).edits <= minEdit2 && (rMatches.get(ri).edits < minEdit2 || rMatches.get(ri).length > tempBarcodeLength2)){
														tempBarcodeEnd2 = rMatches.get(ri).end + 1;
														tempBarcodeLength2 = rMatches.get(ri).length;
														tempEnzymeEnd2 = rEnzymeMatches.get(rEnzymeMatches.size() - 1).end + 1;
														minEdit2 = rMatches.get(ri).edits;
														break;
													}
												}
											}
										}
									}
									
									if(inputFile2 == null || !checkReversedReads || tempEnzymeEnd2 != -1){
										isMatch = true;
										if(matches.get(j).edits <= minEdit && (matches.get(j).edits < minEdit || matches.get(j).length > barcodeLength)){
											barcodeIndex = i;
											barcodeEnd = matches.get(j).end + 1;
											barcodeLength = matches.get(j).length;
											enzymeEnd = enzymeMatches.get(enzymeMatches.size() - 1).end + 1;
											if(inputFile2 != null){
												if(checkReversedReads){
													barcodeEnd2 = tempBarcodeEnd2;
													enzymeEnd2 = tempEnzymeEnd2;
												}else{
													barcodeEnd2 = Math.max(0, Math.min(read.readR[1].length(), barcodeEnd/* - randUMILength*/));
													enzymeEnd2 = Math.max(0, Math.min(read.readR[1].length(), enzymeEnd/* - randUMILength*/));
												}
											}
											minEdit = matches.get(j).edits;
											break;
										}
									}
								}
							}
						}
						if(isMatch){
							if(singleBarcodeMatchOnly && barcodeMatchCount >= 1){
								barcodeIndex = -1;
								break;
							}
							barcodeMatchCount++;
						}
					}
				}else{ //check for barcode in index reads
					for(int i = 0; i < sampleDNAF.size(); i++){
						ArrayList<Match> matches = null;
						if(probB < 0.0){
							matches = UtilMethods.searchWithN(read.readIF[1], 0, Math.min(read.readIF[1].length(), sampleDNAF.get(i).length() +
									(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * sampleDNAF.get(i).length()) : (int)editMaxB) : 0)),
									sampleDNAF.get(i), editMaxB, allowIndelsB, true, minOverlapB, wildcard, barcodePatternF.get(i));
						}else{
							matches = UtilMethods.searchWithProb(read.readIF[1], 0, Math.min(read.readIF[1].length(), sampleDNAF.get(i).length() +
									(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * sampleDNAF.get(i).length()) : (int)editMaxB) : 0)),
									read.readIF[3], sampleDNAF.get(i), null, probB, minOverlapB, wildcard);
						}
						if(!matches.isEmpty()){
							ArrayList<Match> rMatches = null;
							if(inputFile2 != null && checkReversedReads){
								if(probB < 0.0){
									rMatches = UtilMethods.searchWithN(read.readIR[1], randUMILength, Math.min(read.readIR[1].length(), randUMILength +
											(hasReversedBarcode ? sampleDNAR.get(i).length() : sampleDNAF.get(i).length()) + (probB < 0.0 && allowIndelsB ?
													(editMaxB < 0.0 ? (int)(-editMaxB * (hasReversedBarcode ? sampleDNAR.get(i).length() : sampleDNAF.get(i).length())) : (int)editMaxB) : 0)),
											hasReversedBarcode ? sampleDNAR.get(i) : UtilMethods.complement(sampleDNAF.get(i)), editMaxB, allowIndelsB, true, Integer.MAX_VALUE, wildcard, barcodePatternR.get(i));
								}else{
									rMatches = UtilMethods.searchWithProb(read.readIR[1], randUMILength, Math.min(read.readIR[1].length(), randUMILength +
											(hasReversedBarcode ? sampleDNAR.get(i).length() : sampleDNAF.get(i).length()) + (probB < 0.0 && allowIndelsB ?
													(editMaxB < 0.0 ? (int)(-editMaxB * (hasReversedBarcode ? sampleDNAR.get(i).length() : sampleDNAF.get(i).length())) : (int)editMaxB) : 0)),
											read.readIR[3], hasReversedBarcode ? sampleDNAR.get(i) : UtilMethods.complement(sampleDNAF.get(i)), null, probB, minOverlapB, wildcard);
								}
							}
							if((inputFile2 == null || !checkReversedReads || !rMatches.isEmpty()) && matches.get(matches.size() - 1).edits <= minEdit && (matches.get(matches.size() - 1).edits < minEdit || matches.get(matches.size() - 1).end + 1 > barcodeEnd)){
								barcodeIndex = i;
								barcodeEnd = matches.get(matches.size() - 1).end + 1;
								minEdit = matches.get(matches.size() - 1).edits;
							}
							if(singleBarcodeMatchOnly && barcodeMatchCount >= 1){
								barcodeIndex = -1;
								break;
							}
							barcodeMatchCount++;
						}
					}
				}
				
				if(barcodeIndex != -1){ //if barcode and enzyme is found for forwards and reversed reads
					//check if the quality is good enough
					boolean qualityAcceptable;
					if(filterAlgorithm){
						qualityAcceptable = UtilMethods.toError(read.readF[3], 0) <= qualityFilter && (inputFile2 == null || UtilMethods.toError(read.readR[3], 0) <= qualityFilter);
					}else{
						qualityAcceptable = UtilMethods.toQScore(read.readF[3], 0) >= qualityFilter && (inputFile2 == null || UtilMethods.toQScore(read.readR[3], 0) >= qualityFilter);
					}
					
					if(qualityAcceptable){
						//initialize file writers if they are not already initialized
						String newSequence1;
						String newQuality1 = null;
						String newSequence2 = null;
						String newQuality2 = null;
						
						//remove barcode and enzyme
						if(indexFile == null){
							if(removeEnzyme && removeBarRand){
								newSequence1 = read.readF[1].substring(enzymeEnd);
								newQuality1 = read.readF[3].substring(enzymeEnd);
							}else if(removeEnzyme){
								newSequence1 = read.readF[1].substring(0, barcodeEnd/* + randUMILength*/) + read.readF[1].substring(enzymeEnd);
								newQuality1 = read.readF[3].substring(0, barcodeEnd/* + randUMILength*/) + read.readF[3].substring(enzymeEnd);
							}else if(removeBarRand){
								newSequence1 = read.readF[1].substring(barcodeEnd/* + randUMILength*/);
								newQuality1 = read.readF[3].substring(barcodeEnd/* + randUMILength*/);
							}else{
								newSequence1 = read.readF[1];
								newQuality1 = read.readF[3];
							}
						}else{
							newSequence1 = read.readF[1];
							newQuality1 = read.readF[3];
						}
						
						if(inputFile2 != null){
							if(indexFile2 == null){
								if(removeEnzyme && removeBarRand){
									newSequence2 = read.readR[1].substring(enzymeEnd2);
									newQuality2 = read.readR[3].substring(enzymeEnd2);
								}else if(removeEnzyme){
									newSequence2 = read.readR[1].substring(0, barcodeEnd2) + read.readR[1].substring(enzymeEnd2);
									newQuality2 = read.readR[3].substring(0, barcodeEnd2) + read.readR[3].substring(enzymeEnd2);
								}else if(removeBarRand){
									newSequence2 = read.readR[1].substring(barcodeEnd2);
									newQuality2 = read.readR[3].substring(barcodeEnd2);
								}else{
									newSequence2 = read.readR[1];
									newQuality2 = read.readR[3];
								}
							}else{
								newSequence2 = read.readR[1];
								newQuality2 = read.readR[3];
							}
						}
						
						String[] temp;
						String[] qualityTrimmed;
						int trimmedQuality = 0;
						int removedAdapter = 0;
						
						//trim N
						temp = UtilMethods.trimN(newSequence1, newQuality1, trimNPercent);
						newSequence1 = temp[0];
						newQuality1 = temp[1];
						
						//quality trim
						if(trimAlgorithm){
							temp = UtilMethods.qualityTrim2(newSequence1, newQuality1, qualityTrimQScore1, true, qualityTrimLength);
							qualityTrimmed = UtilMethods.qualityTrim2(temp[0], temp[1], qualityTrimQScore2, false, qualityTrimLength);
						}else{
							temp = UtilMethods.qualityTrim1(newSequence1, newQuality1, qualityTrimQScore1, true, qualityTrimLength);
							qualityTrimmed = UtilMethods.qualityTrim1(temp[0], temp[1], qualityTrimQScore2, false, qualityTrimLength);
						}
						if(newSequence1.length() != qualityTrimmed[0].length()){
							trimmedQuality++;
						}
						//remove adapters
						String[] removedAdapters = UtilMethods.removeAdapters(qualityTrimmed[0], qualityTrimmed[1], adaptersF, editMaxA, minOverlapA, maxOffsetA, allowIndelsA, probA, wildcard, adapterPatternF);
						if(qualityTrimmed[0].length() != removedAdapters[0].length()){
							removedAdapter++;
						}
						newSequence1 = removedAdapters[0];
						newQuality1 = removedAdapters[1];
						
						//do the same for reversed reads
						if(inputFile2 != null){
							temp = UtilMethods.trimN(newSequence2, newQuality2, trimNPercent);
							newSequence2 = temp[0];
							newQuality2 = temp[1];
							
							if(trimAlgorithm){
								temp = UtilMethods.qualityTrim2(newSequence2, newQuality2, qualityTrimQScore1, true, qualityTrimLength);
								qualityTrimmed = UtilMethods.qualityTrim2(temp[0], temp[1], qualityTrimQScore2, false, qualityTrimLength);
							}else{
								temp = UtilMethods.qualityTrim1(newSequence2, newQuality2, qualityTrimQScore1, true, qualityTrimLength);
								qualityTrimmed = UtilMethods.qualityTrim1(temp[0], temp[1], qualityTrimQScore2, false, qualityTrimLength);
							}
							if(newSequence2.length() != qualityTrimmed[0].length()){
								trimmedQuality++;
							}
							removedAdapters = UtilMethods.removeAdapters(qualityTrimmed[0], qualityTrimmed[1], adaptersR, editMaxA, minOverlapA, maxOffsetA, allowIndelsA, probA, wildcard, adapterPatternR);
							if(qualityTrimmed[0].length() != removedAdapters[0].length()){
								removedAdapter++;
							}
							newSequence2 = removedAdapters[0];
							newQuality2 = removedAdapters[1];
						}
						
						boolean mergedReads = false;
						//merge paired-end reads
						if(mergePairedEnds && inputFile2 != null){
							String[] merged = UtilMethods.mergeReads(newSequence1, newQuality1, newSequence2, newQuality2, editMaxM, probM, wildcard);
							if(merged[0].length() != newSequence1.length() + newSequence2.length()){
								mergedReads = true;
							}
							newSequence1 = merged[0];
							newQuality1 = merged[1];
						}
						
						//check if reads are trimmed and if length is too long or short
						if((!removeUntrimmedReads || trimmedQuality == (inputFile2 == null ? 1 : 2)) && (!removeNoAdapterReads || removedAdapter == (inputFile2 == null ? 1 : 2)) && (!removeNoMergeReads || mergedReads) &&
								minLength <= newSequence1.length() && newSequence1.length() <= maxLength && (inputFile2 == null || mergePairedEnds || (minLength <= newSequence2.length() && newSequence2.length() <= maxLength))){
							//statistics
							DNACounts[barcodeIndex][0].add(inputFile2 == null ? 1 : 2);
							DNACounts[barcodeIndex][1].add(read.readF[1].length() + (inputFile2 == null ? 0 : read.readR[1].length()));
							DNACounts[barcodeIndex][3].add(removedAdapter);
							totalRemovedAdapters.add(removedAdapter);
							DNACounts[barcodeIndex][2].add(mergedReads ? 2 : 0);
							totalReadsMerged.add(mergedReads ? 2 : 0);
							DNACounts[barcodeIndex][4].add(trimmedQuality);
							totalQualityTrimmed.add(trimmedQuality);
							
							synchronized(locks[barcodeIndex]){
								try{
									Strings strings = new Strings();
									if(writers1[barcodeIndex] == null){
										if(!outputGZIP && !removeFirstDup && !removeBestDup){
											writers1[barcodeIndex] = new BufferedWriter(new FileWriter(outputDir + "sample_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R1.fastq"), BUFFER_SIZE);
											strings.add(outputDir + "sample_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R1.fastq");
											if(indexFile != null){
												writers3[barcodeIndex] = new BufferedWriter(new FileWriter(outputDir + "index_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R1.fastq"), BUFFER_SIZE);
												strings.add(outputDir + "index_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R1.fastq");
											}
										}else{
											writers1[barcodeIndex] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(((removeFirstDup || removeBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "sample_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R1.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
											strings.add(((removeFirstDup || removeBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "sample_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R1.fastq.gz");
											if(indexFile != null){
												writers3[barcodeIndex] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(((removeFirstDup || removeBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "index_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R1.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
												strings.add(((removeFirstDup || removeBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "index_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R1.fastq.gz");
											}
										}
									}
									if(!mergePairedEnds && inputFile2 != null && writers2[barcodeIndex] == null){
										if(!outputGZIP && !removeFirstDup && !removeBestDup){
											writers2[barcodeIndex] = new BufferedWriter(new FileWriter(outputDir + "sample_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R2.fastq"), BUFFER_SIZE);
											strings.add(outputDir + "sample_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R2.fastq");
											if(indexFile2 != null){
												writers4[barcodeIndex] = new BufferedWriter(new FileWriter(outputDir + "index_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R2.fastq"), BUFFER_SIZE);
												strings.add(outputDir + "index_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R2.fastq");
											}
										}else{
											writers2[barcodeIndex] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(((removeFirstDup || removeBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "sample_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R2.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
											strings.add(((removeFirstDup || removeBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "sample_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R2.fastq.gz");
											if(indexFile2 != null){
												writers4[barcodeIndex] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(((removeFirstDup || removeBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "index_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R2.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
												strings.add(((removeFirstDup || removeBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "index_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R2.fastq.gz");
											}
										}
									}
									if(strings.size() > 0)
										files.offer(strings);
									
									//print to the whatever file the read belongs to
									writers1[barcodeIndex].write(read.readF[0]);
									writers1[barcodeIndex].newLine();
									writers1[barcodeIndex].write(newSequence1);
									writers1[barcodeIndex].newLine();
									writers1[barcodeIndex].write(description2);
									writers1[barcodeIndex].newLine();
									writers1[barcodeIndex].write(newQuality1);
									writers1[barcodeIndex].newLine();
									
									if(indexFile != null){
										writers3[barcodeIndex].write(read.readIF[0]);
										writers3[barcodeIndex].newLine();
										writers3[barcodeIndex].write(read.readIF[1]);
										writers3[barcodeIndex].newLine();
										writers3[barcodeIndex].write(description2);
										writers3[barcodeIndex].newLine();
										writers3[barcodeIndex].write(read.readIF[3]);
										writers3[barcodeIndex].newLine();
									}
									
									if(!mergePairedEnds && inputFile2 != null){
										writers2[barcodeIndex].write(read.readR[0]);
										writers2[barcodeIndex].newLine();
										writers2[barcodeIndex].write(newSequence2);
										writers2[barcodeIndex].newLine();
										writers2[barcodeIndex].write(description2);
										writers2[barcodeIndex].newLine();
										writers2[barcodeIndex].write(newQuality2);
										writers2[barcodeIndex].newLine();
										
										if(indexFile2 != null){
											writers4[barcodeIndex].write(read.readIR[0]);
											writers4[barcodeIndex].newLine();
											writers4[barcodeIndex].write(read.readIR[1]);
											writers4[barcodeIndex].newLine();
											writers4[barcodeIndex].write(description2);
											writers4[barcodeIndex].newLine();
											writers4[barcodeIndex].write(read.readIR[3]);
											writers4[barcodeIndex].newLine();
										}
									}
								}catch(Exception e){
									UtilMethods.defaultExceptionHandler(logWriter, e);
								}
							}
							flag = true;
						}
					}
				}
			}
			if(!flag){
				//statistics
				DNACounts[DNACounts.length - 1][0].add(inputFile2 == null ? 1 : 2);
				DNACounts[DNACounts.length - 1][1].add(read.readF[1].length() + (inputFile2 == null ? 0 : read.readR[1].length()));
				undeterminedDNA.add(inputFile2 == null ? 1 : 2);
				
				synchronized(locks[locks.length - 1]){
					try{
						if(writers1[writers1.length - 1] == null){ //initialize writer for the undetermined file if it has not already been initialized
							if(outputGZIP){
								writers1[writers1.length - 1] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "sample_undetermined_R1.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
								if(inputFile2 != null)
									undeterminedWriterR[0] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "sample_undetermined_R2.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
								if(indexFile != null){
									writers3[writers3.length - 1] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "index_undetermined_R1.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
									if(inputFile2 != null)
										undeterminedWriterR[1] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "index_undetermined_R2.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
								}
							}else{
								writers1[writers1.length - 1] = new BufferedWriter(new FileWriter(outputDir + "sample_undetermined_R1.fastq"), BUFFER_SIZE);
								if(inputFile2 != null)
									undeterminedWriterR[0] = new BufferedWriter(new FileWriter(outputDir + "sample_undetermined_R2.fastq"), BUFFER_SIZE);
								if(indexFile != null){
									writers3[writers3.length - 1] = new BufferedWriter(new FileWriter(outputDir + "index_undetermined_R1.fastq"), BUFFER_SIZE);
									if(inputFile2 != null)
										undeterminedWriterR[1] = new BufferedWriter(new FileWriter(outputDir + "index_undetermined_R2.fastq"), BUFFER_SIZE);
								}
							}
							generatedUndeterminedFile = true;
						}
						
						//print to undetermined file
						writers1[writers1.length - 1].write(read.readF[0]);
						writers1[writers1.length - 1].newLine();
						writers1[writers1.length - 1].write(tempF1);
						writers1[writers1.length - 1].newLine();
						writers1[writers1.length - 1].write(description2);
						writers1[writers1.length - 1].newLine();
						writers1[writers1.length - 1].write(tempF2);
						writers1[writers1.length - 1].newLine();
						if(indexFile != null){
							writers3[writers3.length - 1].write(read.readIF[0]);
							writers3[writers3.length - 1].newLine();
							writers3[writers3.length - 1].write(read.readIF[1]);
							writers3[writers3.length - 1].newLine();
							writers3[writers3.length - 1].write(description2);
							writers3[writers3.length - 1].newLine();
							writers3[writers3.length - 1].write(read.readIF[3]);
							writers3[writers3.length - 1].newLine();
						}
						if(inputFile2 != null){
							undeterminedWriterR[0].write(read.readR[0]);
							undeterminedWriterR[0].newLine();
							undeterminedWriterR[0].write(tempR1);
							undeterminedWriterR[0].newLine();
							undeterminedWriterR[0].write(description2);
							undeterminedWriterR[0].newLine();
							undeterminedWriterR[0].write(tempR2);
							undeterminedWriterR[0].newLine();
							if(indexFile2 != null){
								undeterminedWriterR[1].write(read.readIR[0]);
								undeterminedWriterR[1].newLine();
								undeterminedWriterR[1].write(read.readIR[1]);
								undeterminedWriterR[1].newLine();
								undeterminedWriterR[1].write(description2);
								undeterminedWriterR[1].newLine();
								undeterminedWriterR[1].write(read.readIR[3]);
								undeterminedWriterR[1].newLine();
							}
						}
					}catch(Exception e){
						UtilMethods.defaultExceptionHandler(logWriter, e);
					}
				}
			}
		});
		
		//close everything
		reader1.close();
		if(indexFile != null)
			reader3.close();
		for(int i = 0; i < writers1.length; i++){
			if(writers1[i] != null)
				writers1[i].close();
			if(indexFile != null && writers3[i] != null)
				writers3[i].close();
		}
		if(inputFile2 != null){
			reader2.close();
			if(undeterminedWriterR[0] != null)
				undeterminedWriterR[0].close();
			if(indexFile2 != null){
				reader4.close();
				if(undeterminedWriterR[1] != null)
					undeterminedWriterR[1].close();
			}
			if(!mergePairedEnds){
				for(int i = 0; i < writers2.length; i++){
					if(writers2[i] != null)
						writers2[i].close();
					if(indexFile2 != null && writers4[i] != null)
						writers4[i].close();
				}
			}
		}
		
		//update stats
		for(int i = 0; i < sampleDNAF.size(); i++){
			EnumMap<Stat, String> map = stats.get(sampleMapF.get(sampleDNAF.get(i)));
			map.put(Stat.SEQUENCE_COUNT, DECIMAL_FORMAT.format(DNACounts[i][0].sum()));
			map.put(Stat.SEQUENCE_PERCENT, DECIMAL_FORMAT.format((double)DNACounts[i][0].sum() / (double)totalDNAProcessed.sum()));
			map.put(Stat.MERGED_COUNT, DECIMAL_FORMAT.format(DNACounts[i][2].sum()));
			map.put(Stat.MERGED_PERCENT, DECIMAL_FORMAT.format((double)DNACounts[i][2].sum() / (double)totalDNAProcessed.sum()));
			map.put(Stat.REMOVEADAPTER_COUNT, DECIMAL_FORMAT.format(DNACounts[i][3].sum()));
			map.put(Stat.REMOVEADAPTER_PERCENT, DECIMAL_FORMAT.format((double)DNACounts[i][3].sum() / (double)DNACounts[i][0].sum()));
			map.put(Stat.QUALITYTRIM_COUNT, DECIMAL_FORMAT.format(DNACounts[i][4].sum()));
			map.put(Stat.QUALITYTRIM_PERCENT, DECIMAL_FORMAT.format((double)DNACounts[i][4].sum() / (double)DNACounts[i][0].sum()));
			map.put(Stat.BASE_COUNT, DECIMAL_FORMAT.format(DNACounts[i][1].sum()));
			map.put(Stat.BASE_PERCENT, DECIMAL_FORMAT.format((double)DNACounts[i][1].sum() / (double)baseCount.sum()));
		}
		stats.get("Undetermined").put(Stat.SEQUENCE_COUNT, DECIMAL_FORMAT.format(DNACounts[DNACounts.length - 1][0].sum()));
		stats.get("Undetermined").put(Stat.SEQUENCE_PERCENT, DECIMAL_FORMAT.format((double)DNACounts[DNACounts.length - 1][0].sum() / (double)totalDNAProcessed.sum()));
		stats.get("Undetermined").put(Stat.BASE_COUNT, DECIMAL_FORMAT.format(DNACounts[DNACounts.length - 1][1].sum()));
		stats.get("Undetermined").put(Stat.BASE_PERCENT, DECIMAL_FORMAT.format((double)DNACounts[DNACounts.length - 1][1].sum() / (double)baseCount.sum()));
		
		return files;
	}
	
	//replaceOriginal = true will replace original file
	//replaceOriginal = false will generate new file in folder outDir with different name
	private long deduplicate(String readPath1, String readPath2, String indexPath1, String indexPath2, String outDir, boolean replaceOriginal, boolean standalone, String dupPath1, String dupPath2) throws Exception{
		long removed = 0L;
		
		File originalFile1 = new File(readPath1);
		int pos = originalFile1.getName().indexOf(".");
		String justName1 = pos > 0 ? originalFile1.getName().substring(0, pos) : originalFile1.getName();
		File outputFile1 = new File((replaceOriginal ? (originalFile1.getParentFile().getAbsolutePath() + File.separatorChar) : outputDir) + justName1 + "_dedup.fastq" + (outputGZIP ? ".gz" : ""));
		
		File originalFile3 = new File(indexPath1);
		pos = originalFile3.getName().indexOf(".");
		String justName3 = pos > 0 ? originalFile3.getName().substring(0, pos) : originalFile3.getName();
		File outputFile3 = new File((replaceOriginal ? (originalFile3.getParentFile().getAbsolutePath() + File.separatorChar) : outputDir) + justName3 + "_dedup.fastq" + (outputGZIP ? ".gz" : ""));
		
		File originalFile2 = null;
		String justName2 = null;
		File outputFile2 = null;
		
		File originalFile4 = null;
		String justName4 = null;
		File outputFile4 = null;
		
		if(readPath2 != null){
			originalFile2 = new File(readPath2);
			pos = originalFile2.getName().indexOf(".");
			justName2 = pos > 0 ? originalFile2.getName().substring(0, pos) : originalFile2.getName();
			outputFile2 = new File((replaceOriginal ? (originalFile2.getParentFile().getAbsolutePath() + File.separatorChar) : outputDir) + justName2 + "_dedup.fastq" + (outputGZIP ? ".gz" : ""));
			
			originalFile4 = new File(indexPath2);
			pos = originalFile4.getName().indexOf(".");
			justName4 = pos > 0 ? originalFile4.getName().substring(0, pos) : originalFile4.getName();
			outputFile4 = new File((replaceOriginal ? (originalFile4.getParentFile().getAbsolutePath() + File.separatorChar) : outputDir) + justName4 + "_dedup.fastq" + (outputGZIP ? ".gz" : ""));
		}
		
		BufferedReader reader1;
		BufferedReader reader2 = null;
		BufferedReader reader3;
		BufferedReader reader4 = null;
		if(!standalone || inputGZIP){
			reader1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(readPath1), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
			reader3 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(indexPath1), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
			if(readPath2 != null){
				reader2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(readPath2), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
				reader4 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(indexPath2), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
			}
		}else{
			reader1 = new BufferedReader(new FileReader(readPath1), BUFFER_SIZE);
			reader3 = new BufferedReader(new FileReader(indexPath1), BUFFER_SIZE);
			if(readPath2 != null){
				reader2 = new BufferedReader(new FileReader(readPath2), BUFFER_SIZE);
				reader4 = new BufferedReader(new FileReader(indexPath2), BUFFER_SIZE);
			}
		}
		
		if(removeFirstDup){
			BufferedWriter writer1;
			BufferedWriter writer2 = null;
			BufferedWriter writer3;
			BufferedWriter writer4 = null;
			if(outputGZIP){
				writer1 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile1.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
				writer3 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile3.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
				if(readPath2 != null){
					writer2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile2.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
					writer4 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile4.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
				}
			}else{
				writer1 = new BufferedWriter(new FileWriter(outputFile1.getAbsolutePath()), BUFFER_SIZE);
				writer3 = new BufferedWriter(new FileWriter(outputFile3.getAbsolutePath()), BUFFER_SIZE);
				if(readPath2 != null){
					writer2 = new BufferedWriter(new FileWriter(outputFile2.getAbsolutePath()), BUFFER_SIZE);
					writer4 = new BufferedWriter(new FileWriter(outputFile4.getAbsolutePath()), BUFFER_SIZE);
				}
			}
			BufferedWriter dupWriter1 = null;
			BufferedWriter dupWriter2 = null;
			HashSet<BitSet> set = new HashSet<BitSet>();
			String[] lines1 = new String[4];
			String[] lines2 = new String[4];
			String[] lines3 = new String[4];
			String[] lines4 = new String[4];
			while((lines1[0] = reader1.readLine()) != null && (lines1[1] = reader1.readLine()) != null &&
					(lines1[2] = reader1.readLine()) != null && (lines1[3] = reader1.readLine()) != null &&
					(lines3[0] = reader3.readLine()) != null && (lines3[1] = reader3.readLine()) != null &&
					(lines3[2] = reader3.readLine()) != null && (lines3[3] = reader3.readLine()) != null &&
					(readPath2 == null || ((lines2[0] = reader2.readLine()) != null && (lines2[1] = reader2.readLine()) != null &&
					(lines2[2] = reader2.readLine()) != null && (lines2[3] = reader2.readLine()) != null &&
					(lines4[0] = reader4.readLine()) != null && (lines4[1] = reader4.readLine()) != null &&
					(lines4[2] = reader4.readLine()) != null && (lines4[3] = reader4.readLine()) != null))){
				if(standalone){
					totalDNAProcessed.add(readPath2 == null ? 1 : 2);
					baseCount.add(lines1[1].length() + (readPath2 == null ? 0 : lines2[1].length()));
					if(!parallel && totalDNAProcessed.sum() % printProcessedInterval == 0){
						logWriter.println("Reads Processed So Far: " + DECIMAL_FORMAT.format(totalDNAProcessed.sum()));
						logWriter.println("Run Time So Far: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
						logWriter.println("Current Memory Usage (GB): " + DECIMAL_FORMAT.format(UtilMethods.currentMemoryUsage()));
						logWriter.flush();
					}
				}
				if(lines3[1].length() < randUMILength)
					continue;
				BitSet key = UtilMethods.toBit(lines3[1].substring(lines3[1].length() - randUMILength, lines3[1].length()));
				if(set.contains(key)){ //if a read with the same UMI has already been encountered, then the current read is removed
					duplicatesRemoved += readPath2 == null ? 1 : 2;
					removed += readPath2 == null ? 1 : 2;
					if(duplicatesRemoved % printDuplicateInterval == 0){
						logWriter.println("Deduplicated So Far: " + DECIMAL_FORMAT.format(duplicatesRemoved));
						logWriter.println("Run Time So Far: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
						logWriter.println("Current Memory Usage (GB): " + DECIMAL_FORMAT.format(UtilMethods.currentMemoryUsage()));
						logWriter.flush();
					}
					if(saveDup){ //print removed duplicate to another file if needed
						if(dupWriter1 == null){ //initialize duplicate writer if it is not initialized
							dupWriter1 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(dupPath1), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
							dupWriter2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(dupPath2), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
							generatedDupFiles += 2;
						}
						dupWriter1.write(lines1[0]);
						dupWriter1.newLine();
						dupWriter1.write(lines1[1]);
						dupWriter1.newLine();
						dupWriter1.write(description2);
						dupWriter1.newLine();
						dupWriter1.write(lines1[3]);
						dupWriter1.newLine();
						
						dupWriter2.write(lines3[0]);
						dupWriter2.newLine();
						dupWriter2.write(lines3[1]);
						dupWriter2.newLine();
						dupWriter2.write(description2);
						dupWriter2.newLine();
						dupWriter2.write(lines3[3]);
						dupWriter2.newLine();
						
						if(readPath2 != null){
							dupWriter1.write(lines2[0]);
							dupWriter1.newLine();
							dupWriter1.write(lines2[1]);
							dupWriter1.newLine();
							dupWriter1.write(description2);
							dupWriter1.newLine();
							dupWriter1.write(lines2[3]);
							dupWriter1.newLine();
							
							dupWriter2.write(lines4[0]);
							dupWriter2.newLine();
							dupWriter2.write(lines4[1]);
							dupWriter2.newLine();
							dupWriter2.write(description2);
							dupWriter2.newLine();
							dupWriter2.write(lines4[3]);
							dupWriter2.newLine();
						}
					}
				}else{ //if the UMI is not encountered before, then the read is saved
					writer1.write(lines1[0]);
					writer1.newLine();
					writer1.write(lines1[1]);
					writer1.newLine();
					writer1.write(description2);
					writer1.newLine();
					writer1.write(lines1[3]);
					writer1.newLine();
					
					writer3.write(lines3[0]);
					writer3.newLine();
					writer3.write(lines3[1]);
					writer3.newLine();
					writer3.write(description2);
					writer3.newLine();
					writer3.write(lines3[3]);
					writer3.newLine();
					
					if(readPath2 != null){
						writer2.write(lines2[0]);
						writer2.newLine();
						writer2.write(lines2[1]);
						writer2.newLine();
						writer2.write(description2);
						writer2.newLine();
						writer2.write(lines2[3]);
						writer2.newLine();
						
						writer4.write(lines4[0]);
						writer4.newLine();
						writer4.write(lines4[1]);
						writer4.newLine();
						writer4.write(description2);
						writer4.newLine();
						writer4.write(lines4[3]);
						writer4.newLine();
					}
					
					set.add(key);
				}
			}
			writer1.close();
			writer3.close();
			if(readPath2 != null){
				writer2.close();
				writer4.close();
			}
			if(dupWriter1 != null){
				dupWriter1.close();
				dupWriter2.close();
			}
		}else if(removeBestDup){
			BufferedWriter dupWriter1 = null;
			BufferedWriter dupWriter2 = null;
			
			HashMap<BitSet, CompressedRead> map = new HashMap<BitSet, CompressedRead>();
			String[] lines1 = new String[4];
			String[] lines2 = new String[4];
			String[] lines3 = new String[4];
			String[] lines4 = new String[4];
			while((lines1[0] = reader1.readLine()) != null && (lines1[1] = reader1.readLine()) != null &&
					(lines1[2] = reader1.readLine()) != null && (lines1[3] = reader1.readLine()) != null &&
					(lines3[0] = reader3.readLine()) != null && (lines3[1] = reader3.readLine()) != null &&
					(lines3[2] = reader3.readLine()) != null && (lines3[3] = reader3.readLine()) != null &&
					(readPath2 == null || ((lines2[0] = reader2.readLine()) != null && (lines2[1] = reader2.readLine()) != null &&
					(lines2[2] = reader2.readLine()) != null && (lines2[3] = reader2.readLine()) != null &&
					(lines4[0] = reader4.readLine()) != null && (lines4[1] = reader4.readLine()) != null &&
					(lines4[2] = reader4.readLine()) != null && (lines4[3] = reader4.readLine()) != null))){
				if(standalone){
					totalDNAProcessed.add(readPath2 == null ? 1 : 2);
					baseCount.add(lines1[1].length() + (readPath2 == null ? 0 : lines2[1].length()));
					if(!parallel && totalDNAProcessed.sum() % printProcessedInterval == 0){
						logWriter.println("Reads Processed So Far: " + DECIMAL_FORMAT.format(totalDNAProcessed.sum()));
						logWriter.println("Run Time So Far: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
						logWriter.println("Current Memory Usage (GB): " + DECIMAL_FORMAT.format(UtilMethods.currentMemoryUsage()));
						logWriter.flush();
					}
				}
				if(lines3[1].length() < randUMILength)
					continue;
				BitSet key = UtilMethods.toBit(lines3[1].substring(lines3[1].length() - randUMILength, lines3[1].length()));
				double currQuality = UtilMethods.toError(lines1[3] + (readPath2 == null ? "" : lines2[3]), 0);
				if(map.containsKey(key)){
					duplicatesRemoved += readPath2 == null ? 1 : 2;
					removed += readPath2 == null ? 1 : 2;
					if(duplicatesRemoved % printDuplicateInterval == 0){
						logWriter.println("Deduplicated So Far: " + DECIMAL_FORMAT.format(duplicatesRemoved));
						logWriter.println("Run Time So Far: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
						logWriter.println("Current Memory Usage (GB): " + DECIMAL_FORMAT.format(UtilMethods.currentMemoryUsage()));
						logWriter.flush();
					}
					CompressedRead other = map.get(key);
					if(other.error > currQuality){ //if current read quality is better
						if(saveDup){
							if(dupWriter1 == null){ //initialize writer for writing duplicates
								dupWriter1 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(dupPath1), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
								dupWriter2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(dupPath2), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
								generatedDupFiles += 2;
							}
							//print other read to duplicates
							dupWriter1.write(other.descriptionFRead);
							dupWriter1.newLine();
							dupWriter1.write(UtilMethods.toSequence(other.sequenceFRead));
							dupWriter1.newLine();
							dupWriter1.write(description2);
							dupWriter1.newLine();
							dupWriter1.write(UtilMethods.byteArrayToQScore(other.qualityFRead));
							dupWriter1.newLine();
							
							dupWriter2.write(other.descriptionFIndex);
							dupWriter2.newLine();
							dupWriter2.write(UtilMethods.toSequence(other.sequenceFIndex));
							dupWriter2.newLine();
							dupWriter2.write(description2);
							dupWriter2.newLine();
							dupWriter2.write(UtilMethods.byteArrayToQScore(other.qualityFIndex));
							dupWriter2.newLine();
							
							if(readPath2 != null){
								dupWriter1.write(other.descriptionRRead);
								dupWriter1.newLine();
								dupWriter1.write(UtilMethods.toSequence(other.sequenceRRead));
								dupWriter1.newLine();
								dupWriter1.write(description2);
								dupWriter1.newLine();
								dupWriter1.write(UtilMethods.byteArrayToQScore(other.qualityRRead));
								dupWriter1.newLine();
								
								dupWriter2.write(other.descriptionRIndex);
								dupWriter2.newLine();
								dupWriter2.write(UtilMethods.toSequence(other.sequenceRIndex));
								dupWriter2.newLine();
								dupWriter2.write(description2);
								dupWriter2.newLine();
								dupWriter2.write(UtilMethods.byteArrayToQScore(other.qualityRIndex));
								dupWriter2.newLine();
							}
						}
						//save current read in HashMap
						if(readPath2 == null)
							map.put(key, new CompressedRead(currQuality, lines1[0], lines1[1], lines1[3], lines3[0], lines3[1], lines3[3]));
						else
							map.put(key, new CompressedRead(currQuality, lines1[0], lines1[1], lines1[3], lines3[0], lines3[1], lines3[3], lines2[0], lines2[1], lines2[3], lines4[0], lines4[1], lines4[3]));
					}else{ //remove current read
						if(saveDup){
							if(dupWriter1 == null){
								dupWriter1 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(dupPath1), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
								dupWriter2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(dupPath2), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
								generatedDupFiles += 2;
							}
							dupWriter1.write(lines1[0]);
							dupWriter1.newLine();
							dupWriter1.write(lines1[1]);
							dupWriter1.newLine();
							dupWriter1.write(description2);
							dupWriter1.newLine();
							dupWriter1.write(lines1[3]);
							dupWriter1.newLine();
							
							dupWriter2.write(lines3[0]);
							dupWriter2.newLine();
							dupWriter2.write(lines3[1]);
							dupWriter2.newLine();
							dupWriter2.write(description2);
							dupWriter2.newLine();
							dupWriter2.write(lines3[3]);
							dupWriter2.newLine();
							
							if(readPath2 != null){
								dupWriter1.write(lines2[0]);
								dupWriter1.newLine();
								dupWriter1.write(lines2[1]);
								dupWriter1.newLine();
								dupWriter1.write(description2);
								dupWriter1.newLine();
								dupWriter1.write(lines2[3]);
								dupWriter1.newLine();
								
								dupWriter2.write(lines4[0]);
								dupWriter2.newLine();
								dupWriter2.write(lines4[1]);
								dupWriter2.newLine();
								dupWriter2.write(description2);
								dupWriter2.newLine();
								dupWriter2.write(lines4[3]);
								dupWriter2.newLine();
							}
						}
					}
				}else{ //if the current UMI was not encountered before
					if(readPath2 == null)
						map.put(key, new CompressedRead(currQuality, lines1[0], lines1[1], lines1[3], lines3[0], lines3[1], lines3[3]));
					else
						map.put(key, new CompressedRead(currQuality, lines1[0], lines1[1], lines1[3], lines3[0], lines3[1], lines3[3], lines2[0], lines2[1], lines2[3], lines4[0], lines4[1], lines4[3]));
				}
			}
			reader1.close();
			reader3.close();
			if(readPath2 != null){
				reader2.close();
				reader4.close();
			}
			
			BufferedWriter writer1;
			BufferedWriter writer2 = null;
			BufferedWriter writer3;
			BufferedWriter writer4 = null;
			if(outputGZIP){
				writer1 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile1.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
				writer3 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile3.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
				if(readPath2 != null){
					writer2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile2.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
					writer4 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile4.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
				}
			}else{
				writer1 = new BufferedWriter(new FileWriter(outputFile1.getAbsolutePath()), BUFFER_SIZE);
				writer3 = new BufferedWriter(new FileWriter(outputFile3.getAbsolutePath()), BUFFER_SIZE);
				if(readPath2 != null){
					writer2 = new BufferedWriter(new FileWriter(outputFile2.getAbsolutePath()), BUFFER_SIZE);
					writer4 = new BufferedWriter(new FileWriter(outputFile4.getAbsolutePath()), BUFFER_SIZE);
				}
			}
			
			//go through the best reads and print them
			for(CompressedRead val : map.values()){
				writer1.write(val.descriptionFRead);
				writer1.newLine();
				writer1.write(UtilMethods.toSequence(val.sequenceFRead));
				writer1.newLine();
				writer1.write(description2);
				writer1.newLine();
				writer1.write(UtilMethods.byteArrayToQScore(val.qualityFRead));
				writer1.newLine();
				
				writer3.write(val.descriptionFIndex);
				writer3.newLine();
				writer3.write(UtilMethods.toSequence(val.sequenceFIndex));
				writer3.newLine();
				writer3.write(description2);
				writer3.newLine();
				writer3.write(UtilMethods.byteArrayToQScore(val.qualityFIndex));
				writer3.newLine();
				
				if(readPath2 != null){
					writer2.write(val.descriptionRRead);
					writer2.newLine();
					writer2.write(UtilMethods.toSequence(val.sequenceRRead));
					writer2.newLine();
					writer2.write(description2);
					writer2.newLine();
					writer2.write(UtilMethods.byteArrayToQScore(val.qualityRRead));
					writer2.newLine();
					
					writer4.write(val.descriptionRIndex);
					writer4.newLine();
					writer4.write(UtilMethods.toSequence(val.sequenceRIndex));
					writer4.newLine();
					writer4.write(description2);
					writer4.newLine();
					writer4.write(UtilMethods.byteArrayToQScore(val.qualityRIndex));
					writer4.newLine();
				}
			}
			writer1.close();
			writer3.close();
			if(readPath2 != null){
				writer2.close();
				writer4.close();
			}
		}
		
		if(replaceOriginal){ //replace the original files
			originalFile1.delete();
			originalFile3.delete();
			outputFile1.renameTo(new File(originalFile1.getParentFile().getAbsolutePath() + File.separatorChar + justName1 + (outputGZIP ? ".fastq.gz" : ".fastq")));
			outputFile3.renameTo(new File(originalFile3.getParentFile().getAbsolutePath() + File.separatorChar + justName3 + (outputGZIP ? ".fastq.gz" : ".fastq")));
			if(readPath2 != null){
				originalFile2.delete();
				originalFile4.delete();
				outputFile2.renameTo(new File(originalFile2.getParentFile().getAbsolutePath() + File.separatorChar + justName2 + (outputGZIP ? ".fastq.gz" : ".fastq")));
				outputFile4.renameTo(new File(originalFile4.getParentFile().getAbsolutePath() + File.separatorChar + justName4 + (outputGZIP ? ".fastq.gz" : ".fastq")));
			}
		}
		
		return removed;
	}
	
	private void pairMergeFiles() throws Exception{
		BufferedReader reader1;
		BufferedReader reader2;
		if(inputGZIP){
			reader1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFile.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
			reader2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFile2.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
		}else{
			reader1 = new BufferedReader(new FileReader(inputFile), BUFFER_SIZE);
			reader2 = new BufferedReader(new FileReader(inputFile2), BUFFER_SIZE);
		}
		BufferedWriter writer;
		int pos = inputFile.getName().indexOf(".");
		String justName = pos > 0 ? inputFile.getName().substring(0, pos) : inputFile.getName();
		if(outputGZIP){
			writer = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + justName + "_pairmerged.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
		}else{
			writer = new BufferedWriter(new FileWriter(outputDir + justName + "_pairmerged.fastq"), BUFFER_SIZE);
		}
		
		Object lock = new Object();
		
		StreamSupport.stream(new ReadSpliterator<Read>(splitBatchSize, reader1, reader2, null, null), parallel).forEach((read) -> {
			totalDNAProcessed.add(2);
			baseCount.add(read.readF[1].length() + read.readR[1].length());
			if(!parallel && totalDNAProcessed.sum() % printProcessedInterval == 0){
				logWriter.println("Reads Processed So Far: " + DECIMAL_FORMAT.format(totalDNAProcessed.sum()));
				logWriter.println("Run Time So Far: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
				logWriter.println("Current Memory Usage (GB): " + DECIMAL_FORMAT.format(UtilMethods.currentMemoryUsage()));
				logWriter.flush();
			}
			//merge the two lines
			String[] merged = UtilMethods.mergeReads(read.readF[1], read.readF[3], read.readR[1], read.readR[3], editMaxM, probM, wildcard);
			if(!removeNoMergeReads || merged[0].length() != read.readF[1].length() + read.readR[1].length()){
				if(merged[0].length() != read.readF[1].length() + read.readR[1].length())
					totalReadsMerged.add(2);
				
				synchronized(lock){
					try{
						//print to file
						writer.write(read.readF[0]);
						writer.newLine();
						writer.write(merged[0]);
						writer.newLine();
						writer.write(description2);
						writer.newLine();
						writer.write(merged[1]);
						writer.newLine();
					}catch(Exception e){
						UtilMethods.defaultExceptionHandler(logWriter, e);
					}
				}
			}else{
				undeterminedDNA.add(2);
			}
		});
		
		reader1.close();
		reader2.close();
		writer.close();
	}
	
	private void filterFile() throws Exception{
		BufferedReader reader;
		if(inputGZIP){
			reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFile.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
		}else{
			reader = new BufferedReader(new FileReader(inputFile), BUFFER_SIZE);
		}
		BufferedWriter writer;
		File outputFile;
		int pos = inputFile.getName().indexOf(".");
		String justName = pos > 0 ? inputFile.getName().substring(0, pos) : inputFile.getName();
		if(outputGZIP){
			writer = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + justName + "_filtered.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
			outputFile = new File(outputDir + justName + "_filtered.fastq.gz");
		}else{
			writer = new BufferedWriter(new FileWriter(outputDir + justName + "_filtered.fastq"), BUFFER_SIZE);
			outputFile = new File(outputDir + justName + "_filtered.fastq");
		}
		
		Object lock = new Object();
		
		StreamSupport.stream(new ReadSpliterator<Read>(splitBatchSize, reader, null, null, null), parallel).forEach((read) -> {
			totalDNAProcessed.add(1);
			baseCount.add(read.readF[1].length());
			if(!parallel && totalDNAProcessed.sum() % printProcessedInterval == 0){
				logWriter.println("Reads Processed So Far: " + DECIMAL_FORMAT.format(totalDNAProcessed.sum()));
				logWriter.println("Run Time So Far: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
				logWriter.println("Current Memory Usage (GB): " + DECIMAL_FORMAT.format(UtilMethods.currentMemoryUsage()));
				logWriter.flush();
			}
			
			//check if quality is good enough
			boolean qualityAcceptable;
			if(filterAlgorithm){
				qualityAcceptable = UtilMethods.toError(read.readF[3], 0) <= qualityFilter;
			}else{
				qualityAcceptable = UtilMethods.toQScore(read.readF[3], 0) >= qualityFilter;
			}
			if(qualityAcceptable && (removeDNAWithNPercent > 1.0 || (removeDNAWithNPercent >= 0.0 && UtilMethods.percentN(read.readF[1]) <= removeDNAWithNPercent) || (removeDNAWithNPercent < 0.0 && UtilMethods.countN(read.readF[1]) <= 0))){
				String[] temp;
				String[] qualityTrimmed;
				
				//trim 'N' bp
				temp = UtilMethods.trimN(read.readF[1], read.readF[3], trimNPercent);
				read.readF[1] = temp[0];
				read.readF[3] = temp[1];
				
				int trimmedQuality = 0;
				int removedAdapter = 0;
				
				//quality trim
				if(trimAlgorithm){
					temp = UtilMethods.qualityTrim2(read.readF[1], read.readF[3], qualityTrimQScore1, true, qualityTrimLength);
					qualityTrimmed = UtilMethods.qualityTrim2(temp[0], temp[1], qualityTrimQScore2, false, qualityTrimLength);
				}else{
					temp = UtilMethods.qualityTrim1(read.readF[1], read.readF[3], qualityTrimQScore1, true, qualityTrimLength);
					qualityTrimmed = UtilMethods.qualityTrim1(temp[0], temp[1], qualityTrimQScore2, false, qualityTrimLength);
				}
				if(read.readF[1].length() != qualityTrimmed[0].length()){
					trimmedQuality++;
				}
				//remove adapters
				String[] removedAdapters = UtilMethods.removeAdapters(qualityTrimmed[0], qualityTrimmed[1], adaptersF, editMaxA, minOverlapA, maxOffsetA, allowIndelsA, probA, wildcard, adapterPatternF);
				if(qualityTrimmed[0].length() != removedAdapters[0].length()){
					removedAdapter++;
				}
				read.readF[1] = removedAdapters[0];
				read.readF[3] = removedAdapters[1];
				
				if((!removeUntrimmedReads || trimmedQuality == 1) && (!removeNoAdapterReads || removedAdapter == 1) &&
						read.readF[1].length() >= minLength && read.readF[1].length() <= maxLength){
					totalQualityTrimmed.add(trimmedQuality);
					totalRemovedAdapters.add(removedAdapter);
					
					synchronized(lock){
						try{
							//print to file
							writer.write(read.readF[0]);
							writer.newLine();
							writer.write(read.readF[1]);
							writer.newLine();
							writer.write(description2);
							writer.newLine();
							writer.write(read.readF[3]);
							writer.newLine();
						}catch(Exception e){
							UtilMethods.defaultExceptionHandler(logWriter, e);
						}
					}
				}else{
					undeterminedDNA.add(1);
				}
			}else{
				undeterminedDNA.add(1);
			}
		});
		
		reader.close();
		writer.close();
		
		if(replaceOriginal){
			inputFile.delete();
			outputFile.renameTo(new File(inputFile.getParentFile().getAbsolutePath() + File.separatorChar + justName + (outputGZIP ? ".fastq.gz" : ".fastq")));
		}
	}
	
	public static void main(String[] args) throws Exception{
		Thread.setDefaultUncaughtExceptionHandler(new Thread.UncaughtExceptionHandler(){
			@Override
			public void uncaughtException(Thread t, Throwable e){ //default exception handler
				UtilMethods.defaultExceptionHandler(logWriter, e);
			}
		});
		
		//read commands and parse them, then call the constructor to start the processing
		//updates the global variables that contain all the options
		if(args[0].equals("--demultiplex") || args[0].equals("-d")){
			boolean clearDir = false;
			for(int i = 1; i < args.length; i++){
				if(args[i].equals("-dF")){
					removeFirstDup = true;
					removeBestDup = false;
				}else if(args[i].equals("-dB")){
					removeBestDup = true;
					removeFirstDup = false;
				}else if(args[i].equals("--umi")){
					randUMILength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-N")){
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						removeDNAWithNPercent = Double.parseDouble(args[++i]);
					}else{
						removeDNAWithNPercent = -1.0;
					}
				}else if(args[i].equals("-fB")){
					maxOffsetB = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-fA")){
					maxOffsetA = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-Q")){
					qualityFilter = Double.parseDouble(args[++i]);
				}else if(args[i].equals("-r")){
					inputFile = new File(args[++i]);
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						inputFile2 = new File(args[++i]);
					}
				}else if(args[i].equals("-i")){
					indexFile = new File(args[++i]);
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						indexFile2 = new File(args[++i]);
					}
					randUMILength = 12;
				}else if(args[i].equals("-s")){
					sampleFile = new File(args[++i]);
				}else if(args[i].equals("-o")){
					outputDir = args[++i];
					if(!outputDir.endsWith(File.separator) && !outputDir.endsWith("/") && !outputDir.endsWith("\\")){
						outputDir += File.separator;
					}
				}else if(args[i].equals("-gz")){
					outputGZIP = true;
				}else if(args[i].equals("--cleardir")){
					clearDir = true;
				}else if(args[i].equals("--savetemp")){
					saveTemp = true;
				}else if(args[i].equals("-kB")){
					removeBarRand = false;
				}else if(args[i].equals("-kE")){
					removeEnzyme = false;
				}else if(args[i].equals("-eA")){
					double d = Double.parseDouble(args[++i]);
					editMaxA = d < 1.0 ? -d : d;
				}else if(args[i].equals("-eB")){
					double d = Double.parseDouble(args[++i]);
					editMaxB = d < 1.0 ? -d : d;
				}else if(args[i].equals("-eM")){
					double d = Double.parseDouble(args[++i]);
					editMaxM = d < 1.0 ? -d : d;
				}else if(args[i].equals("--savedup")){
					saveDup = true;
				}else if(args[i].equals("--printprocessed")){
					printProcessedInterval = Long.parseLong(args[++i]);
				}else if(args[i].equals("--printduplicate")){
					printDuplicateInterval = Long.parseLong(args[++i]);
				}else if(args[i].equals("-m")){
					mergePairedEnds = true;
				}else if(args[i].equals("--min")){
					minLength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--max")){
					maxLength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-a")){
					while(i + 1 < args.length && !args[i + 1].startsWith("-")){
						if(args[i + 1].startsWith("^")){
							adaptersF.add(new Adapter(args[i + 1].substring(1), true, true));
						}else{
							adaptersF.add(new Adapter(args[i + 1], true, false));
						}
						i++;
					}
				}else if(args[i].equals("-A")){
					while(i + 1 < args.length && !args[i + 1].startsWith("-")){
						if(args[i + 1].endsWith("$")){
							adaptersF.add(new Adapter(args[i + 1].substring(0, args[i + 1].length() - 1), false, true));
						}else{
							adaptersF.add(new Adapter(args[i + 1], false, false));
						}
						i++;
					}
				}else if(args[i].equals("-z")){
					while(i + 1 < args.length && !args[i + 1].startsWith("-")){
						if(args[i + 1].startsWith("^")){
							adaptersR.add(new Adapter(args[i + 1].substring(1), true, true));
						}else{
							adaptersR.add(new Adapter(args[i + 1], true, false));
						}
						i++;
					}
				}else if(args[i].equals("-Z")){
					while(i + 1 < args.length && !args[i + 1].startsWith("-")){
						if(args[i + 1].endsWith("$")){
							adaptersR.add(new Adapter(args[i + 1].substring(0, args[i + 1].length() - 1), false, true));
						}else{
							adaptersR.add(new Adapter(args[i + 1], false, false));
						}
						i++;
					}
				}else if(args[i].equals("-oA")){
					minOverlapA = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-oB")){
					minOverlapB = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--altqtrim")){
					trimAlgorithm = true;
					qualityTrimLength = Integer.MAX_VALUE;
				}else if(args[i].equals("-q")){
					if(i + 2 >= args.length || args[i + 2].startsWith("-")){
						qualityTrimQScore2 = Integer.parseInt(args[++i]);
					}else{
						qualityTrimQScore1 = Integer.parseInt(args[++i]);
						qualityTrimQScore2 = Integer.parseInt(args[++i]);
					}
				}else if(args[i].equals("--qtrimlen")){
					qualityTrimLength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-n")){
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						trimNPercent = Double.parseDouble(args[++i]);
					}else{
						trimNPercent = 0.5;
					}
				}else if(args[i].equals("-iB")){
					allowIndelsB = true;
				}else if(args[i].equals("-iA")){
					allowIndelsA = true;
				}else if(args[i].equals("--altqfilter")){
					filterAlgorithm = true;
				}else if(args[i].equals("-R")){
					checkReversedReads = true;
				}else if(args[i].equals("-w")){
					wildcard = true;
				}else if(args[i].equals("-S")){
					singleBarcodeMatchOnly = true;
				}else if(args[i].equals("-nQ")){
					removeUntrimmedReads = true;
				}else if(args[i].equals("-nA")){
					removeNoAdapterReads = true;
				}else if(args[i].equals("-nM")){
					removeNoMergeReads = true;
				}else if(args[i].equals("-p")){
					parallel = true;
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						splitBatchSize = Integer.parseInt(args[++i]);
					}
				}else if(args[i].equals("-pB")){
					probB = 0.5;
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						probB = Double.parseDouble(args[++i]);
					}
				}else if(args[i].equals("-pA")){
					probA = 0.5;
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						probA = Double.parseDouble(args[++i]);
					}
				}else if(args[i].equals("-pM")){
					probM = 0.5;
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						probM = Double.parseDouble(args[++i]);
					}
				}else{
					throw new Exception("Command not supported: " + args[i]);
				}
			}
			boolean isDirClear = false;
			if(outputDir == null){
				outputDir = inputFile.getParent();
				if(!outputDir.endsWith(File.separator) && !outputDir.endsWith("/") && !outputDir.endsWith("\\")){
					outputDir += File.separator;
				}
			}else{
				File f = new File(outputDir);
				if(f.exists() && clearDir){
					UtilMethods.deleteFolder(f);
					isDirClear = true;
				}
				if(!f.exists())
					f.mkdirs();
			}
			File f2 = new File(outputDir + "temp" + File.separatorChar);
			if(!f2.exists())
				f2.mkdirs();
			if(saveDup && (removeBestDup || removeFirstDup)){
				File f3 = new File(outputDir + "dup" + File.separatorChar);
				if(!f3.exists())
					f3.mkdirs();
			}
			if(inputFile.getAbsolutePath().toLowerCase().endsWith(".gz") || inputFile.getAbsolutePath().toLowerCase().endsWith(".gzip")){
				inputGZIP = true;
			}
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outputDir + "FastQParse_" + Mode.DEMULTIPLEX.description2 + ".log"), BUFFER_SIZE_LOG));
			Date date = new Date();
			logWriter.println("Started on: " + DATE_FORMAT.format(date));
			System.out.println("Started on: " + DATE_FORMAT.format(date));
			if(isDirClear){
				logWriter.println("Cleared Directory: " + outputDir);
			}
			new FastQParseMain(Mode.DEMULTIPLEX);
			
			date = new Date();
			logWriter.println("Ended on: " + DATE_FORMAT.format(date));
			System.out.println("Ended on: " + DATE_FORMAT.format(date));
		}else if(args[0].equals("--dedup") || args[0].equals("-D")){
			for(int i = 1; i < args.length; i++){
				if(args[i].equals("-dF")){
					removeFirstDup = true;
					removeBestDup = false;
				}else if(args[i].equals("-dB")){
					removeBestDup = true;
					removeFirstDup = false;
				}else if(args[i].equals("-r")){
					inputFile = new File(args[++i]);
				}else if(args[i].equals("-i")){
					indexFile = new File(args[++i]);
					randUMILength = 12;
				}else if(args[i].equals("-o")){
					outputDir = args[++i];
					if(!outputDir.endsWith(File.separator) && !outputDir.endsWith("/") && !outputDir.endsWith("\\")){
						outputDir += File.separator;
					}
				}else if(args[i].equals("-O")){
					replaceOriginal = true;
				}else if(args[i].equals("-gz")){
					outputGZIP = true;
				}else if(args[i].equals("--savedup")){
					saveDup = true;
				}else if(args[i].equals("--umi")){
					randUMILength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--printprocessed")){
					printProcessedInterval = Long.parseLong(args[++i]);
				}else if(args[i].equals("--printduplicate")){
					printDuplicateInterval = Long.parseLong(args[++i]);
				}else{
					throw new Exception("Command not supported: " + args[i]);
				}
			}
			if(outputDir == null){
				outputDir = inputFile.getParent();
				if(!outputDir.endsWith(File.separator) && !outputDir.endsWith("/") && !outputDir.endsWith("\\")){
					outputDir += File.separator;
				}
			}else{
				File f = new File(outputDir);
				if(!f.exists())
					f.mkdirs();
			}
			if(inputFile.getAbsolutePath().toLowerCase().endsWith(".gz") || inputFile.getAbsolutePath().toLowerCase().endsWith(".gzip")){
				inputGZIP = true;
			}
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outputDir + "FastQParse_" + Mode.DEDUP.description2 + ".log"), BUFFER_SIZE_LOG));
			Date date = new Date();
			logWriter.println("Started on: " + DATE_FORMAT.format(date));
			System.out.println("Started on: " + DATE_FORMAT.format(date));
			if(!removeFirstDup && !removeBestDup)
				removeFirstDup = true;
			new FastQParseMain(Mode.DEDUP);
			
			date = new Date();
			logWriter.println("Ended on: " + DATE_FORMAT.format(date));
			System.out.println("Ended on: " + DATE_FORMAT.format(date));
		}else if(args[0].equals("--pairmerge") || args[0].equals("-m")){
			for(int i = 1; i < args.length; i++){
				if(args[i].equals("-r")){
					inputFile = new File(args[++i]);
					inputFile2 = new File(args[++i]);
				}else if(args[i].equals("-o")){
					outputDir = args[++i];
					if(!outputDir.endsWith(File.separator) && !outputDir.endsWith("/") && !outputDir.endsWith("\\")){
						outputDir += File.separator;
					}
				}else if(args[i].equals("-gz")){
					outputGZIP = true;
				}else if(args[i].equals("-eM")){
					double d = Double.parseDouble(args[++i]);
					editMaxM = d < 1.0 ? -d : d;
				}else if(args[i].equals("--printprocessed")){
					printProcessedInterval = Long.parseLong(args[++i]);
				}else if(args[i].equals("-w")){
					wildcard = true;
				}else if(args[i].equals("-nM")){
					removeNoMergeReads = true;
				}else if(args[i].equals("-p")){
					parallel = true;
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						splitBatchSize = Integer.parseInt(args[++i]);
					}
				}else if(args[i].equals("-pM")){
					probM = 0.5;
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						probM = Double.parseDouble(args[++i]);
					}
				}else{
					throw new Exception("Command not supported: " + args[i]);
				}
			}
			if(outputDir == null){
				outputDir = inputFile.getParent();
				if(!outputDir.endsWith(File.separator) && !outputDir.endsWith("/") && !outputDir.endsWith("\\")){
					outputDir += File.separator;
				}
			}else{
				File f = new File(outputDir);
				if(!f.exists())
					f.mkdirs();
			}
			if(inputFile.getAbsolutePath().toLowerCase().endsWith(".gz") || inputFile.getAbsolutePath().toLowerCase().endsWith(".gzip")){
				inputGZIP = true;
			}
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outputDir + "FastQParse_" + Mode.PAIRMERGE.description2 + ".log"), BUFFER_SIZE_LOG));
			Date date = new Date();
			logWriter.println("Started on: " + DATE_FORMAT.format(date));
			System.out.println("Started on: " + DATE_FORMAT.format(date));
			new FastQParseMain(Mode.PAIRMERGE);
			
			date = new Date();
			logWriter.println("Ended on: " + DATE_FORMAT.format(date));
			System.out.println("Ended on: " + DATE_FORMAT.format(date));
		}else if(args[0].equals("--filter") || args[0].equals("-f")){
			for(int i = 1; i < args.length; i++){
				if(args[i].equals("-r")){
					inputFile = new File(args[++i]);
				}else if(args[i].equals("-o")){
					outputDir = args[++i];
					if(!outputDir.endsWith(File.separator) && !outputDir.endsWith("/") && !outputDir.endsWith("\\")){
						outputDir += File.separator;
					}
				}else if(args[i].equals("-gz")){
					outputGZIP = true;
				}else if(args[i].equals("-eA")){
					double d = Double.parseDouble(args[++i]);
					editMaxA = d < 1.0 ? -d : d;
				}else if(args[i].equals("--printprocessed")){
					printProcessedInterval = Long.parseLong(args[++i]);
				}else if(args[i].equals("-Q")){
					qualityFilter = Double.parseDouble(args[++i]);
				}else if(args[i].equals("-fA")){
					maxOffsetA = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-O")){
					replaceOriginal = true;
				}else if(args[i].equals("--min")){
					minLength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--max")){
					maxLength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-a")){
					while(i + 1 < args.length && !args[i + 1].startsWith("-")){
						if(args[i + 1].startsWith("^")){
							adaptersF.add(new Adapter(args[i + 1].substring(1), true, true));
						}else{
							adaptersF.add(new Adapter(args[i + 1], true, false));
						}
						i++;
					}
				}else if(args[i].equals("-A")){
					while(i + 1 < args.length && !args[i + 1].startsWith("-")){
						if(args[i + 1].endsWith("$")){
							adaptersF.add(new Adapter(args[i + 1].substring(0, args[i + 1].length() - 1), false, true));
						}else{
							adaptersF.add(new Adapter(args[i + 1], false, false));
						}
						i++;
					}
				}else if(args[i].equals("-oA")){
					minOverlapA = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--altqtrim")){
					trimAlgorithm = true;
					qualityTrimLength = Integer.MAX_VALUE;
				}else if(args[i].equals("-q")){
					if(i + 2 >= args.length || args[i + 2].startsWith("-")){
						qualityTrimQScore2 = Integer.parseInt(args[++i]);
					}else{
						qualityTrimQScore1 = Integer.parseInt(args[++i]);
						qualityTrimQScore2 = Integer.parseInt(args[++i]);
					}
				}else if(args[i].equals("--qtrimlen")){
					qualityTrimLength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-N")){
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						removeDNAWithNPercent = Double.parseDouble(args[++i]);
					}else{
						removeDNAWithNPercent = -1.0;
					}
				}else if(args[i].equals("-n")){
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						trimNPercent = Double.parseDouble(args[++i]);
					}else{
						trimNPercent = 0.5;
					}
				}else if(args[i].equals("-iA")){
					allowIndelsA = true;
				}else if(args[i].equals("--altqfilter")){
					filterAlgorithm = true;
				}else if(args[i].equals("-w")){
					wildcard = true;
				}else if(args[i].equals("-nQ")){
					removeUntrimmedReads = true;
				}else if(args[i].equals("-nA")){
					removeNoAdapterReads = true;
				}else if(args[i].equals("-p")){
					parallel = true;
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						splitBatchSize = Integer.parseInt(args[++i]);
					}
				}else if(args[i].equals("-pA")){
					probA = 0.5;
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						probA = Double.parseDouble(args[++i]);
					}
				}else{
					throw new Exception("Command not supported: " + args[i]);
				}
			}
			if(outputDir == null){
				outputDir = inputFile.getParent();
				if(!outputDir.endsWith(File.separator) && !outputDir.endsWith("/") && !outputDir.endsWith("\\")){
					outputDir += File.separator;
				}
			}else{
				File f = new File(outputDir);
				if(!f.exists())
					f.mkdirs();
			}
			if(inputFile.getAbsolutePath().toLowerCase().endsWith(".gz") || inputFile.getAbsolutePath().toLowerCase().endsWith(".gzip")){
				inputGZIP = true;
			}
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outputDir + "FastQParse_" + Mode.FILTER.description2 + ".log"), BUFFER_SIZE_LOG));
			Date date = new Date();
			logWriter.println("Started on: " + DATE_FORMAT.format(date));
			System.out.println("Started on: " + DATE_FORMAT.format(date));
			new FastQParseMain(Mode.FILTER);
			
			date = new Date();
			logWriter.println("Ended on: " + DATE_FORMAT.format(date));
			System.out.println("Ended on: " + DATE_FORMAT.format(date));
		}else if(args[0].equals("--simulate")){
			for(int i = 1; i < args.length; i++){
				if(args[i].equals("-o")){
					outputDir = args[++i];
					if(!outputDir.endsWith(File.separator) && !outputDir.endsWith("/") && !outputDir.endsWith("\\")){
						outputDir += File.separator;
					}
				}else if(args[i].equals("-s")){
					sampleFile = new File(args[++i]);
				}else if(args[i].equals("-eB")){
					editMaxB = Double.parseDouble(args[++i]);
				}else if(args[i].equals("-eA")){
					editMaxA = Double.parseDouble(args[++i]);
				}else if(args[i].equals("-eM")){
					editMaxM = Double.parseDouble(args[++i]);
				}else if(args[i].equals("-A")){
					while(i + 1 < args.length && !args[i + 1].startsWith("-")){
						if(args[i + 1].endsWith("$")){
							adaptersF.add(new Adapter(args[i + 1].substring(0, args[i + 1].length() - 1), false, true));
						}else{
							adaptersF.add(new Adapter(args[i + 1], false, false));
						}
						i++;
					}
				}else if(args[i].equals("-Z")){
					while(i + 1 < args.length && !args[i + 1].startsWith("-")){
						if(args[i + 1].endsWith("$")){
							adaptersR.add(new Adapter(args[i + 1].substring(0, args[i + 1].length() - 1), false, true));
						}else{
							adaptersR.add(new Adapter(args[i + 1], false, false));
						}
						i++;
					}
				}else if(args[i].equals("--umi")){
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						randUMILength = Integer.parseInt(args[++i]);
					}else{
						randUMILength = 12;
					}
					simUMI = true;
				}else if(args[i].equals("--reversed")){
					simReversed = true;
				}else if(args[i].equals("--merging")){
					simMerging = true;
					simReversed = true;
				}else if(args[i].equals("--iter")){
					simIter = Long.parseLong(args[++i]);
				}else if(args[i].equals("--len")){
					simReadLength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-iB")){
					allowIndelsB = true;
				}else if(args[i].equals("-iA")){
					allowIndelsA = true;
				}else if(args[i].equals("-oA")){
					minOverlapA = Integer.parseInt(args[++i]);
				}else{
					throw new Exception("Command not supported: " + args[i]);
				}
			}
			
			File f = new File(outputDir);
			if(!f.exists())
				f.mkdirs();
			
			File f2 = new File(outputDir + "results" + File.separatorChar);
			if(!f2.exists())
				f2.mkdirs();
			
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outputDir + "FastQParse_" + Mode.SIM_READS.description2 + ".log"), BUFFER_SIZE_LOG));
			
			new FastQParseMain(Mode.SIM_READS);
		}else{
			throw new Exception("Command not supported: " + args[0]);
		}
		if(logWriter != null){
			logWriter.close();
		}
	}
}
