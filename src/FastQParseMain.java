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
import java.util.Arrays;
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
	private static final int BUFFER_SIZE = 16384; //buffer size for buffered reader/writer
	private static final int BUFFER_SIZE_GZIP = 16384; //buffer size for gzip stream
	private static final int BUFFER_SIZE_LOG = 8192; //buffer size for log/stats files
	private static final String description2 = "+"; //the second description (3rd line in each read)
	
	private static PrintWriter logWriter; //prints to log file
	
	private static File inputFileF; //.fastq file with all data
	private static File inputFileR; //.fastq file with data for paired-ends sequencing
	private static File sampleInfoFile; //.txt file with sample barcode sequences
	private static File indexFileF; //index file for UMI
	private static File indexFileR; //index file for paired-ends sequencing
	private static String outputDir; //directory where output files go
	
	private static int randUMILength; //length of random UMI sequence (could be 0)
	
	private static boolean keepFirstDup = false; //remove duplicates that are not the first
	private static boolean keepBestDup = false; //remove duplicates that are not the best
	
	private static int maxOffsetB = 0; //tries offsets from 0 to maxOffset when matching barcode/enzyme
	private static int maxOffsetA = 0; //tries offsets from 0 to maxOffset when matching adapter
	
	private static double maxNPercent = 2.0; //if the percentage of DNA that is N is greater than this value then it will be undetermined
	
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
	
	private static int minOverlapA = 3; //minimum overlap when matching adapters
	private static int minOverlapB = Integer.MAX_VALUE; //minimum overlap when matching barcodes
	private static int minOverlapM = 3; //minimum overlap when matching for merging location
	
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
	private static boolean removeUnmergedReads = false;
	
	private static boolean parallel = false; //use parallel streams or not
	private static int splitBatchSize = 1024; //batch size when using ReadSpliterator
	
	private static boolean simReversed = false; //generate simulated data
	private static boolean simUMI = false;
	private static boolean simMerging = false;
	private static int simReadLength = 100;
	private static long simIter = 1000;
	
	private static double probB = -1.0; //prior probabilities for probability based matching
	private static double probA = -1.0;
	private static double probM = -1.0;
	
	private enum Mode{ //features
		PROCESS("Process Reads", "process"), SIM_READS("Generate Simulated Reads", "sim");
		
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
	
	//non-static variables that are modified when the program is running
	private HashMap<String, String> sampleMapF = new HashMap<String, String>(); //maps barcode sequence to sample name
	private ArrayList<String> sampleBarcodeF = new ArrayList<String>(); //names saved twice to keep them in order
	private ArrayList<String> sampleBarcodeR = new ArrayList<String>();
	
	private ArrayList<ArrayList<String>> sampleEnzymesF = new ArrayList<ArrayList<String>>();
	private ArrayList<ArrayList<String>> sampleEnzymesR = new ArrayList<ArrayList<String>>();
	
	private ArrayList<HashMap<Character, BitVector>> barcodePatternsF;
	private ArrayList<HashMap<Character, BitVector>> barcodePatternsR;
	private ArrayList<ArrayList<HashMap<Character, BitVector>>> enzymePatternsF;
	private ArrayList<ArrayList<HashMap<Character, BitVector>>> enzymePatternsR;
	private ArrayList<HashMap<Character, BitVector>> adapterPatternsF;
	private ArrayList<HashMap<Character, BitVector>> adapterPatternsR;
	
	private LongAdder totalReadsProcessed = new LongAdder(); //total DNA reads
	private long duplicatesRemoved = 0L; //total duplicates removed
	private LongAdder undeterminedReads = new LongAdder(); //total undetermined DNA
	private LongAdder totalRemovedAdapters = new LongAdder(); //total DNA with adapters removed
	private LongAdder totalQualityTrimmed = new LongAdder(); //total DNA that has been quality trimmed
	private LongAdder baseCount = new LongAdder(); //total base count
	private LongAdder totalReadsMerged = new LongAdder(); //total number of reads that are actually merged
	
	private boolean hasReversedBarcode = false;
	
	private int generatedDupFiles = 0; //total generated duplicate files
	private boolean generatedUndeterminedFile = false;
	
	private static String newEnzymes = ""; //enzymes added through a command, only used for printing to log file
	
	private long startTime;
	
	private HashMap<String, EnumMap<Stat, String>> stats = new HashMap<String, EnumMap<Stat, String>>(); //stats for each sample
	
	//constructor handles printing to log and calling the methods to process the files
	public FastQParseMain(Mode mode) throws Exception{
		logWriter.println("Mode: " + mode.description1);
		System.out.println("Mode: " + mode.description1);
		logWriter.println();
		
		if(mode == Mode.PROCESS){
			if(sampleInfoFile == null){
				logWriter.println("Demultiplex: false");
				logWriter.println();
			}else{
				logWriter.println("Demultiplex: true");
				logWriter.println();
				
				logWriter.println("Sample Info File: " + sampleInfoFile.getAbsolutePath());
				logWriter.println("Processing Sample File...");
				logWriter.println();
				logWriter.flush();
				
				readSample(); //read sample file
				
				logWriter.println("Sample File Processing Finished.");
				logWriter.println("Number of Samples: " + sampleBarcodeF.size());
				logWriter.println("Has Reversed Barcodes: " + hasReversedBarcode);
				logWriter.println();
			}
			
			logWriter.println("Addition Enzymes: " + newEnzymes);
			logWriter.println();
			
			logWriter.println("Forwards Read File: " + inputFileF.getAbsolutePath());
			if(indexFileF != null)
				logWriter.println("Forwards Index File: " + indexFileF.getAbsolutePath());
			if(inputFileR != null){
				logWriter.println("Reversed Read File: " + inputFileR.getAbsolutePath());
				if(indexFileR != null)
					logWriter.println("Reversed Index File: " + indexFileR.getAbsolutePath());
			}
			logWriter.println("Is Input GZIPPED: " + inputGZIP);
			logWriter.println();
			
			logWriter.println("Output Directory: " + outputDir);
			logWriter.println("Is Output GZIPPED: " + outputGZIP);
			logWriter.println();
			
			logWriter.println("Parallel: " + parallel);
			logWriter.println("Parallel Batch Size: " + splitBatchSize);
			logWriter.println();
			
			logWriter.println("Length of Random UMI: " + randUMILength);
			logWriter.println("Only Keep First Duplicate: " + keepFirstDup);
			logWriter.println("Only Keep Best Duplicate: " + keepBestDup);
			logWriter.println("Save Duplicates in Separate Files: " + saveDup);
			logWriter.println("Save Temporary Files: " + saveTemp);
			logWriter.println();
			
			if(editMaxB < 0.0)
				logWriter.println("Edit Percentage for Barcodes/Enzymes: " + -editMaxB);
			else
				logWriter.println("Max Edits for Barcodes/Enzymes: " + editMaxB);
			logWriter.println("Allow Insertions and Deletions for Barcodes/Enzymes: " + allowIndelsB);
			if(probB < 0.0)
				logWriter.println("Probability Based Matching for Barcodes: false");
			else
				logWriter.println("Prior Probability for Barcode Matching: " + probB);
			logWriter.println("Minimum Barcode Overlap: " + minOverlapB);
			logWriter.println("Maximum Offset For Barcodes/Enzymes: " + maxOffsetB);
			logWriter.println("Maximum Offset For Adapters: " + maxOffsetA);
			logWriter.println("Trim Barcode: " + removeBarRand);
			logWriter.println("Trim Enzyme: " + removeEnzyme);
			logWriter.println("Check Reversed Reads for Enzymes and Barcodes: " + checkReversedReads);
			logWriter.println("Remove Reads With Multiple Barcode Matches: " + singleBarcodeMatchOnly);
			logWriter.println();
			
			logWriter.println("Merge Paired-End Reads: " + mergePairedEnds);
			if(editMaxM < 0.0)
				logWriter.println("Edit Percentage for Merging Paired-Ends: " + -editMaxM);
			else
				logWriter.println("Max Edits for Merging Paired-Ends: " + editMaxM);
			if(probM < 0.0)
				logWriter.println("Probability Based Matching for Merging Paired End Reads: false");
			else
				logWriter.println("Prior Probability for Paired End Merging: " + probM);
			logWriter.println("Minimum Merge Location Overlap: " + minOverlapM);
			logWriter.println("Remove Reads That Are Not Merged: " + removeUnmergedReads);
			logWriter.println();
			
			String adaptersFString = adaptersF.toString();
			logWriter.println("Forwards Read Adapters: " + adaptersFString.substring(1, adaptersFString.length() - 1));
			if(inputFileR != null){
				String adaptersRString = adaptersR.toString();
				logWriter.println("Reversed Read Adapters: " + adaptersRString.substring(1, adaptersRString.length() - 1));
			}
			if(editMaxA < 0.0)
				logWriter.println("Edit Percentage for Adapters: " + -editMaxA);
			else
				logWriter.println("Max Edits for Adapters: " + editMaxA);
			logWriter.println("Allow Insertions and Deletions for Adapters: " + allowIndelsA);
			if(probA < 0.0)
				logWriter.println("Probability Based Matching for Adapters: false");
			else
				logWriter.println("Prior Probability for Adapter Matching: " + probA);
			logWriter.println("Minimum Adapter Overlap: " + minOverlapA);
			logWriter.println("Remove Reads That Do Not Contain Any Adapters: " + removeNoAdapterReads);
			logWriter.println();
			
			logWriter.println("Allow Wildcard Characters: " + wildcard);
			logWriter.println();
			
			logWriter.println("Quality Trim Algorithm: " + (trimAlgorithm ? "Sum" : "Local Average"));
			logWriter.println("5' Quality Trim Score Threshold: " + qualityTrimQScore1);
			logWriter.println("3' Quality Trim Score Threshold: " + qualityTrimQScore2);
			logWriter.println("Quality Trim Length: " + qualityTrimLength);
			logWriter.println("Remove Reads That Are Not Quality Trimmed: " + removeUntrimmedReads);
			logWriter.println();
			
			if(trimNPercent > 1.0)
				logWriter.println("Trim Leading and Trailing N: false");
			else
				logWriter.println("Percentage to Trim Leading and Trailing N By: " + trimNPercent);
			logWriter.println();
			
			if(maxNPercent >= 0 && maxNPercent <= 1.0)
				logWriter.println("Remove Reads With % of N Greater Than: " + maxNPercent);
			else if(maxNPercent < 0)
				logWriter.println("Remove Reads With at Least 1 N");
			else
				logWriter.println("Remove Reads With N: false");
			logWriter.println();
			
			logWriter.println("Quality Filter Threshold: " + qualityFilter);
			logWriter.println("Quality Filter Algorithm: " + (filterAlgorithm ? "Error Sum" : "Average"));
			logWriter.println();
			
			logWriter.println("Minimum Read Length: " + minLength);
			logWriter.println("Maximum Read Length: " + maxLength);
			logWriter.println();
			
			logWriter.println("Interval to Print Processed Reads: " + DECIMAL_FORMAT.format(printProcessedInterval));
			logWriter.println("Interval to Print Duplicate Reads: " + DECIMAL_FORMAT.format(printDuplicateInterval));
			logWriter.println();
			
			logWriter.println("Processing...");
			logWriter.println();
			logWriter.flush();
			
			//initialize the stats map
			EnumMap<Stat, String> map;
			if(sampleInfoFile == null){
				map = new EnumMap<Stat, String>(Stat.class);
				map.put(Stat.SEQUENCE_FORWARDS, "");
				if(hasReversedBarcode)
					map.put(Stat.SEQUENCE_REVERSED, "");
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
				stats.put("Kept", map);
			}else{
				for(int i = 0; i < sampleBarcodeF.size(); i++){
					map = new EnumMap<Stat, String>(Stat.class);
					map.put(Stat.SEQUENCE_FORWARDS, sampleBarcodeF.get(i));
					if(hasReversedBarcode)
						map.put(Stat.SEQUENCE_REVERSED, sampleBarcodeR.get(i));
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
					stats.put(sampleMapF.get(sampleBarcodeF.get(i)), map);
				}
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
			
			if(sampleInfoFile != null && probB < 0.0){ //generate pattern masks for edit based matching
				barcodePatternsF = new ArrayList<HashMap<Character, BitVector>>();
				for(int i = 0; i < sampleBarcodeF.size(); i++){
					barcodePatternsF.add(UtilMethods.genPatternMasks(sampleBarcodeF.get(i), allowIndelsB, wildcard));
				}
				
				enzymePatternsF = new ArrayList<ArrayList<HashMap<Character, BitVector>>>();
				for(int i = 0; i < sampleEnzymesF.size(); i++){
					enzymePatternsF.add(new ArrayList<HashMap<Character, BitVector>>());
					for(int j = 0; j < sampleEnzymesF.get(i).size(); j++){
						enzymePatternsF.get(i).add(UtilMethods.genPatternMasks(sampleEnzymesF.get(i).get(j), allowIndelsB, wildcard));
					}
				}
				
				if(inputFileR != null && checkReversedReads){
					barcodePatternsR = new ArrayList<HashMap<Character, BitVector>>();
					for(int i = 0; i < sampleBarcodeR.size(); i++){
						barcodePatternsR.add(UtilMethods.genPatternMasks(sampleBarcodeR.get(i), allowIndelsB, wildcard));
					}
					
					enzymePatternsR = new ArrayList<ArrayList<HashMap<Character, BitVector>>>();
					for(int i = 0; i < sampleEnzymesR.size(); i++){
						enzymePatternsR.add(new ArrayList<HashMap<Character, BitVector>>());
						for(int j = 0; j < sampleEnzymesR.get(i).size(); j++){
							enzymePatternsR.get(i).add(UtilMethods.genPatternMasks(sampleEnzymesR.get(i).get(j), allowIndelsB, wildcard));
						}
					}
				}
			}
			
			if(probA < 0.0){
				adapterPatternsF = new ArrayList<HashMap<Character, BitVector>>();
				for(int i = 0; i < adaptersF.size(); i++){
					adapterPatternsF.add(UtilMethods.genPatternMasks(adaptersF.get(i).isStart ? adaptersF.get(i).str : UtilMethods.reverse(adaptersF.get(i).str), allowIndelsA, wildcard));
				}
				
				if(inputFileR != null){
					adapterPatternsR = new ArrayList<HashMap<Character, BitVector>>();
					for(int i = 0; i < adaptersR.size(); i++){
						adapterPatternsR.add(UtilMethods.genPatternMasks(adaptersR.get(i).isStart ? adaptersR.get(i).str : UtilMethods.reverse(adaptersR.get(i).str), allowIndelsA, wildcard));
					}
				}
			}
			
			startTime = System.currentTimeMillis();
			int generatedSampleFiles = 0;
			
			//deduplicate after demultiplex
			if((keepFirstDup || keepBestDup) && indexFileF != null){ //assume that index files are provided for deduplicating for simpler logic
				ConcurrentLinkedQueue<Strings> files = processFiles();
				
				logWriter.println();
				logWriter.println("Processing Completed.");
				logWriter.println();
				logWriter.println("Run Time So Far: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
				logWriter.println();
				logWriter.println("Deduplicating...");
				logWriter.flush();
				
				for(Strings str : files){
					File file1 = new File(str.get(0));
					int pos = file1.getName().indexOf(".");
					String justNameF = pos > 0 ? file1.getName().substring(0, pos) : file1.getName();
					
					File file2 = new File(str.get(1));
					pos = file2.getName().indexOf(".");
					String justNameR = pos > 0 ? file2.getName().substring(0, pos) : file2.getName();
					
					logWriter.println();
					logWriter.println("Deduplicating File: " + justNameF);
					logWriter.flush();
					//deduplicate
					long removed = deduplicate(str.get(0), inputFileR == null ? null : str.get(2), str.get(1), inputFileR == null ? null : str.get(3),
							outputDir, outputDir + "dup" + File.separatorChar + justNameF + "_dup.fastq.gz", outputDir + "dup" + File.separatorChar + justNameR + "_dup.fastq.gz");
					for(int j = 0; j < sampleBarcodeF.size(); j++){ //calculate duduplicated reads count
						if(justNameF.contains(sampleMapF.get(sampleBarcodeF.get(j)))){
							map = stats.get(sampleMapF.get(sampleBarcodeF.get(j)));
							map.put(Stat.DEDUP_COUNT, DECIMAL_FORMAT.format(removed));
							map.put(Stat.DEDUP_PERCENT, DECIMAL_FORMAT.format((double)removed / Double.parseDouble(map.get(Stat.SEQUENCE_COUNT).replace(",", ""))));
							break;
						}
					}
					generatedSampleFiles += str.size();
					System.gc();
				}
				if(generatedUndeterminedFile)
					generatedSampleFiles += (inputFileR == null ? 1 : 2) + (indexFileF == null ? 0 : (inputFileR == null ? 1 : 2));
				logWriter.println();
				logWriter.println("Deduplicating Completed.");
				if(!saveTemp){
					File tempFile = new File(outputDir + "temp" + File.separatorChar);
					if(tempFile.exists())
						UtilMethods.deleteFolder(tempFile);
				}
			}else{ //no deduplicating after demultiplex
				ConcurrentLinkedQueue<Strings> files = processFiles();
				
				for(Strings str : files)
					generatedSampleFiles += str.size();
				if(generatedUndeterminedFile)
					generatedSampleFiles += (inputFileR == null ? 1 : 2) + (indexFileF == null ? 0 : (inputFileR == null ? 1 : 2));
				
				logWriter.println();
				logWriter.println("Processing Completed.");
				File tempFile = new File(outputDir + "temp" + File.separatorChar);
				if(tempFile.exists())
					UtilMethods.deleteFolder(tempFile);
			}
			
			logWriter.println();
			logWriter.println("Total Run Time: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
			logWriter.println();
			logWriter.println("Number of Lines Processed: " + DECIMAL_FORMAT.format(totalReadsProcessed.sum() * 4));
			logWriter.println("Number of Reads Processed: " + DECIMAL_FORMAT.format(totalReadsProcessed.sum()));
			logWriter.println("Number of Generated Sample Files: " + generatedSampleFiles);
			logWriter.println("Number of Generated Duplicate Files: " + generatedDupFiles);
			
			//fill in info for total stats
			map = new EnumMap<Stat, String>(Stat.class);
			map.put(Stat.SEQUENCE_FORWARDS, "");
			if(hasReversedBarcode)
				map.put(Stat.SEQUENCE_REVERSED, "");
			map.put(Stat.SEQUENCE_COUNT, DECIMAL_FORMAT.format(totalReadsProcessed.sum()));
			map.put(Stat.SEQUENCE_PERCENT, DECIMAL_FORMAT.format(1));
			map.put(Stat.MERGED_COUNT, DECIMAL_FORMAT.format(totalReadsMerged.sum()));
			map.put(Stat.MERGED_PERCENT, DECIMAL_FORMAT.format((double)totalReadsMerged.sum() / (double)totalReadsProcessed.sum()));
			map.put(Stat.DEDUP_COUNT, DECIMAL_FORMAT.format(duplicatesRemoved));
			map.put(Stat.DEDUP_PERCENT, DECIMAL_FORMAT.format((double)duplicatesRemoved / (double)totalReadsProcessed.sum()));
			map.put(Stat.REMOVEADAPTER_COUNT, DECIMAL_FORMAT.format(totalRemovedAdapters.sum()));
			map.put(Stat.REMOVEADAPTER_PERCENT, DECIMAL_FORMAT.format((double)totalRemovedAdapters.sum() / (double)totalReadsProcessed.sum()));
			map.put(Stat.QUALITYTRIM_COUNT, DECIMAL_FORMAT.format(totalQualityTrimmed.sum()));
			map.put(Stat.QUALITYTRIM_PERCENT, DECIMAL_FORMAT.format((double)totalQualityTrimmed.sum() / (double)totalReadsProcessed.sum()));
			map.put(Stat.BASE_COUNT, DECIMAL_FORMAT.format(baseCount.sum()));
			map.put(Stat.BASE_PERCENT, DECIMAL_FORMAT.format(1));
			stats.put("Total", map);
			
			//print stats file with cool padding
			int width = 25;
			int space = 1;
			
			PrintWriter statsWriter = new PrintWriter(new BufferedWriter(new FileWriter(outputDir + "FastQParse_" + mode.description2 + ".stats"), BUFFER_SIZE_LOG));
			StringBuilder builder = new StringBuilder();
			String colTitle1 = sampleInfoFile == null ? "File" : "Sample";
			builder.append(colTitle1);
			for(int i = 0; i < width - colTitle1.length(); i++){
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
			if(sampleInfoFile == null){
				String[] rows = {"Kept", "Undetermined", "Total"};
				for(String row : rows){
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
			}else{
				for(int i = 0; i < sampleBarcodeF.size() + 2; i++){
					String row = i == sampleBarcodeF.size() ? "Undetermined" : (i == sampleBarcodeF.size() + 1 ? "Total" : sampleMapF.get(sampleBarcodeF.get(i)));
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
			}
			statsWriter.close();
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
		BufferedReader reader = new BufferedReader(new FileReader(sampleInfoFile), BUFFER_SIZE_LOG);
		
		String tempLine;
		while((tempLine = reader.readLine()) != null){
			String[] line = tempLine.split("\\s+");
			sampleMapF.put(line[1].toUpperCase(), line[0]);
			sampleBarcodeF.add(line[1].toUpperCase());
			if(line.length > 4){ //if there are reversed barcodes
				sampleBarcodeR.add(line[4].toUpperCase());
				hasReversedBarcode = true;
			}
			sampleEnzymesF.add(EnzymeList.enzymes.get(line[2].toUpperCase()));
			if(line.length > 3){ //if there are reversed enzymes
				sampleEnzymesR.add(EnzymeList.enzymes.get(line[3].toUpperCase()));
			}else{ //reversed enzymes will be forwards enzymes
				sampleEnzymesR.add(EnzymeList.enzymes.get(line[2].toUpperCase()));
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
			
			if(r.nextFloat() >= 0.5f){ //have merge location within edit range
				simMerged += 2;
				mergeF = UtilMethods.randSeq(r, minOverlapM + r.nextInt(simReadLength - minOverlapM + 1));
				mergeR = UtilMethods.randEdit(r, UtilMethods.reverseComplement(mergeF), r.nextInt((int)editMaxM + 1), false);
			}else{
				if(r.nextFloat() >= 0.5f){ //merge length too short
					mergeF = UtilMethods.randSeq(r, r.nextInt(minOverlapM - 1) + 1);
					mergeR = UtilMethods.randEdit(r, UtilMethods.reverseComplement(mergeF), r.nextInt(Math.min((int)editMaxM + 1, mergeF.length())), false);
				}else{ //have merge location outside edit range
					mergeF = UtilMethods.randSeq(r, minOverlapM + r.nextInt(simReadLength - minOverlapM + 1));
					mergeR = UtilMethods.randEdit(r, UtilMethods.reverseComplement(mergeF), (int)editMaxM + r.nextInt(Math.min((int)editMaxM * 2 + 1, mergeF.length()) - (int)editMaxM), false);
				}
			}
			
			String seqF = UtilMethods.randSeq(r, simReadLength - mergeF.length());
			String readF = seqF + mergeF;
			String qualF = UtilMethods.randQuality(r, simReadLength);
			
			String seqR = UtilMethods.randSeq(r, simReadLength - mergeR.length());
			String readR = seqR + mergeR;
			String qualR = UtilMethods.randQuality(r, simReadLength);
			
			String mergedReadF = seqF + mergeF + UtilMethods.reverseComplement(seqR);
			String mergedQualF = UtilMethods.makeStr('A', mergedReadF.length()); //the real quality values don't really matter
			
			inWriterF.write("@SIMULATED READ");
			inWriterF.newLine();
			inWriterF.write(readF);
			inWriterF.newLine();
			inWriterF.write(description2);
			inWriterF.newLine();
			inWriterF.write(qualF);
			inWriterF.newLine();
			
			inWriterR.write("@SIMULATED READ");
			inWriterR.newLine();
			inWriterR.write(readR);
			inWriterR.newLine();
			inWriterR.write(description2);
			inWriterR.newLine();
			inWriterR.write(qualR);
			inWriterR.newLine();
			
			outWriterF.write("@SIMULATED READ");
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
			outWriterF[i] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_sample_" + sampleMapF.get(sampleBarcodeF.get(i)) + "_R1.fastq.gz"))));
		}
		
		if(simReversed){
			inWriterR = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "simulated_reads_R2.fastq.gz"))));
			outWriterR = new BufferedWriter[sampleMapF.size()];
			undeterminedWriterR = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_sample_undetermined_R2.fastq.gz"))));
			
			for(int i = 0; i < sampleMapF.size(); i++){
				outWriterR[i] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_sample_" + sampleMapF.get(sampleBarcodeF.get(i)) + "_R2.fastq.gz"))));
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
			boolean correctR = !checkReversedReads || r.nextFloat() >= 0.1f;
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
				readF = UtilMethods.randEdit(r, sampleBarcodeF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB) +
						UtilMethods.randEdit(r, sampleEnzymesF.get(barcode).get(r.nextInt(sampleEnzymesF.get(barcode).size())), r.nextInt((int)editMaxB + 1), allowIndelsB) + seqF;
				
				if(simReversed){
					seqR = UtilMethods.randSeq(r, simReadLength);
					readR = UtilMethods.randEdit(r, hasReversedBarcode ? sampleBarcodeR.get(barcode) : sampleBarcodeF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB) +
							UtilMethods.randEdit(r, sampleEnzymesR.get(barcode).get(r.nextInt(sampleEnzymesR.get(barcode).size())), r.nextInt((int)editMaxB + 1), allowIndelsB) + seqR;
				}
			}else if(correctF){ //forwards is correct, reversed is not
				simUndetermined += 2;
				undetermined = true;
				
				seqF = UtilMethods.randSeq(r, simReadLength);
				readF = UtilMethods.randEdit(r, sampleBarcodeF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB) +
						UtilMethods.randEdit(r, sampleEnzymesF.get(barcode).get(r.nextInt(sampleEnzymesF.get(barcode).size())), r.nextInt((int)editMaxB + 1), allowIndelsB) + seqF;
				
				seqR = UtilMethods.randSeq(r, simReadLength);
				String randBarcode;
				String randEnzyme;
				if(r.nextFloat() < 0.5f){ //error in barcode
					randBarcode = UtilMethods.randEdit(r, hasReversedBarcode ? sampleBarcodeR.get(barcode) : sampleBarcodeF.get(barcode),
							(int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, (hasReversedBarcode ? sampleBarcodeR.get(barcode) : sampleBarcodeF.get(barcode)).length()) - (int)editMaxB), allowIndelsB);
					randEnzyme = UtilMethods.randEdit(r, sampleEnzymesR.get(barcode).get(r.nextInt(sampleEnzymesR.get(barcode).size())), r.nextInt((int)editMaxB + 1), allowIndelsB);
				}else{ //error in enzyme
					randBarcode = UtilMethods.randEdit(r, hasReversedBarcode ? sampleBarcodeR.get(barcode) : sampleBarcodeF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB);
					int enzyme = r.nextInt(sampleEnzymesR.get(barcode).size());
					randEnzyme = UtilMethods.randEdit(r, sampleEnzymesR.get(barcode).get(enzyme), (int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, sampleEnzymesR.get(barcode).get(enzyme).length()) - (int)editMaxB), allowIndelsB);
				}
				readR = randBarcode + randEnzyme + seqR;
			}else if(simReversed && correctR){ //reversed is correct, forwards is not
				simUndetermined += 2;
				undetermined = true;
				
				seqF = UtilMethods.randSeq(r, simReadLength);
				String randBarcode;
				String randEnzyme;
				if(r.nextFloat() < 0.5f){ //error in barcode
					randBarcode = UtilMethods.randEdit(r, sampleBarcodeF.get(barcode), (int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, sampleBarcodeF.get(barcode).length()) - (int)editMaxB), allowIndelsB);
					randEnzyme = UtilMethods.randEdit(r, sampleEnzymesF.get(barcode).get(r.nextInt(sampleEnzymesF.get(barcode).size())), r.nextInt((int)editMaxB + 1), allowIndelsB);
				}else{ //error in enzyme
					randBarcode = UtilMethods.randEdit(r, sampleBarcodeF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB);
					int enzyme = r.nextInt(sampleEnzymesF.get(barcode).size());
					randEnzyme = UtilMethods.randEdit(r, sampleEnzymesF.get(barcode).get(enzyme), (int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, sampleEnzymesF.get(barcode).get(enzyme).length()) - (int)editMaxB), allowIndelsB);
				}
				readF = randBarcode + randEnzyme + seqF;
				
				seqR = UtilMethods.randSeq(r, simReadLength);
				readR = UtilMethods.randEdit(r, hasReversedBarcode ? sampleBarcodeR.get(barcode) : sampleBarcodeF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB) +
						UtilMethods.randEdit(r, sampleEnzymesR.get(barcode).get(r.nextInt(sampleEnzymesR.get(barcode).size())), r.nextInt((int)editMaxB + 1), allowIndelsB) + seqR;
			}else{ //forwards is not correct, reversed could be correct
				simUndetermined += simReversed ? 2 : 1;
				undetermined = true;
				
				seqF = UtilMethods.randSeq(r, simReadLength);
				String randBarcode;
				String randEnzyme;
				if(r.nextFloat() < 0.5f){ //error in barcode
					randBarcode = UtilMethods.randEdit(r, sampleBarcodeF.get(barcode), (int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, sampleBarcodeF.get(barcode).length()) - (int)editMaxB), allowIndelsB);
					randEnzyme = UtilMethods.randEdit(r, sampleEnzymesF.get(barcode).get(r.nextInt(sampleEnzymesF.get(barcode).size())), r.nextInt((int)editMaxB + 1), allowIndelsB);
				}else{ //error in enzyme
					randBarcode = UtilMethods.randEdit(r, sampleBarcodeF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB);
					int enzyme = r.nextInt(sampleEnzymesF.get(barcode).size());
					randEnzyme = UtilMethods.randEdit(r, sampleEnzymesF.get(barcode).get(enzyme), (int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, sampleEnzymesF.get(barcode).get(enzyme).length()) - (int)editMaxB), allowIndelsB);
				}
				readF = randBarcode + randEnzyme + seqF;
				
				if(simReversed){
					seqR = UtilMethods.randSeq(r, simReadLength);
					if(r.nextFloat() < 0.5f){ //error in barcode
						randBarcode = UtilMethods.randEdit(r, hasReversedBarcode ? sampleBarcodeR.get(barcode) : sampleBarcodeF.get(barcode),
								(int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, (hasReversedBarcode ? sampleBarcodeR.get(barcode) : sampleBarcodeF.get(barcode)).length()) - (int)editMaxB), allowIndelsB);
						randEnzyme = UtilMethods.randEdit(r, sampleEnzymesR.get(barcode).get(r.nextInt(sampleEnzymesR.get(barcode).size())), r.nextInt((int)editMaxB + 1), allowIndelsB);
					}else{ //error in enzyme
						randBarcode = UtilMethods.randEdit(r, hasReversedBarcode ? sampleBarcodeR.get(barcode) : sampleBarcodeF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB);
						int enzyme = r.nextInt(sampleEnzymesR.get(barcode).size());
						randEnzyme = UtilMethods.randEdit(r, sampleEnzymesR.get(barcode).get(enzyme), (int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, sampleEnzymesR.get(barcode).get(enzyme).length()) - (int)editMaxB), allowIndelsB);
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
			
			inWriterF.write("@SIMULATED READ");
			inWriterF.newLine();
			inWriterF.write(readF);
			inWriterF.newLine();
			inWriterF.write(description2);
			inWriterF.newLine();
			inWriterF.write(UtilMethods.makeStr('A', readF.length()));
			inWriterF.newLine();
			
			if(undetermined){
				undeterminedWriterF.write("@SIMULATED READ");
				undeterminedWriterF.newLine();
				undeterminedWriterF.write(readF);
				undeterminedWriterF.newLine();
				undeterminedWriterF.write(description2);
				undeterminedWriterF.newLine();
				undeterminedWriterF.write(UtilMethods.makeStr('A', readF.length()));
				undeterminedWriterF.newLine();
			}else{
				outWriterF[barcode].write("@SIMULATED READ");
				outWriterF[barcode].newLine();
				outWriterF[barcode].write(seqF);
				outWriterF[barcode].newLine();
				outWriterF[barcode].write(description2);
				outWriterF[barcode].newLine();
				outWriterF[barcode].write(UtilMethods.makeStr('A', seqF.length()));
				outWriterF[barcode].newLine();
			}
			
			if(simReversed){
				inWriterR.write("@SIMULATED READ");
				inWriterR.newLine();
				inWriterR.write(readR);
				inWriterR.newLine();
				inWriterR.write(description2);
				inWriterR.newLine();
				inWriterR.write(UtilMethods.makeStr('A', readR.length()));
				inWriterR.newLine();
				
				if(undetermined){
					undeterminedWriterR.write("@SIMULATED READ");
					undeterminedWriterR.newLine();
					undeterminedWriterR.write(readR);
					undeterminedWriterR.newLine();
					undeterminedWriterR.write(description2);
					undeterminedWriterR.newLine();
					undeterminedWriterR.write(UtilMethods.makeStr('A', readR.length()));
					undeterminedWriterR.newLine();
				}else{
					outWriterR[barcode].write("@SIMULATED READ");
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
			logWriter.println(sampleMapF.get(sampleBarcodeF.get(i)) + "\t" + DECIMAL_FORMAT.format(simCounts[i][0]) + "\t" + DECIMAL_FORMAT.format(simCounts[i][1]));
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
			outWriterF[i] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_sample_" + sampleMapF.get(sampleBarcodeF.get(i)) + "_R1.fastq.gz"))));
			outWriterF2[i] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_index_" + sampleMapF.get(sampleBarcodeF.get(i)) + "_R1.fastq.gz"))));
		}
		
		if(simReversed){
			inWriterR = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "simulated_reads_R2.fastq.gz"))));
			inWriterR2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "simulated_indexes_R2.fastq.gz"))));
			outWriterR = new BufferedWriter[sampleMapF.size()];
			outWriterR2 = new BufferedWriter[sampleMapF.size()];
			undeterminedWriterR = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_sample_undetermined_R2.fastq.gz"))));
			undeterminedWriterR2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_index_undetermined_R2.fastq.gz"))));
			
			for(int i = 0; i < sampleMapF.size(); i++){
				outWriterR[i] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_sample_" + sampleMapF.get(sampleBarcodeF.get(i)) + "_R2.fastq.gz"))));
				outWriterR2[i] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "results" + File.separatorChar + "sim_index_" + sampleMapF.get(sampleBarcodeF.get(i)) + "_R2.fastq.gz"))));
			}
		}
		
		long[][] simCounts = new long[sampleMapF.size()][1];
		long simUndetermined = 0;
		
		for(long i = 0; i < simIter; i++){
			Random r = new Random();
			
			int barcode = r.nextInt(sampleMapF.size());
			boolean correctF = r.nextFloat() >= 0.1f;
			boolean correctR = !checkReversedReads || r.nextFloat() >= 0.1f;
			boolean undetermined = false;
			
			String readF = null;
			String indexF = null;
			String readR = null;
			String indexR = null;
			
			if(correctF && (!simReversed || correctR)){
				simCounts[barcode][0] += simReversed ? 2 : 1;
				
				readF = UtilMethods.randSeq(r, simReadLength);
				indexF = UtilMethods.randEdit(r, sampleBarcodeF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB) + UtilMethods.randSeq(r, randUMILength);
				
				if(simReversed){
					readR = UtilMethods.randSeq(r, simReadLength);
					indexR = UtilMethods.randEdit(r, hasReversedBarcode ? sampleBarcodeR.get(barcode) : sampleBarcodeF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB) + UtilMethods.randSeq(r, randUMILength);
				}
			}else if(correctF){
				simUndetermined += 2;
				undetermined = true;
				
				readF = UtilMethods.randSeq(r, simReadLength);
				indexF = UtilMethods.randEdit(r, sampleBarcodeF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB) + UtilMethods.randSeq(r, randUMILength);
				
				readR = UtilMethods.randSeq(r, simReadLength);
				indexR = UtilMethods.randEdit(r, hasReversedBarcode ? sampleBarcodeR.get(barcode) : sampleBarcodeF.get(barcode),
						(int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, (hasReversedBarcode ? sampleBarcodeR.get(barcode) : sampleBarcodeF.get(barcode)).length()) - (int)editMaxB), allowIndelsB) + UtilMethods.randSeq(r, randUMILength);
			}else if(simReversed && correctR){
				simUndetermined += 2;
				undetermined = true;
				
				readF = UtilMethods.randSeq(r, simReadLength);
				indexF = UtilMethods.randEdit(r, sampleBarcodeF.get(barcode), (int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, sampleBarcodeF.get(barcode).length()) - (int)editMaxB), allowIndelsB) + UtilMethods.randSeq(r, randUMILength);
				
				readR = UtilMethods.randSeq(r, simReadLength);
				indexR = UtilMethods.randEdit(r, hasReversedBarcode ? sampleBarcodeR.get(barcode) : sampleBarcodeF.get(barcode), r.nextInt((int)editMaxB + 1), allowIndelsB) + UtilMethods.randSeq(r, randUMILength);
			}else{
				simUndetermined += simReversed ? 2 : 1;
				undetermined = true;
				
				readF = UtilMethods.randSeq(r, simReadLength);
				indexF = UtilMethods.randEdit(r, sampleBarcodeF.get(barcode), (int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, sampleBarcodeF.get(barcode).length()) - (int)editMaxB), allowIndelsB) + UtilMethods.randSeq(r, randUMILength);
				
				readR = UtilMethods.randSeq(r, simReadLength);
				indexR = UtilMethods.randEdit(r, hasReversedBarcode ? sampleBarcodeR.get(barcode) : sampleBarcodeF.get(barcode),
						(int)editMaxB + r.nextInt(Math.min((int)editMaxB * 2 + 1, (hasReversedBarcode ? sampleBarcodeR.get(barcode) : sampleBarcodeF.get(barcode)).length()) - (int)editMaxB), allowIndelsB) + UtilMethods.randSeq(r, randUMILength);
			}
			
			inWriterF.write("@SIMULATED READ");
			inWriterF.newLine();
			inWriterF.write(readF);
			inWriterF.newLine();
			inWriterF.write(description2);
			inWriterF.newLine();
			inWriterF.write(UtilMethods.makeStr('A', readF.length()));
			inWriterF.newLine();
			
			inWriterF2.write("@SIMULATED INDEX");
			inWriterF2.newLine();
			inWriterF2.write(indexF);
			inWriterF2.newLine();
			inWriterF2.write(description2);
			inWriterF2.newLine();
			inWriterF2.write(UtilMethods.makeStr('A', indexF.length()));
			inWriterF2.newLine();
			
			if(undetermined){
				undeterminedWriterF.write("@SIMULATED READ");
				undeterminedWriterF.newLine();
				undeterminedWriterF.write(readF);
				undeterminedWriterF.newLine();
				undeterminedWriterF.write(description2);
				undeterminedWriterF.newLine();
				undeterminedWriterF.write(UtilMethods.makeStr('A', readF.length()));
				undeterminedWriterF.newLine();
				
				undeterminedWriterF2.write("@SIMULATED INDEX");
				undeterminedWriterF2.newLine();
				undeterminedWriterF2.write(indexF);
				undeterminedWriterF2.newLine();
				undeterminedWriterF2.write(description2);
				undeterminedWriterF2.newLine();
				undeterminedWriterF2.write(UtilMethods.makeStr('A', indexF.length()));
				undeterminedWriterF2.newLine();
			}else{
				outWriterF[barcode].write("@SIMULATED READ");
				outWriterF[barcode].newLine();
				outWriterF[barcode].write(readF);
				outWriterF[barcode].newLine();
				outWriterF[barcode].write(description2);
				outWriterF[barcode].newLine();
				outWriterF[barcode].write(UtilMethods.makeStr('A', readF.length()));
				outWriterF[barcode].newLine();
				
				outWriterF2[barcode].write("@SIMULATED INDEX");
				outWriterF2[barcode].newLine();
				outWriterF2[barcode].write(indexF);
				outWriterF2[barcode].newLine();
				outWriterF2[barcode].write(description2);
				outWriterF2[barcode].newLine();
				outWriterF2[barcode].write(UtilMethods.makeStr('A', indexF.length()));
				outWriterF2[barcode].newLine();
			}
			
			if(simReversed){
				inWriterR.write("@SIMULATED READ");
				inWriterR.newLine();
				inWriterR.write(readR);
				inWriterR.newLine();
				inWriterR.write(description2);
				inWriterR.newLine();
				inWriterR.write(UtilMethods.makeStr('A', readR.length()));
				inWriterR.newLine();
				
				inWriterR2.write("@SIMULATED INDEX");
				inWriterR2.newLine();
				inWriterR2.write(indexR);
				inWriterR2.newLine();
				inWriterR2.write(description2);
				inWriterR2.newLine();
				inWriterR2.write(UtilMethods.makeStr('A', indexR.length()));
				inWriterR2.newLine();
				
				if(undetermined){
					undeterminedWriterR.write("@SIMULATED READ");
					undeterminedWriterR.newLine();
					undeterminedWriterR.write(readR);
					undeterminedWriterR.newLine();
					undeterminedWriterR.write(description2);
					undeterminedWriterR.newLine();
					undeterminedWriterR.write(UtilMethods.makeStr('A', readR.length()));
					undeterminedWriterR.newLine();
					
					undeterminedWriterR2.write("@SIMULATED INDEX");
					undeterminedWriterR2.newLine();
					undeterminedWriterR2.write(indexR);
					undeterminedWriterR2.newLine();
					undeterminedWriterR2.write(description2);
					undeterminedWriterR2.newLine();
					undeterminedWriterR2.write(UtilMethods.makeStr('A', indexR.length()));
					undeterminedWriterR2.newLine();
				}else{
					outWriterR[barcode].write("@SIMULATED READ");
					outWriterR[barcode].newLine();
					outWriterR[barcode].write(readR);
					outWriterR[barcode].newLine();
					outWriterR[barcode].write(description2);
					outWriterR[barcode].newLine();
					outWriterR[barcode].write(UtilMethods.makeStr('A', readR.length()));
					outWriterR[barcode].newLine();
					
					outWriterR2[barcode].write("@SIMULATED INDEX");
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
			logWriter.println(sampleMapF.get(sampleBarcodeF.get(i)) + "\t" + DECIMAL_FORMAT.format(simCounts[i][0]));
		}
		logWriter.println("Undetermined\t" + DECIMAL_FORMAT.format(simUndetermined));
		logWriter.println("Total\t" + DECIMAL_FORMAT.format(simReversed ? simIter * 2 : simIter));
	}
	
	//the ConcurrentLinkedQueue contains strings that represent the output files
	private ConcurrentLinkedQueue<Strings> processFiles() throws Exception{
		BufferedReader readerF; //forwards input
		BufferedReader readerR = null; //reversed input
		BufferedReader readerIF = null; //forwards index input
		BufferedReader readerIR = null; //reversed index input
		if(inputGZIP){
			readerF = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFileF), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
			if(indexFileF != null)
				readerIF = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(indexFileF), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
			if(inputFileR != null){
				readerR = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFileR), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
				if(indexFileR != null)
					readerIR = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(indexFileR), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
			}
		}else{
			readerF = new BufferedReader(new FileReader(inputFileF), BUFFER_SIZE);
			if(indexFileF != null)
				readerIF = new BufferedReader(new FileReader(indexFileF), BUFFER_SIZE);
			if(inputFileR != null){
				readerR = new BufferedReader(new FileReader(inputFileR), BUFFER_SIZE);
				if(indexFileR != null)
					readerIR = new BufferedReader(new FileReader(indexFileR), BUFFER_SIZE);
			}
		}
		
		//if not demultiplexing, then the array sizes should be 1
		//all writers need to be initialized right here or else Java will complain about them changing in the lambda statement
		BufferedWriter[] writersF = new BufferedWriter[sampleInfoFile == null ? 1 : sampleMapF.size()]; //forwards output
		BufferedWriter[] writersR = !mergePairedEnds && inputFileR != null ? new BufferedWriter[sampleInfoFile == null ? 1 : sampleMapF.size()] : null; //reversed output
		BufferedWriter[] writersIF = indexFileF != null ? new BufferedWriter[sampleInfoFile == null ? 1 : sampleMapF.size()] : null; //forwards index output
		BufferedWriter[] writersIR = !mergePairedEnds && inputFileR != null && indexFileR != null ? new BufferedWriter[sampleInfoFile == null ? 1 : sampleMapF.size()] : null; //reversed index output
		
		BufferedWriter[] undeterminedWriterF = new BufferedWriter[2]; //for forwards and forwards index undetermined reads
		BufferedWriter[] undeterminedWriterR = inputFileR != null ? new BufferedWriter[2] : null; //for reversed and reversed index undetermined reads
		
		ConcurrentLinkedQueue<Strings> files = new ConcurrentLinkedQueue<Strings>();
		//LongAdders to keep track of stats
		LongAdder[][] readCounts = new LongAdder[(sampleInfoFile == null ? 1 : sampleMapF.size()) + 1][5];
		//0 = read count
		//1 = base pair count
		//2 = number of merged reads
		//3 = number of reads with removed adapters
		//4 = number of quality trimmed reads
		//note that forwards and reversed reads will contribute to the counts
		for(int i = 0; i < readCounts.length; i++){
			for(int j = 0; j < readCounts[i].length; j++){
				readCounts[i][j] = new LongAdder();
			}
		}
		//Objects are used as locks
		//each sample file gets its own lock so multiple threads can write to different samples at the same time
		Object[] locks = new Object[(sampleInfoFile == null ? 1 : sampleMapF.size()) + 1];
		for(int i = 0; i < locks.length; i++){
			locks[i] = new Object();
		}
		
		StreamSupport.stream(new ReadSpliterator<Read>(splitBatchSize, readerF, readerR, readerIF, readerIR), parallel).forEach((read) -> {
			totalReadsProcessed.add(inputFileR == null ? 1 : 2);
			
			boolean valid = false; //undetermined or not
			
			//if this is supported in parallel then the logWriter needs to be synchronized
			//easier and faster to just not print
			if(!parallel && totalReadsProcessed.sum() % printProcessedInterval == 0){
				logWriter.println("Reads Processed So Far: " + DECIMAL_FORMAT.format(totalReadsProcessed));
				logWriter.println("Run Time So Far: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
				logWriter.println("Current Memory Usage (GB): " + DECIMAL_FORMAT.format(UtilMethods.currentMemoryUsage()));
				logWriter.flush();
			}
			
			//save the quality and sequence lines so the undetermined reads are not modified
			String tempSequenceF = read.readF[1];
			String tempQualityF = read.readF[3];
			String tempSequenceR = null;
			String tempQualityR = null;
			if(inputFileR != null){
				tempSequenceR = read.readR[1];
				tempQualityR = read.readR[3];
			}
			
			//check the percentage of N
			if(maxNPercent > 1.0 || (maxNPercent >= 0.0 && UtilMethods.percentN(read.readF[1]) <= maxNPercent && (inputFileR == null || UtilMethods.percentN(read.readR[1]) <= maxNPercent)) ||
					(maxNPercent < 0.0 && UtilMethods.countN(read.readF[1]) <= 0 && (inputFileR == null || UtilMethods.countN(read.readR[1]) <= 0))){
				int barcodeIndex = -1;
				int barcodeEnd = -1;
				int barcodeLength = 0;
				int enzymeEnd = -1;
				int barcodeEnd2 = -1;
				int enzymeEnd2 = -1;
				int minEdit = Integer.MAX_VALUE;
				int barcodeMatchCount = 0; //used for checking if a read matches multiple barcodes
				
				if(sampleInfoFile != null){ //only if there are barcodes
					if(indexFileF == null){ //check for enzyme and barcode in reads, in the following order: forwards barcode, forwards enzyme, reversed barcode, reversed enzyme
						for(int i = 0; i < sampleBarcodeF.size(); i++){ //a barcode match is a barcode match if there is an enzyme after it
							ArrayList<Match> matches = null;
							if(probB < 0.0){
								matches = UtilMethods.searchWithN(read.readF[1], 0, Math.min(read.readF[1].length(), maxOffsetB + sampleBarcodeF.get(i).length() +
										(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * sampleBarcodeF.get(i).length()) : (int)editMaxB) : 0)), sampleBarcodeF.get(i), editMaxB, allowIndelsB, false, minOverlapB, wildcard, barcodePatternsF.get(i));
							}else{
								matches = UtilMethods.searchWithProb(read.readF[1], 0, Math.min(read.readF[1].length(), maxOffsetB + sampleBarcodeF.get(i).length() +
										(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * sampleBarcodeF.get(i).length()) : (int)editMaxB) : 0)), read.readF[3], sampleBarcodeF.get(i), null, probB, minOverlapB, wildcard);
							}
							boolean isMatch = false;
							for(int j = 0; j < matches.size(); j++){
								for(int k = 0; k < sampleEnzymesF.get(i).size(); k++){
									ArrayList<Match> enzymeMatches = null;
									if(probB < 0.0){
										enzymeMatches = UtilMethods.searchWithN(read.readF[1], matches.get(j).end + 1/* + randUMILength*/, Math.min(read.readF[1].length(), maxOffsetB + matches.get(j).end + 1 + /*randUMILength + */sampleEnzymesF.get(i).get(k).length() +
												(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * sampleEnzymesF.get(i).get(k).length()) : (int)editMaxB) : 0)), sampleEnzymesF.get(i).get(k), editMaxB, allowIndelsB, true, Integer.MAX_VALUE, wildcard, enzymePatternsF.get(i).get(k));
									}else{
										enzymeMatches = UtilMethods.searchWithProb(read.readF[1], matches.get(j).end + 1/* + randUMILength*/, Math.min(read.readF[1].length(), maxOffsetB + matches.get(j).end + 1 + /*randUMILength + */sampleEnzymesF.get(i).get(k).length() +
												(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * sampleEnzymesF.get(i).get(k).length()) : (int)editMaxB) : 0)), read.readF[3], sampleEnzymesF.get(i).get(k), null, probB, minOverlapB, wildcard);
									}
									if(!enzymeMatches.isEmpty()){
										int tempBarcodeEnd2 = -1;
										int tempBarcodeLength2 = 0;
										int tempEnzymeEnd2 = -1;
										if(inputFileR != null && checkReversedReads){
											int minEdit2 = Integer.MAX_VALUE;
											ArrayList<Match> rMatches = null; //check for the reversed barcode and enzyme
											if(probB < 0.0){
												rMatches = UtilMethods.searchWithN(read.readR[1], 0, Math.min(read.readR[1].length(), maxOffsetB + (hasReversedBarcode ? sampleBarcodeR.get(i).length() : sampleBarcodeF.get(i).length()) +
														(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * (hasReversedBarcode ? sampleBarcodeR.get(i).length() : sampleBarcodeF.get(i).length())) : (int)editMaxB) : 0)),
														hasReversedBarcode ? sampleBarcodeR.get(i) : UtilMethods.complement(sampleBarcodeF.get(i)), editMaxB, allowIndelsB, false, minOverlapB, wildcard, barcodePatternsR.get(i));
											}else{
												rMatches = UtilMethods.searchWithProb(read.readR[1], 0, Math.min(read.readR[1].length(), maxOffsetB + (hasReversedBarcode ? sampleBarcodeR.get(i).length() : sampleBarcodeF.get(i).length()) +
														(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * (hasReversedBarcode ? sampleBarcodeR.get(i).length() : sampleBarcodeF.get(i).length())) : (int)editMaxB) : 0)), read.readR[3],
														hasReversedBarcode ? sampleBarcodeR.get(i) : UtilMethods.complement(sampleBarcodeF.get(i)), null, probB, minOverlapB, wildcard);
											}
											for(int ri = 0; ri < rMatches.size(); ri++){
												for(int rj = 0; rj < sampleEnzymesR.get(i).size(); rj++){
													ArrayList<Match> rEnzymeMatches = null;
													if(probB < 0.0){
														rEnzymeMatches = UtilMethods.searchWithN(read.readR[1], rMatches.get(ri).end + 1/* + randUMILength*/, Math.min(read.readR[1].length(), maxOffsetB + rMatches.get(ri).end + 1 + /*randUMILength + */sampleEnzymesR.get(i).get(rj).length() +
																(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * sampleEnzymesR.get(i).get(rj).length()) : (int)editMaxB) : 0)), sampleEnzymesR.get(i).get(rj), editMaxB, allowIndelsB, true, Integer.MAX_VALUE, wildcard, enzymePatternsR.get(i).get(rj));
													}else{
														rEnzymeMatches = UtilMethods.searchWithProb(read.readR[1], rMatches.get(ri).end + 1/* + randUMILength*/, Math.min(read.readR[1].length(), maxOffsetB + rMatches.get(ri).end + 1 + /*randUMILength + */sampleEnzymesR.get(i).get(rj).length() +
																(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * sampleEnzymesR.get(i).get(rj).length()) : (int)editMaxB) : 0)), read.readR[3], sampleEnzymesR.get(i).get(rj), null, probB, minOverlapB, wildcard);
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
										
										if(inputFileR == null || !checkReversedReads || tempEnzymeEnd2 != -1){
											isMatch = true;
											if(matches.get(j).edits <= minEdit && (matches.get(j).edits < minEdit || matches.get(j).length > barcodeLength)){
												//if all checks are passed, then the current barcodes/enzymes are saved as the matched barcodes/enzymes
												//this state may be reached by multiple different match locations or different barcodes/enzymes, but the best location is chosen
												barcodeIndex = i;
												barcodeEnd = matches.get(j).end + 1;
												barcodeLength = matches.get(j).length;
												enzymeEnd = enzymeMatches.get(enzymeMatches.size() - 1).end + 1;
												if(inputFileR != null){
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
						for(int i = 0; i < sampleBarcodeF.size(); i++){
							ArrayList<Match> matches = null;
							if(probB < 0.0){
								matches = UtilMethods.searchWithN(read.readIF[1], 0, Math.min(read.readIF[1].length(), sampleBarcodeF.get(i).length() +
										(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * sampleBarcodeF.get(i).length()) : (int)editMaxB) : 0)),
										sampleBarcodeF.get(i), editMaxB, allowIndelsB, true, minOverlapB, wildcard, barcodePatternsF.get(i));
							}else{
								matches = UtilMethods.searchWithProb(read.readIF[1], 0, Math.min(read.readIF[1].length(), sampleBarcodeF.get(i).length() +
										(probB < 0.0 && allowIndelsB ? (editMaxB < 0.0 ? (int)(-editMaxB * sampleBarcodeF.get(i).length()) : (int)editMaxB) : 0)),
										read.readIF[3], sampleBarcodeF.get(i), null, probB, minOverlapB, wildcard);
							}
							if(!matches.isEmpty()){
								ArrayList<Match> rMatches = null; //check reversed reads for barcode matches
								if(inputFileR != null && checkReversedReads){
									if(probB < 0.0){
										rMatches = UtilMethods.searchWithN(read.readIR[1], randUMILength, Math.min(read.readIR[1].length(), randUMILength +
												(hasReversedBarcode ? sampleBarcodeR.get(i).length() : sampleBarcodeF.get(i).length()) + (probB < 0.0 && allowIndelsB ?
														(editMaxB < 0.0 ? (int)(-editMaxB * (hasReversedBarcode ? sampleBarcodeR.get(i).length() : sampleBarcodeF.get(i).length())) : (int)editMaxB) : 0)),
												hasReversedBarcode ? sampleBarcodeR.get(i) : UtilMethods.complement(sampleBarcodeF.get(i)), editMaxB, allowIndelsB, true, Integer.MAX_VALUE, wildcard, barcodePatternsR.get(i));
									}else{
										rMatches = UtilMethods.searchWithProb(read.readIR[1], randUMILength, Math.min(read.readIR[1].length(), randUMILength +
												(hasReversedBarcode ? sampleBarcodeR.get(i).length() : sampleBarcodeF.get(i).length()) + (probB < 0.0 && allowIndelsB ?
														(editMaxB < 0.0 ? (int)(-editMaxB * (hasReversedBarcode ? sampleBarcodeR.get(i).length() : sampleBarcodeF.get(i).length())) : (int)editMaxB) : 0)),
												read.readIR[3], hasReversedBarcode ? sampleBarcodeR.get(i) : UtilMethods.complement(sampleBarcodeF.get(i)), null, probB, minOverlapB, wildcard);
									}
								}
								if((inputFileR == null || !checkReversedReads || !rMatches.isEmpty()) && matches.get(matches.size() - 1).edits <= minEdit && (matches.get(matches.size() - 1).edits < minEdit || matches.get(matches.size() - 1).end + 1 > barcodeEnd)){
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
				}
				
				if(sampleInfoFile == null || barcodeIndex != -1){ //no sample file or if barcodes and enzymes were found in the read
					//check if the quality is good enough
					boolean qualityAcceptable;
					if(filterAlgorithm){
						qualityAcceptable = UtilMethods.toError(read.readF[3], 0) <= qualityFilter && (inputFileR == null || UtilMethods.toError(read.readR[3], 0) <= qualityFilter);
					}else{
						qualityAcceptable = UtilMethods.toQScore(read.readF[3], 0) >= qualityFilter && (inputFileR == null || UtilMethods.toQScore(read.readR[3], 0) >= qualityFilter);
					}
					
					if(qualityAcceptable){
						//remove barcode and enzyme
						if(sampleInfoFile != null){
							if(indexFileF == null){
								if(removeEnzyme && removeBarRand){
									read.readF[1] = read.readF[1].substring(enzymeEnd);
									read.readF[3] = read.readF[3].substring(enzymeEnd);
								}else if(removeEnzyme){
									read.readF[1] = read.readF[1].substring(0, barcodeEnd/* + randUMILength*/) + read.readF[1].substring(enzymeEnd);
									read.readF[3] = read.readF[3].substring(0, barcodeEnd/* + randUMILength*/) + read.readF[3].substring(enzymeEnd);
								}else if(removeBarRand){
									read.readF[1] = read.readF[1].substring(barcodeEnd/* + randUMILength*/);
									read.readF[3] = read.readF[3].substring(barcodeEnd/* + randUMILength*/);
								}
							}
							
							if(inputFileR != null){
								if(indexFileR == null){
									if(removeEnzyme && removeBarRand){
										read.readR[1] = read.readR[1].substring(enzymeEnd2);
										read.readR[3] = read.readR[3].substring(enzymeEnd2);
									}else if(removeEnzyme){
										read.readR[1] = read.readR[1].substring(0, barcodeEnd2) + read.readR[1].substring(enzymeEnd2);
										read.readR[3] = read.readR[3].substring(0, barcodeEnd2) + read.readR[3].substring(enzymeEnd2);
									}else if(removeBarRand){
										read.readR[1] = read.readR[1].substring(barcodeEnd2);
										read.readR[3] = read.readR[3].substring(barcodeEnd2);
									}
								}
							}
						}
						
						//these temp arrays are needed to check if a read is trimmed
						String[] temp;
						String[] qualityTrimmed;
						
						int trimmedQuality = 0;
						int removedAdapter = 0;
						
						//trim N
						temp = UtilMethods.trimN(read.readF[1], read.readF[3], trimNPercent);
						read.readF[1] = temp[0];
						read.readF[3] = temp[1];
						
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
						String[] removedAdapters = UtilMethods.removeAdapters(qualityTrimmed[0], qualityTrimmed[1], adaptersF, editMaxA, minOverlapA, maxOffsetA, allowIndelsA, probA, wildcard, adapterPatternsF);
						if(qualityTrimmed[0].length() != removedAdapters[0].length()){
							removedAdapter++;
						}
						read.readF[1] = removedAdapters[0];
						read.readF[3] = removedAdapters[1];
						
						//do the same for reversed reads
						if(inputFileR != null){
							temp = UtilMethods.trimN(read.readR[1], read.readR[3], trimNPercent);
							read.readR[1] = temp[0];
							read.readR[3] = temp[1];
							
							if(trimAlgorithm){
								temp = UtilMethods.qualityTrim2(read.readR[1], read.readR[3], qualityTrimQScore1, true, qualityTrimLength);
								qualityTrimmed = UtilMethods.qualityTrim2(temp[0], temp[1], qualityTrimQScore2, false, qualityTrimLength);
							}else{
								temp = UtilMethods.qualityTrim1(read.readR[1], read.readR[3], qualityTrimQScore1, true, qualityTrimLength);
								qualityTrimmed = UtilMethods.qualityTrim1(temp[0], temp[1], qualityTrimQScore2, false, qualityTrimLength);
							}
							if(read.readR[1].length() != qualityTrimmed[0].length()){
								trimmedQuality++;
							}
							removedAdapters = UtilMethods.removeAdapters(qualityTrimmed[0], qualityTrimmed[1], adaptersR, editMaxA, minOverlapA, maxOffsetA, allowIndelsA, probA, wildcard, adapterPatternsR);
							if(qualityTrimmed[0].length() != removedAdapters[0].length()){
								removedAdapter++;
							}
							read.readR[1] = removedAdapters[0];
							read.readR[3] = removedAdapters[1];
						}
						
						boolean mergedReads = false;
						//merge paired-end reads
						if(mergePairedEnds && inputFileR != null){
							String[] merged = UtilMethods.mergeReads(read.readF[1], read.readF[3], read.readR[1], read.readR[3], editMaxM, probM, minOverlapM, wildcard);
							if(merged[0].length() != read.readF[1].length() + read.readR[1].length()){
								mergedReads = true;
							}
							read.readF[1] = merged[0];
							read.readF[3] = merged[1];
						}
						
						//check if reads are trimmed or merged and if length is too long or short
						if((!removeUntrimmedReads || trimmedQuality == (inputFileR == null ? 1 : 2)) && (!removeNoAdapterReads || removedAdapter == (inputFileR == null ? 1 : 2)) && (!removeUnmergedReads || mergedReads) &&
								minLength <= read.readF[1].length() && read.readF[1].length() <= maxLength && (inputFileR == null || mergePairedEnds || (minLength <= read.readR[1].length() && read.readR[1].length() <= maxLength))){
							//statistics
							readCounts[sampleInfoFile == null ? 0 : barcodeIndex][0].add(inputFileR == null ? 1 : 2);
							
							readCounts[sampleInfoFile == null ? 0 : barcodeIndex][1].add(read.readF[1].length() + (inputFileR == null || mergePairedEnds ? 0 : read.readR[1].length()));
							baseCount.add(read.readF[1].length() + (inputFileR == null || mergePairedEnds ? 0 : read.readR[1].length()));
							
							readCounts[sampleInfoFile == null ? 0 : barcodeIndex][2].add(mergedReads ? 2 : 0);
							totalReadsMerged.add(mergedReads ? 2 : 0);
							
							readCounts[sampleInfoFile == null ? 0 : barcodeIndex][3].add(removedAdapter);
							totalRemovedAdapters.add(removedAdapter);
							
							readCounts[sampleInfoFile == null ? 0 : barcodeIndex][4].add(trimmedQuality);
							totalQualityTrimmed.add(trimmedQuality);
							
							synchronized(locks[sampleInfoFile == null ? 0 : barcodeIndex]){ //synchronized output per file, so no two threads can write to the same file at the same time
								try{
									Strings strings = new Strings();
									//initialize file writers if they are not already initialized
									//this lazy initialization will not create files for samples whose barcodes do not exist in any reads
									if(sampleInfoFile == null){
										if(writersF[0] == null){ //forwards
											if(!outputGZIP && !keepFirstDup && !keepBestDup){
												writersF[0] = new BufferedWriter(new FileWriter(outputDir + "reads_kept_R1.fastq"), BUFFER_SIZE);
												strings.add(outputDir + "reads_kept_R1.fastq");
												
												if(indexFileF != null){
													writersIF[0] = new BufferedWriter(new FileWriter(outputDir + "index_kept_R1.fastq"), BUFFER_SIZE);
													strings.add(outputDir + "index_kept_R1.fastq");
												}
											}else{
												writersF[0] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(((keepFirstDup || keepBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "reads_kept_R1.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
												strings.add(((keepFirstDup || keepBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "reads_kept_R1.fastq.gz");
												
												if(indexFileF != null){
													writersIF[0] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(((keepFirstDup || keepBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "index_kept_R1.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
													strings.add(((keepFirstDup || keepBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "index_kept_R1.fastq.gz");
												}
											}
										}
										
										if(!mergePairedEnds && inputFileR != null && writersR[0] == null){ //reversed
											if(!outputGZIP && !keepFirstDup && !keepBestDup){
												writersR[0] = new BufferedWriter(new FileWriter(outputDir + "reads_kept_R2.fastq"), BUFFER_SIZE);
												strings.add(outputDir + "reads_kept_R2.fastq");
												
												if(indexFileR != null){
													writersIR[0] = new BufferedWriter(new FileWriter(outputDir + "index_kept_R2.fastq"), BUFFER_SIZE);
													strings.add(outputDir + "index_kept_R2.fastq");
												}
											}else{
												writersR[0] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(((keepFirstDup || keepBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "reads_kept_R2.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
												strings.add(((keepFirstDup || keepBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "reads_kept_R2.fastq.gz");
												
												if(indexFileR != null){
													writersIR[0] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(((keepFirstDup || keepBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "index_kept_R2.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
													strings.add(((keepFirstDup || keepBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "index_kept_R2.fastq.gz");
												}
											}
										}
									}else{
										if(writersF[barcodeIndex] == null){ //forwards
											if(!outputGZIP && !keepFirstDup && !keepBestDup){
												writersF[barcodeIndex] = new BufferedWriter(new FileWriter(outputDir + "sample_" + sampleMapF.get(sampleBarcodeF.get(barcodeIndex)) + "_R1.fastq"), BUFFER_SIZE);
												strings.add(outputDir + "sample_" + sampleMapF.get(sampleBarcodeF.get(barcodeIndex)) + "_R1.fastq");
												
												if(indexFileF != null){
													writersIF[barcodeIndex] = new BufferedWriter(new FileWriter(outputDir + "index_" + sampleMapF.get(sampleBarcodeF.get(barcodeIndex)) + "_R1.fastq"), BUFFER_SIZE);
													strings.add(outputDir + "index_" + sampleMapF.get(sampleBarcodeF.get(barcodeIndex)) + "_R1.fastq");
												}
											}else{
												writersF[barcodeIndex] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(((keepFirstDup || keepBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "sample_" + sampleMapF.get(sampleBarcodeF.get(barcodeIndex)) + "_R1.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
												strings.add(((keepFirstDup || keepBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "sample_" + sampleMapF.get(sampleBarcodeF.get(barcodeIndex)) + "_R1.fastq.gz");
												
												if(indexFileF != null){
													writersIF[barcodeIndex] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(((keepFirstDup || keepBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "index_" + sampleMapF.get(sampleBarcodeF.get(barcodeIndex)) + "_R1.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
													strings.add(((keepFirstDup || keepBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "index_" + sampleMapF.get(sampleBarcodeF.get(barcodeIndex)) + "_R1.fastq.gz");
												}
											}
										}
										
										if(!mergePairedEnds && inputFileR != null && writersR[barcodeIndex] == null){ //reversed
											if(!outputGZIP && !keepFirstDup && !keepBestDup){
												writersR[barcodeIndex] = new BufferedWriter(new FileWriter(outputDir + "sample_" + sampleMapF.get(sampleBarcodeF.get(barcodeIndex)) + "_R2.fastq"), BUFFER_SIZE);
												strings.add(outputDir + "sample_" + sampleMapF.get(sampleBarcodeF.get(barcodeIndex)) + "_R2.fastq");
												
												if(indexFileR != null){
													writersIR[barcodeIndex] = new BufferedWriter(new FileWriter(outputDir + "index_" + sampleMapF.get(sampleBarcodeF.get(barcodeIndex)) + "_R2.fastq"), BUFFER_SIZE);
													strings.add(outputDir + "index_" + sampleMapF.get(sampleBarcodeF.get(barcodeIndex)) + "_R2.fastq");
												}
											}else{
												writersR[barcodeIndex] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(((keepFirstDup || keepBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "sample_" + sampleMapF.get(sampleBarcodeF.get(barcodeIndex)) + "_R2.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
												strings.add(((keepFirstDup || keepBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "sample_" + sampleMapF.get(sampleBarcodeF.get(barcodeIndex)) + "_R2.fastq.gz");
												
												if(indexFileR != null){
													writersIR[barcodeIndex] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(((keepFirstDup || keepBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "index_" + sampleMapF.get(sampleBarcodeF.get(barcodeIndex)) + "_R2.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
													strings.add(((keepFirstDup || keepBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "index_" + sampleMapF.get(sampleBarcodeF.get(barcodeIndex)) + "_R2.fastq.gz");
												}
											}
										}
									}
									
									if(strings.size() > 0) //if a new writer was created
										files.offer(strings);
									
									//print to the whatever file the read belongs to
									writersF[sampleInfoFile == null ? 0 : barcodeIndex].write(read.readF[0]);
									writersF[sampleInfoFile == null ? 0 : barcodeIndex].newLine();
									writersF[sampleInfoFile == null ? 0 : barcodeIndex].write(read.readF[1]);
									writersF[sampleInfoFile == null ? 0 : barcodeIndex].newLine();
									writersF[sampleInfoFile == null ? 0 : barcodeIndex].write(description2);
									writersF[sampleInfoFile == null ? 0 : barcodeIndex].newLine();
									writersF[sampleInfoFile == null ? 0 : barcodeIndex].write(read.readF[3]);
									writersF[sampleInfoFile == null ? 0 : barcodeIndex].newLine();
									
									if(indexFileF != null){
										writersIF[sampleInfoFile == null ? 0 : barcodeIndex].write(read.readIF[0]);
										writersIF[sampleInfoFile == null ? 0 : barcodeIndex].newLine();
										writersIF[sampleInfoFile == null ? 0 : barcodeIndex].write(read.readIF[1]);
										writersIF[sampleInfoFile == null ? 0 : barcodeIndex].newLine();
										writersIF[sampleInfoFile == null ? 0 : barcodeIndex].write(description2);
										writersIF[sampleInfoFile == null ? 0 : barcodeIndex].newLine();
										writersIF[sampleInfoFile == null ? 0 : barcodeIndex].write(read.readIF[3]);
										writersIF[sampleInfoFile == null ? 0 : barcodeIndex].newLine();
									}
									
									if(!mergePairedEnds && inputFileR != null){
										writersR[sampleInfoFile == null ? 0 : barcodeIndex].write(read.readR[0]);
										writersR[sampleInfoFile == null ? 0 : barcodeIndex].newLine();
										writersR[sampleInfoFile == null ? 0 : barcodeIndex].write(read.readR[1]);
										writersR[sampleInfoFile == null ? 0 : barcodeIndex].newLine();
										writersR[sampleInfoFile == null ? 0 : barcodeIndex].write(description2);
										writersR[sampleInfoFile == null ? 0 : barcodeIndex].newLine();
										writersR[sampleInfoFile == null ? 0 : barcodeIndex].write(read.readR[3]);
										writersR[sampleInfoFile == null ? 0 : barcodeIndex].newLine();
										
										if(indexFileR != null){
											writersIR[sampleInfoFile == null ? 0 : barcodeIndex].write(read.readIR[0]);
											writersIR[sampleInfoFile == null ? 0 : barcodeIndex].newLine();
											writersIR[sampleInfoFile == null ? 0 : barcodeIndex].write(read.readIR[1]);
											writersIR[sampleInfoFile == null ? 0 : barcodeIndex].newLine();
											writersIR[sampleInfoFile == null ? 0 : barcodeIndex].write(description2);
											writersIR[sampleInfoFile == null ? 0 : barcodeIndex].newLine();
											writersIR[sampleInfoFile == null ? 0 : barcodeIndex].write(read.readIR[3]);
											writersIR[sampleInfoFile == null ? 0 : barcodeIndex].newLine();
										}
									}
								}catch(Exception e){
									UtilMethods.defaultExceptionHandler(logWriter, e);
								}
							}
							valid = true;
						}
					}
				}
			}
			
			if(!valid){
				//statistics
				readCounts[readCounts.length - 1][0].add(inputFileR == null ? 1 : 2);
				undeterminedReads.add(inputFileR == null ? 1 : 2);
				
				readCounts[readCounts.length - 1][1].add(tempSequenceF.length() + (inputFileR == null ? 0 : tempSequenceR.length()));
				baseCount.add(tempSequenceF.length() + (inputFileR == null ? 0 : tempSequenceR.length()));
				
				synchronized(locks[locks.length - 1]){
					try{
						if(undeterminedWriterF[0] == null){ //initialize writers for the undetermined files if they have not already been initialized
							if(outputGZIP){
								undeterminedWriterF[0] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + (sampleInfoFile == null ? "reads" : "sample") + "_undetermined_R1.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
								if(inputFileR != null)
									undeterminedWriterR[0] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + (sampleInfoFile == null ? "reads" : "sample") + "_undetermined_R2.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
								if(indexFileF != null){
									undeterminedWriterF[1] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "index_undetermined_R1.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
									if(inputFileR != null)
										undeterminedWriterR[1] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "index_undetermined_R2.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
								}
							}else{
								undeterminedWriterF[0] = new BufferedWriter(new FileWriter(outputDir + (sampleInfoFile == null ? "reads" : "sample") + "_undetermined_R1.fastq"), BUFFER_SIZE);
								if(inputFileR != null)
									undeterminedWriterR[0] = new BufferedWriter(new FileWriter(outputDir + (sampleInfoFile == null ? "reads" : "sample") + "_undetermined_R2.fastq"), BUFFER_SIZE);
								if(indexFileF != null){
									undeterminedWriterF[1] = new BufferedWriter(new FileWriter(outputDir + "index_undetermined_R1.fastq"), BUFFER_SIZE);
									if(inputFileR != null)
										undeterminedWriterR[1] = new BufferedWriter(new FileWriter(outputDir + "index_undetermined_R2.fastq"), BUFFER_SIZE);
								}
							}
							generatedUndeterminedFile = true;
						}
						
						//print to undetermined file
						undeterminedWriterF[0].write(read.readF[0]);
						undeterminedWriterF[0].newLine();
						undeterminedWriterF[0].write(tempSequenceF);
						undeterminedWriterF[0].newLine();
						undeterminedWriterF[0].write(description2);
						undeterminedWriterF[0].newLine();
						undeterminedWriterF[0].write(tempQualityF);
						undeterminedWriterF[0].newLine();
						
						if(indexFileF != null){
							undeterminedWriterF[1].write(read.readIF[0]);
							undeterminedWriterF[1].newLine();
							undeterminedWriterF[1].write(read.readIF[1]);
							undeterminedWriterF[1].newLine();
							undeterminedWriterF[1].write(description2);
							undeterminedWriterF[1].newLine();
							undeterminedWriterF[1].write(read.readIF[3]);
							undeterminedWriterF[1].newLine();
						}
						
						if(inputFileR != null){
							undeterminedWriterR[0].write(read.readR[0]);
							undeterminedWriterR[0].newLine();
							undeterminedWriterR[0].write(tempSequenceR);
							undeterminedWriterR[0].newLine();
							undeterminedWriterR[0].write(description2);
							undeterminedWriterR[0].newLine();
							undeterminedWriterR[0].write(tempQualityR);
							undeterminedWriterR[0].newLine();
							
							if(indexFileR != null){
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
		readerF.close();
		if(undeterminedWriterF[0] != null)
			undeterminedWriterF[0].close();
		if(indexFileF != null){
			readerIF.close();
			if(undeterminedWriterF[1] != null)
				undeterminedWriterF[1].close();
		}
		for(int i = 0; i < writersF.length; i++){
			if(writersF[i] != null)
				writersF[i].close();
			if(indexFileF != null && writersIF[i] != null)
				writersIF[i].close();
		}
		if(inputFileR != null){
			readerR.close();
			if(undeterminedWriterR != null && undeterminedWriterR[0] != null)
				undeterminedWriterR[0].close();
			if(indexFileR != null){
				readerIR.close();
				if(undeterminedWriterR != null && undeterminedWriterR[1] != null)
					undeterminedWriterR[1].close();
			}
			if(!mergePairedEnds){
				for(int i = 0; i < writersR.length; i++){
					if(writersR[i] != null)
						writersR[i].close();
					if(indexFileR != null && writersIR[i] != null)
						writersIR[i].close();
				}
			}
		}
		
		//update stats
		if(sampleInfoFile == null){
			EnumMap<Stat, String> map = stats.get("Kept");
			map.put(Stat.SEQUENCE_COUNT, DECIMAL_FORMAT.format(readCounts[0][0].sum()));
			map.put(Stat.SEQUENCE_PERCENT, DECIMAL_FORMAT.format((double)readCounts[0][0].sum() / (double)totalReadsProcessed.sum()));
			map.put(Stat.MERGED_COUNT, DECIMAL_FORMAT.format(readCounts[0][2].sum()));
			map.put(Stat.MERGED_PERCENT, DECIMAL_FORMAT.format((double)readCounts[0][2].sum() / (double)totalReadsProcessed.sum()));
			map.put(Stat.REMOVEADAPTER_COUNT, DECIMAL_FORMAT.format(readCounts[0][3].sum()));
			map.put(Stat.REMOVEADAPTER_PERCENT, DECIMAL_FORMAT.format((double)readCounts[0][3].sum() / (double)readCounts[0][0].sum()));
			map.put(Stat.QUALITYTRIM_COUNT, DECIMAL_FORMAT.format(readCounts[0][4].sum()));
			map.put(Stat.QUALITYTRIM_PERCENT, DECIMAL_FORMAT.format((double)readCounts[0][4].sum() / (double)readCounts[0][0].sum()));
			map.put(Stat.BASE_COUNT, DECIMAL_FORMAT.format(readCounts[0][1].sum()));
			map.put(Stat.BASE_PERCENT, DECIMAL_FORMAT.format((double)readCounts[0][1].sum() / (double)baseCount.sum()));
		}else{
			for(int i = 0; i < sampleBarcodeF.size(); i++){
				EnumMap<Stat, String> map = stats.get(sampleMapF.get(sampleBarcodeF.get(i)));
				map.put(Stat.SEQUENCE_COUNT, DECIMAL_FORMAT.format(readCounts[i][0].sum()));
				map.put(Stat.SEQUENCE_PERCENT, DECIMAL_FORMAT.format((double)readCounts[i][0].sum() / (double)totalReadsProcessed.sum()));
				map.put(Stat.MERGED_COUNT, DECIMAL_FORMAT.format(readCounts[i][2].sum()));
				map.put(Stat.MERGED_PERCENT, DECIMAL_FORMAT.format((double)readCounts[i][2].sum() / (double)totalReadsProcessed.sum()));
				map.put(Stat.REMOVEADAPTER_COUNT, DECIMAL_FORMAT.format(readCounts[i][3].sum()));
				map.put(Stat.REMOVEADAPTER_PERCENT, DECIMAL_FORMAT.format((double)readCounts[i][3].sum() / (double)readCounts[i][0].sum()));
				map.put(Stat.QUALITYTRIM_COUNT, DECIMAL_FORMAT.format(readCounts[i][4].sum()));
				map.put(Stat.QUALITYTRIM_PERCENT, DECIMAL_FORMAT.format((double)readCounts[i][4].sum() / (double)readCounts[i][0].sum()));
				map.put(Stat.BASE_COUNT, DECIMAL_FORMAT.format(readCounts[i][1].sum()));
				map.put(Stat.BASE_PERCENT, DECIMAL_FORMAT.format((double)readCounts[i][1].sum() / (double)baseCount.sum()));
			}
		}
		stats.get("Undetermined").put(Stat.SEQUENCE_COUNT, DECIMAL_FORMAT.format(readCounts[readCounts.length - 1][0].sum()));
		stats.get("Undetermined").put(Stat.SEQUENCE_PERCENT, DECIMAL_FORMAT.format((double)readCounts[readCounts.length - 1][0].sum() / (double)totalReadsProcessed.sum()));
		stats.get("Undetermined").put(Stat.BASE_COUNT, DECIMAL_FORMAT.format(readCounts[readCounts.length - 1][1].sum()));
		stats.get("Undetermined").put(Stat.BASE_PERCENT, DECIMAL_FORMAT.format((double)readCounts[readCounts.length - 1][1].sum() / (double)baseCount.sum()));
		
		return files;
	}
	
	private long deduplicate(String readPathF, String readPathR, String indexPathF, String indexPathR, String outDir, String dupPath1, String dupPath2) throws Exception{
		long removed = 0L;
		
		File originalFileF = new File(readPathF);
		int pos = originalFileF.getName().indexOf(".");
		String justNameF = pos > 0 ? originalFileF.getName().substring(0, pos) : originalFileF.getName();
		File outputFileF = new File(outputDir + justNameF + "_dedup.fastq" + (outputGZIP ? ".gz" : ""));
		
		File originalFileIF = new File(indexPathF);
		pos = originalFileIF.getName().indexOf(".");
		String justNameIF = pos > 0 ? originalFileIF.getName().substring(0, pos) : originalFileIF.getName();
		File outputFileIF = new File(outputDir + justNameIF + "_dedup.fastq" + (outputGZIP ? ".gz" : ""));
		
		File originalFileR = null;
		String justNameR = null;
		File outputFileR = null;
		
		File originalFileIR = null;
		String justNameIR = null;
		File outputFileIR = null;
		
		if(readPathR != null){
			originalFileR = new File(readPathR);
			pos = originalFileR.getName().indexOf(".");
			justNameR = pos > 0 ? originalFileR.getName().substring(0, pos) : originalFileR.getName();
			outputFileR = new File(outputDir + justNameR + "_dedup.fastq" + (outputGZIP ? ".gz" : ""));
			
			originalFileIR = new File(indexPathR);
			pos = originalFileIR.getName().indexOf(".");
			justNameIR = pos > 0 ? originalFileIR.getName().substring(0, pos) : originalFileIR.getName();
			outputFileIR = new File(outputDir + justNameIR + "_dedup.fastq" + (outputGZIP ? ".gz" : ""));
		}
		
		BufferedReader readerF;
		BufferedReader readerR = null;
		BufferedReader readerIF;
		BufferedReader readerIR = null;
		readerF = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(readPathF), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
		readerIF = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(indexPathF), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
		if(readPathR != null){
			readerR = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(readPathR), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
			readerIR = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(indexPathR), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
		}
		
		if(keepFirstDup){
			BufferedWriter writerF;
			BufferedWriter writerR = null;
			BufferedWriter writerIF;
			BufferedWriter writerIR = null;
			if(outputGZIP){
				writerF = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFileF.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
				writerIF = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFileIF.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
				if(readPathR != null){
					writerR = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFileR.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
					writerIR = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFileIR.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
				}
			}else{
				writerF = new BufferedWriter(new FileWriter(outputFileF.getAbsolutePath()), BUFFER_SIZE);
				writerIF = new BufferedWriter(new FileWriter(outputFileIF.getAbsolutePath()), BUFFER_SIZE);
				if(readPathR != null){
					writerR = new BufferedWriter(new FileWriter(outputFileR.getAbsolutePath()), BUFFER_SIZE);
					writerIR = new BufferedWriter(new FileWriter(outputFileIR.getAbsolutePath()), BUFFER_SIZE);
				}
			}
			BufferedWriter dupReadWriter = null;
			BufferedWriter dupIndexWriter = null;
			HashSet<BitSet> set = new HashSet<BitSet>();
			String[] linesF = new String[4];
			String[] linesR = new String[4];
			String[] linesIF = new String[4];
			String[] linesIR = new String[4];
			while((linesF[0] = readerF.readLine()) != null && (linesF[1] = readerF.readLine()) != null &&
					(linesF[2] = readerF.readLine()) != null && (linesF[3] = readerF.readLine()) != null &&
					(linesIF[0] = readerIF.readLine()) != null && (linesIF[1] = readerIF.readLine()) != null &&
					(linesIF[2] = readerIF.readLine()) != null && (linesIF[3] = readerIF.readLine()) != null &&
					(readPathR == null || ((linesR[0] = readerR.readLine()) != null && (linesR[1] = readerR.readLine()) != null &&
					(linesR[2] = readerR.readLine()) != null && (linesR[3] = readerR.readLine()) != null &&
					(linesIR[0] = readerIR.readLine()) != null && (linesIR[1] = readerIR.readLine()) != null &&
					(linesIR[2] = readerIR.readLine()) != null && (linesIR[3] = readerIR.readLine()) != null))){
				if(linesIF[1].length() < randUMILength)
					continue;
				BitSet key = UtilMethods.toBit(linesIF[1].substring(linesIF[1].length() - randUMILength, linesIF[1].length()));
				if(set.contains(key)){ //if a read with the same UMI has already been encountered, then the current read is removed
					duplicatesRemoved += readPathR == null ? 1 : 2;
					removed += readPathR == null ? 1 : 2;
					if(duplicatesRemoved % printDuplicateInterval == 0){
						logWriter.println("Deduplicated So Far: " + DECIMAL_FORMAT.format(duplicatesRemoved));
						logWriter.println("Run Time So Far: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
						logWriter.println("Current Memory Usage (GB): " + DECIMAL_FORMAT.format(UtilMethods.currentMemoryUsage()));
						logWriter.flush();
					}
					if(saveDup){ //print removed duplicate to another file if needed
						if(dupReadWriter == null){ //initialize duplicate writer if it is not initialized
							dupReadWriter = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(dupPath1), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
							dupIndexWriter = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(dupPath2), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
							generatedDupFiles += 2;
						}
						dupReadWriter.write(linesF[0]);
						dupReadWriter.newLine();
						dupReadWriter.write(linesF[1]);
						dupReadWriter.newLine();
						dupReadWriter.write(description2);
						dupReadWriter.newLine();
						dupReadWriter.write(linesF[3]);
						dupReadWriter.newLine();
						
						dupIndexWriter.write(linesIF[0]);
						dupIndexWriter.newLine();
						dupIndexWriter.write(linesIF[1]);
						dupIndexWriter.newLine();
						dupIndexWriter.write(description2);
						dupIndexWriter.newLine();
						dupIndexWriter.write(linesIF[3]);
						dupIndexWriter.newLine();
						
						if(readPathR != null){
							dupReadWriter.write(linesR[0]);
							dupReadWriter.newLine();
							dupReadWriter.write(linesR[1]);
							dupReadWriter.newLine();
							dupReadWriter.write(description2);
							dupReadWriter.newLine();
							dupReadWriter.write(linesR[3]);
							dupReadWriter.newLine();
							
							dupIndexWriter.write(linesIR[0]);
							dupIndexWriter.newLine();
							dupIndexWriter.write(linesIR[1]);
							dupIndexWriter.newLine();
							dupIndexWriter.write(description2);
							dupIndexWriter.newLine();
							dupIndexWriter.write(linesIR[3]);
							dupIndexWriter.newLine();
						}
					}
				}else{ //if the UMI is not encountered before, then the read is saved
					writerF.write(linesF[0]);
					writerF.newLine();
					writerF.write(linesF[1]);
					writerF.newLine();
					writerF.write(description2);
					writerF.newLine();
					writerF.write(linesF[3]);
					writerF.newLine();
					
					writerIF.write(linesIF[0]);
					writerIF.newLine();
					writerIF.write(linesIF[1]);
					writerIF.newLine();
					writerIF.write(description2);
					writerIF.newLine();
					writerIF.write(linesIF[3]);
					writerIF.newLine();
					
					if(readPathR != null){
						writerR.write(linesR[0]);
						writerR.newLine();
						writerR.write(linesR[1]);
						writerR.newLine();
						writerR.write(description2);
						writerR.newLine();
						writerR.write(linesR[3]);
						writerR.newLine();
						
						writerIR.write(linesIR[0]);
						writerIR.newLine();
						writerIR.write(linesIR[1]);
						writerIR.newLine();
						writerIR.write(description2);
						writerIR.newLine();
						writerIR.write(linesIR[3]);
						writerIR.newLine();
					}
					
					set.add(key);
				}
			}
			writerF.close();
			writerIF.close();
			if(readPathR != null){
				writerR.close();
				writerIR.close();
			}
			if(dupReadWriter != null){
				dupReadWriter.close();
				dupIndexWriter.close();
			}
		}else if(keepBestDup){
			BufferedWriter dupReadWriter = null;
			BufferedWriter dupIndexWriter = null;
			
			HashMap<BitSet, CompressedRead> map = new HashMap<BitSet, CompressedRead>();
			String[] linesF = new String[4];
			String[] linesR = new String[4];
			String[] linesIF = new String[4];
			String[] linesIR = new String[4];
			while((linesF[0] = readerF.readLine()) != null && (linesF[1] = readerF.readLine()) != null &&
					(linesF[2] = readerF.readLine()) != null && (linesF[3] = readerF.readLine()) != null &&
					(linesIF[0] = readerIF.readLine()) != null && (linesIF[1] = readerIF.readLine()) != null &&
					(linesIF[2] = readerIF.readLine()) != null && (linesIF[3] = readerIF.readLine()) != null &&
					(readPathR == null || ((linesR[0] = readerR.readLine()) != null && (linesR[1] = readerR.readLine()) != null &&
					(linesR[2] = readerR.readLine()) != null && (linesR[3] = readerR.readLine()) != null &&
					(linesIR[0] = readerIR.readLine()) != null && (linesIR[1] = readerIR.readLine()) != null &&
					(linesIR[2] = readerIR.readLine()) != null && (linesIR[3] = readerIR.readLine()) != null))){
				if(linesIF[1].length() < randUMILength)
					continue;
				BitSet key = UtilMethods.toBit(linesIF[1].substring(linesIF[1].length() - randUMILength, linesIF[1].length()));
				double currQuality = UtilMethods.toError(linesF[3] + (readPathR == null ? "" : linesR[3]), 0);
				if(map.containsKey(key)){
					duplicatesRemoved += readPathR == null ? 1 : 2;
					removed += readPathR == null ? 1 : 2;
					if(duplicatesRemoved % printDuplicateInterval == 0){
						logWriter.println("Deduplicated So Far: " + DECIMAL_FORMAT.format(duplicatesRemoved));
						logWriter.println("Run Time So Far: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
						logWriter.println("Current Memory Usage (GB): " + DECIMAL_FORMAT.format(UtilMethods.currentMemoryUsage()));
						logWriter.flush();
					}
					CompressedRead other = map.get(key);
					if(other.error > currQuality){ //if current read quality is better
						if(saveDup){
							if(dupReadWriter == null){ //initialize writer for writing duplicates
								dupReadWriter = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(dupPath1), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
								dupIndexWriter = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(dupPath2), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
								generatedDupFiles += 2;
							}
							//print other read to duplicates
							dupReadWriter.write(other.descriptionFRead);
							dupReadWriter.newLine();
							dupReadWriter.write(UtilMethods.toSequence(other.sequenceFRead));
							dupReadWriter.newLine();
							dupReadWriter.write(description2);
							dupReadWriter.newLine();
							dupReadWriter.write(UtilMethods.byteArrayToQScore(other.qualityFRead));
							dupReadWriter.newLine();
							
							dupIndexWriter.write(other.descriptionFIndex);
							dupIndexWriter.newLine();
							dupIndexWriter.write(UtilMethods.toSequence(other.sequenceFIndex));
							dupIndexWriter.newLine();
							dupIndexWriter.write(description2);
							dupIndexWriter.newLine();
							dupIndexWriter.write(UtilMethods.byteArrayToQScore(other.qualityFIndex));
							dupIndexWriter.newLine();
							
							if(readPathR != null){
								dupReadWriter.write(other.descriptionRRead);
								dupReadWriter.newLine();
								dupReadWriter.write(UtilMethods.toSequence(other.sequenceRRead));
								dupReadWriter.newLine();
								dupReadWriter.write(description2);
								dupReadWriter.newLine();
								dupReadWriter.write(UtilMethods.byteArrayToQScore(other.qualityRRead));
								dupReadWriter.newLine();
								
								dupIndexWriter.write(other.descriptionRIndex);
								dupIndexWriter.newLine();
								dupIndexWriter.write(UtilMethods.toSequence(other.sequenceRIndex));
								dupIndexWriter.newLine();
								dupIndexWriter.write(description2);
								dupIndexWriter.newLine();
								dupIndexWriter.write(UtilMethods.byteArrayToQScore(other.qualityRIndex));
								dupIndexWriter.newLine();
							}
						}
						//save current read in HashMap
						if(readPathR == null)
							map.put(key, new CompressedRead(currQuality, linesF[0], linesF[1], linesF[3], linesIF[0], linesIF[1], linesIF[3]));
						else
							map.put(key, new CompressedRead(currQuality, linesF[0], linesF[1], linesF[3], linesIF[0], linesIF[1], linesIF[3], linesR[0], linesR[1], linesR[3], linesIR[0], linesIR[1], linesIR[3]));
					}else{ //remove current read
						if(saveDup){
							if(dupReadWriter == null){
								dupReadWriter = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(dupPath1), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
								dupIndexWriter = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(dupPath2), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
								generatedDupFiles += 2;
							}
							dupReadWriter.write(linesF[0]);
							dupReadWriter.newLine();
							dupReadWriter.write(linesF[1]);
							dupReadWriter.newLine();
							dupReadWriter.write(description2);
							dupReadWriter.newLine();
							dupReadWriter.write(linesF[3]);
							dupReadWriter.newLine();
							
							dupIndexWriter.write(linesIF[0]);
							dupIndexWriter.newLine();
							dupIndexWriter.write(linesIF[1]);
							dupIndexWriter.newLine();
							dupIndexWriter.write(description2);
							dupIndexWriter.newLine();
							dupIndexWriter.write(linesIF[3]);
							dupIndexWriter.newLine();
							
							if(readPathR != null){
								dupReadWriter.write(linesR[0]);
								dupReadWriter.newLine();
								dupReadWriter.write(linesR[1]);
								dupReadWriter.newLine();
								dupReadWriter.write(description2);
								dupReadWriter.newLine();
								dupReadWriter.write(linesR[3]);
								dupReadWriter.newLine();
								
								dupIndexWriter.write(linesIR[0]);
								dupIndexWriter.newLine();
								dupIndexWriter.write(linesIR[1]);
								dupIndexWriter.newLine();
								dupIndexWriter.write(description2);
								dupIndexWriter.newLine();
								dupIndexWriter.write(linesIR[3]);
								dupIndexWriter.newLine();
							}
						}
					}
				}else{ //if the current UMI was not encountered before
					if(readPathR == null)
						map.put(key, new CompressedRead(currQuality, linesF[0], linesF[1], linesF[3], linesIF[0], linesIF[1], linesIF[3]));
					else
						map.put(key, new CompressedRead(currQuality, linesF[0], linesF[1], linesF[3], linesIF[0], linesIF[1], linesIF[3], linesR[0], linesR[1], linesR[3], linesIR[0], linesIR[1], linesIR[3]));
				}
			}
			readerF.close();
			readerIF.close();
			if(readPathR != null){
				readerR.close();
				readerIR.close();
			}
			
			BufferedWriter writerF;
			BufferedWriter writerR = null;
			BufferedWriter writerIF;
			BufferedWriter writerIR = null;
			if(outputGZIP){
				writerF = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFileF.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
				writerIF = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFileIF.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
				if(readPathR != null){
					writerR = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFileR.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
					writerIR = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFileIR.getAbsolutePath()), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
				}
			}else{
				writerF = new BufferedWriter(new FileWriter(outputFileF.getAbsolutePath()), BUFFER_SIZE);
				writerIF = new BufferedWriter(new FileWriter(outputFileIF.getAbsolutePath()), BUFFER_SIZE);
				if(readPathR != null){
					writerR = new BufferedWriter(new FileWriter(outputFileR.getAbsolutePath()), BUFFER_SIZE);
					writerIR = new BufferedWriter(new FileWriter(outputFileIR.getAbsolutePath()), BUFFER_SIZE);
				}
			}
			
			//go through the best reads and print them
			for(CompressedRead val : map.values()){
				writerF.write(val.descriptionFRead);
				writerF.newLine();
				writerF.write(UtilMethods.toSequence(val.sequenceFRead));
				writerF.newLine();
				writerF.write(description2);
				writerF.newLine();
				writerF.write(UtilMethods.byteArrayToQScore(val.qualityFRead));
				writerF.newLine();
				
				writerIF.write(val.descriptionFIndex);
				writerIF.newLine();
				writerIF.write(UtilMethods.toSequence(val.sequenceFIndex));
				writerIF.newLine();
				writerIF.write(description2);
				writerIF.newLine();
				writerIF.write(UtilMethods.byteArrayToQScore(val.qualityFIndex));
				writerIF.newLine();
				
				if(readPathR != null){
					writerR.write(val.descriptionRRead);
					writerR.newLine();
					writerR.write(UtilMethods.toSequence(val.sequenceRRead));
					writerR.newLine();
					writerR.write(description2);
					writerR.newLine();
					writerR.write(UtilMethods.byteArrayToQScore(val.qualityRRead));
					writerR.newLine();
					
					writerIR.write(val.descriptionRIndex);
					writerIR.newLine();
					writerIR.write(UtilMethods.toSequence(val.sequenceRIndex));
					writerIR.newLine();
					writerIR.write(description2);
					writerIR.newLine();
					writerIR.write(UtilMethods.byteArrayToQScore(val.qualityRIndex));
					writerIR.newLine();
				}
			}
			writerF.close();
			writerIF.close();
			if(readPathR != null){
				writerR.close();
				writerIR.close();
			}
		}
		
		return removed;
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
		if(args[0].equals("-h") || args[0].equals("--help")){
			System.out.println("Help:");
			System.out.println("For more detailed information, visit https://github.com/Daniel-Liu-c0deb0t/FastQParse/wiki/Commands");
			System.out.println("\nOtherwise, here are some common commands:");
			System.out.println("'-r' - Specify one or two files after to either process single end (1 file) or paired end reads (2 files).");
			System.out.println("'-o' - Specify a directory after to indicate the output directory.");
			System.out.println("'-s' - Specify a sample info file that contains barcode information. This will enable demultiplexing.");
			System.out.println("'-m' - Turn on merging paired end reads. Make sure there are two input files, one for forwards reads, and one for reversed reads.");
			System.out.println("'-a' or '-A' or '-z' or '-Z' - Specify adapter sequences after to trim them as forwards 5', forwards 3', reversed 5' or reversed 3' adapters, respectively.");
			System.out.println("'-q' - Specify one number after to quality trim the 3' end using that number as the threshold. Specify two numbers to quality trim both ends, using the first number for the 5' end, and the second number for the 3' end.");
			System.out.println("If specifying more than one parameter per command, separate those parameters with spaces. For example, '-q 10 20' is a valid command.");
		}else if(args[0].equals("--simulate")){
			for(int i = 1; i < args.length; i++){
				if(args[i].equals("-o")){
					outputDir = args[++i];
					if(!outputDir.endsWith(File.separator) && !outputDir.endsWith("/") && !outputDir.endsWith("\\")){
						outputDir += File.separator;
					}
				}else if(args[i].equals("-s")){
					sampleInfoFile = new File(args[++i]);
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
				}else if(args[i].equals("-lB")){
					allowIndelsB = true;
				}else if(args[i].equals("-lA")){
					allowIndelsA = true;
				}else if(args[i].equals("-vA")){
					minOverlapA = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-vM")){
					minOverlapM = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-R")){
					checkReversedReads = true;
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
			boolean clearDir = false;
			for(int i = 0; i < args.length; i++){
				if(args[i].equals("-dF")){
					keepFirstDup = true;
					keepBestDup = false;
				}else if(args[i].equals("-dB")){
					keepBestDup = true;
					keepFirstDup = false;
				}else if(args[i].equals("-u")){
					randUMILength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-N")){
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						maxNPercent = Double.parseDouble(args[++i]);
					}else{
						maxNPercent = -1.0;
					}
				}else if(args[i].equals("-fB")){
					maxOffsetB = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-fA")){
					maxOffsetA = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-Q")){
					qualityFilter = Double.parseDouble(args[++i]);
				}else if(args[i].equals("-r")){
					inputFileF = new File(args[++i]);
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						inputFileR = new File(args[++i]);
					}
				}else if(args[i].equals("-i")){
					indexFileF = new File(args[++i]);
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						indexFileR = new File(args[++i]);
					}
					randUMILength = 12;
				}else if(args[i].equals("-s")){
					sampleInfoFile = new File(args[++i]);
				}else if(args[i].equals("-o")){
					outputDir = args[++i];
					if(!outputDir.endsWith(File.separator) && !outputDir.endsWith("/") && !outputDir.endsWith("\\")){
						outputDir += File.separator;
					}
				}else if(args[i].equals("-gz")){
					outputGZIP = true;
				}else if(args[i].equals("--clear-dir")){
					clearDir = true;
				}else if(args[i].equals("--save-temp")){
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
				}else if(args[i].equals("--save-dup")){
					saveDup = true;
				}else if(args[i].equals("--print-processed")){
					printProcessedInterval = Long.parseLong(args[++i]);
				}else if(args[i].equals("--print-duplicate")){
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
				}else if(args[i].equals("-vA")){
					minOverlapA = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-vB")){
					minOverlapB = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-vM")){
					minOverlapM = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--alt-quality-trim")){
					trimAlgorithm = true;
					qualityTrimLength = Integer.MAX_VALUE;
				}else if(args[i].equals("-q")){
					if(i + 2 >= args.length || args[i + 2].startsWith("-")){
						qualityTrimQScore2 = Integer.parseInt(args[++i]);
					}else{
						qualityTrimQScore1 = Integer.parseInt(args[++i]);
						qualityTrimQScore2 = Integer.parseInt(args[++i]);
					}
				}else if(args[i].equals("--quality-trim-length")){
					qualityTrimLength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("-n")){
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						trimNPercent = Double.parseDouble(args[++i]);
					}else{
						trimNPercent = 0.5;
					}
				}else if(args[i].equals("-lB")){
					allowIndelsB = true;
				}else if(args[i].equals("-lA")){
					allowIndelsA = true;
				}else if(args[i].equals("--alt-quality-filter")){
					filterAlgorithm = true;
				}else if(args[i].equals("-R")){
					checkReversedReads = true;
				}else if(args[i].equals("-w")){
					wildcard = true;
				}else if(args[i].equals("-S")){
					singleBarcodeMatchOnly = true;
				}else if(args[i].equals("-rQ")){
					removeUntrimmedReads = true;
				}else if(args[i].equals("-rA")){
					removeNoAdapterReads = true;
				}else if(args[i].equals("-rM")){
					removeUnmergedReads = true;
				}else if(args[i].equals("-P")){
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
				}else if(args[i].equals("-E")){
					while(i + 2 < args.length && !args[i + 2].startsWith("-")){
						String name = args[++i].toUpperCase();
						String seq = args[++i].toUpperCase();
						if(EnzymeList.enzymes.containsKey(name)){
							EnzymeList.enzymes.get(name).add(seq);
						}else{
							EnzymeList.enzymes.put(name, new ArrayList<String>(Arrays.asList(seq)));
						}
						newEnzymes += name + " = " + seq + ", ";
					}
					newEnzymes.substring(0, newEnzymes.length() - 2);
				}else{
					throw new Exception("Command not supported: " + args[i]);
				}
			}
			boolean isDirClear = false;
			if(outputDir == null){
				outputDir = inputFileF.getParent();
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
			if(saveDup && (keepBestDup || keepFirstDup)){
				File f3 = new File(outputDir + "dup" + File.separatorChar);
				if(!f3.exists())
					f3.mkdirs();
			}
			if(inputFileF.getAbsolutePath().toLowerCase().endsWith(".gz") || inputFileF.getAbsolutePath().toLowerCase().endsWith(".gzip")){
				inputGZIP = true;
			}
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(outputDir + "FastQParse_" + Mode.PROCESS.description2 + ".log"), BUFFER_SIZE_LOG));
			Date date = new Date();
			logWriter.println("Started on: " + DATE_FORMAT.format(date));
			System.out.println("Started on: " + DATE_FORMAT.format(date));
			if(isDirClear){
				logWriter.println("Cleared Directory: " + outputDir);
			}
			new FastQParseMain(Mode.PROCESS);
			
			date = new Date();
			logWriter.println("Ended on: " + DATE_FORMAT.format(date));
			System.out.println("Ended on: " + DATE_FORMAT.format(date));
		}
		if(logWriter != null){
			logWriter.close();
		}
	}
}
