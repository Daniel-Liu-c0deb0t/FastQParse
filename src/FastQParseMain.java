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
	
	private static ArrayList<String> constEnzymesF; //constant enzyme sequence
	private static ArrayList<String> constEnzymesR; //constant enzyme sequence
	private static int randUMILength; //length of random UMI sequence (could be 0)
	private static boolean removeFirstDup = false; //remove duplicates that are not the first
	private static boolean removeBestDup = false; //remove duplicates that are not the best
	private static boolean replaceOriginal = false; //replace original file
	private static int maxOffset = 0; //tries offsets from 0 to maxOffset when matching barcode/enzyme
	private static double removeDNAWithNPercent = 2.0; //if the percentage of DNA that is N is greater than this value then it will be undetermined
	private static double qualityFilter = 0.0; //quality filter threshold
	private static boolean filterAlgorithm = false; //false = average, true = expected error
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
	private static boolean trimAlgorithm = false; //false = 1st method, true = 2nd method
	private static int qualityTrimQScore1 = 0; //quality under this will be trimmed in quality trimming (from start)
	private static int qualityTrimQScore2 = 0; //quality under this will be trimmed in quality trimming (from end)
	private static int qualityTrimLength = 1; //length that is needed for quality trim algorithm 1 and 2
	private static boolean mergeAlgorithm = false; //false = 1st method, true = 2nd method
	private static double trimNPercent = 2.0; //trim leading and trailing 'N' percentage
	private static boolean allowIndels = false; //allow insertions and deletions
	private static boolean adapterAlgorithm = false; //false = 1st method, true = 2nd method
	private static boolean checkReversedReads = false; //whether to check for barcodes/enzymes in reversed reads
	private static boolean wildcard = true; //whether to check for undetermined bp 'N'
	private static boolean singleBarcodeMatchOnly = false; //whether to only allow one barcode match
	private static boolean removeUntrimmedReads = false;
	private static boolean removeNoAdapterReads = false;
	private static boolean removeNoMergeReads = false;
	
	private enum Mode{ //features
		DEMULTIPLEX("Demultiplex", "demultiplex"), DEDUP("Deduplicate Reads", "dedup"),
		PAIRMERGE("Merge Paired-End Reads", "pairmerge"), FILTER("Filter Reads", "filter");
		
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
	private long totalDNAProcessed = 0L; //total DNA read
	private long duplicatesRemoved = 0L; //total duplicates removed
	private long undeterminedDNA = 0L; //total undetermined DNA
	private long totalRemovedAdapters = 0L; //total DNA with adapters removed
	private long totalQualityTrimmed = 0L; //total DNA that has been quality trimmed
	private long baseCount = 0L; //total base count
	private long totalReadsMerged = 0L; //total number of reads that are actually merged
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
			logWriter.println("Maximum Right Offset: " + maxOffset);
			if(removeDNAWithNPercent >= 0 && removeDNAWithNPercent <= 1.0)
				logWriter.println("Remove Reads With % of N Greater Than: " + removeDNAWithNPercent);
			else if(removeDNAWithNPercent < 0)
				logWriter.println("Remove Reads With at Least 1 N");
			else
				logWriter.println("Remove Reads With N: false");
			logWriter.println("Quality Filter Threshold: " + qualityFilter);
			logWriter.println("Quality Filter Algorithm: " + (filterAlgorithm ? "2" : "1"));
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
			logWriter.println("Allow Insertions and Deletions: " + allowIndels);
			logWriter.println("Save Duplicates in Separate Files: " + saveDup);
			logWriter.println("Interval to Print Processed Reads: " + DECIMAL_FORMAT.format(printProcessedInterval));
			logWriter.println("Interval to Print Duplicate Reads: " + DECIMAL_FORMAT.format(printDuplicateInterval));
			logWriter.println("Minimum Read Length: " + minLength);
			logWriter.println("Maximum Read Length: " + maxLength);
			logWriter.println("Merge Paired-End Reads: " + mergePairedEnds);
			logWriter.println("Merge Algorithm: " + (mergeAlgorithm ? "2" : "1"));
			String adaptersFString = adaptersF.toString();
			logWriter.println("Forwards Read Adapters: " + adaptersFString.substring(1, adaptersFString.length() - 1));
			if(inputFile2 != null){
				String adaptersRString = adaptersR.toString();
				logWriter.println("Reversed Read Adapters: " + adaptersRString.substring(1, adaptersRString.length() - 1));
			}
			logWriter.println("Minimum Adapter Overlap: " + minOverlapA);
			logWriter.println("Minimum Barcode Overlap: " + minOverlapB);
			logWriter.println("Remove Adapter Algorithm: " + (adapterAlgorithm ? "2" : "1"));
			logWriter.println("Quality Trim Algorithm: " + (trimAlgorithm ? "2" : "1"));
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
			
			startTime = System.currentTimeMillis();
			int generatedSampleFiles = 0;
			
			//deduplicate after demultiplex
			if(removeFirstDup || removeBestDup){ //assume that index files are provided for deduplicating for simpler logic
				ArrayList<String> files = demultiplexFile(); //demultiplex
				logWriter.println();
				logWriter.println("Demultiplex Completed.");
				logWriter.println();
				logWriter.println("Run Time So Far: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
				logWriter.println();
				logWriter.println("Deduplicating...");
				logWriter.flush();
				for(int i = 0; i < files.size(); i += inputFile2 == null ? 2 : 4){
					File file1 = new File(files.get(i));
					int pos = file1.getName().indexOf(".");
					String justName1 = pos > 0 ? file1.getName().substring(0, pos) : file1.getName();
					
					File file2 = new File(files.get(i + 1));
					pos = file2.getName().indexOf(".");
					String justName2 = pos > 0 ? file2.getName().substring(0, pos) : file2.getName();
					
					logWriter.println();
					logWriter.println("Deduplicating File: " + justName1);
					logWriter.flush();
					//deduplicate
					long removed = deduplicate(files.get(i), inputFile2 == null ? null : files.get(i + 2), files.get(i + 1), inputFile2 == null ? null : files.get(i + 3),
							outputDir, false, false, outputDir + "dup" + File.separatorChar + justName1 + "_dup.fastq.gz", outputDir + "dup" + File.separatorChar + justName2 + "_dup.fastq.gz");
					for(int j = 0; j < sampleDNAF.size(); j++){ //calculate duduplicated reads count
						if(justName1.contains(sampleMapF.get(sampleDNAF.get(j)))){
							map = stats.get(sampleMapF.get(sampleDNAF.get(j)));
							map.put(Stat.DEDUP_COUNT, DECIMAL_FORMAT.format(removed));
							map.put(Stat.DEDUP_PERCENT, DECIMAL_FORMAT.format((double)removed / Double.parseDouble(map.get(Stat.SEQUENCE_COUNT).replace(",", ""))));
							break;
						}
					}
					System.gc();
				}
				generatedSampleFiles = files.size();
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
				generatedSampleFiles = demultiplexFile().size(); //demultiplex
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
			logWriter.println("Number of Lines Processed: " + DECIMAL_FORMAT.format(totalDNAProcessed * 4));
			logWriter.println("Number of Reads Processed: " + DECIMAL_FORMAT.format(totalDNAProcessed));
			logWriter.println("Number of Generated Sample Files: " + generatedSampleFiles);
			logWriter.println("Number of Generated Duplicate Files: " + generatedDupFiles);
			
			//fill in info for total stats
			map = new EnumMap<Stat, String>(Stat.class);
			map.put(Stat.SEQUENCE_FORWARDS, "");
			if(hasReversedBarcode)
				map.put(Stat.SEQUENCE_REVERSED, "");
			map.put(Stat.SEQUENCE_COUNT, DECIMAL_FORMAT.format(totalDNAProcessed));
			map.put(Stat.SEQUENCE_PERCENT, DECIMAL_FORMAT.format(1));
			map.put(Stat.MERGED_COUNT, DECIMAL_FORMAT.format(totalReadsMerged));
			map.put(Stat.MERGED_PERCENT, DECIMAL_FORMAT.format((double)totalReadsMerged / (double)totalDNAProcessed));
			map.put(Stat.DEDUP_COUNT, DECIMAL_FORMAT.format(duplicatesRemoved));
			map.put(Stat.DEDUP_PERCENT, DECIMAL_FORMAT.format((double)duplicatesRemoved / (double)totalDNAProcessed));
			map.put(Stat.REMOVEADAPTER_COUNT, DECIMAL_FORMAT.format(totalRemovedAdapters));
			map.put(Stat.REMOVEADAPTER_PERCENT, DECIMAL_FORMAT.format((double)totalRemovedAdapters / (double)totalDNAProcessed));
			map.put(Stat.QUALITYTRIM_COUNT, DECIMAL_FORMAT.format(totalQualityTrimmed));
			map.put(Stat.QUALITYTRIM_PERCENT, DECIMAL_FORMAT.format((double)totalQualityTrimmed / (double)totalDNAProcessed));
			map.put(Stat.BASE_COUNT, DECIMAL_FORMAT.format(baseCount));
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
			logWriter.println("Number of Lines Processed: " + DECIMAL_FORMAT.format(totalDNAProcessed * 4));
			logWriter.println("Number of Reads Processed: " + DECIMAL_FORMAT.format(totalDNAProcessed));
			logWriter.println("Deduplicated Count: " + DECIMAL_FORMAT.format(duplicatesRemoved));
			logWriter.println("% Deduplicated: " + DECIMAL_FORMAT.format((double)duplicatesRemoved / (double)totalDNAProcessed));
			logWriter.println("Total Number of BP: " + DECIMAL_FORMAT.format(baseCount));
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
			logWriter.println("Merge Algorithm: " + (mergeAlgorithm ? "2" : "1"));
			logWriter.println("Interval to Print Processed Reads: " + DECIMAL_FORMAT.format(printProcessedInterval));
			logWriter.println("Interval to Print Duplicate Reads: " + DECIMAL_FORMAT.format(printDuplicateInterval));
			logWriter.println("Use Wildcard Characters: " + wildcard);
			logWriter.println("Remove Reads That Are Not Merged: " + removeNoMergeReads);
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
			logWriter.println("Number of Lines Processed: " + DECIMAL_FORMAT.format(totalDNAProcessed * 4));
			logWriter.println("Number of Reads Processed: " + DECIMAL_FORMAT.format(totalDNAProcessed));
			logWriter.println("Number of Reads Merged: " + DECIMAL_FORMAT.format(totalReadsMerged));
			logWriter.println("% Merged: " + DECIMAL_FORMAT.format((double)totalReadsMerged / (double)totalDNAProcessed));
			logWriter.println("Number of Reads Removed: " + DECIMAL_FORMAT.format(undeterminedDNA));
			logWriter.println("% Reads Removed: " + DECIMAL_FORMAT.format((double)undeterminedDNA / (double)totalDNAProcessed));
			logWriter.println("Total Number of BP: " + DECIMAL_FORMAT.format(baseCount));
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
			logWriter.println("Allow Insertions and Deletions: " + allowIndels);
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
			logWriter.println("Quality Filter Algorithm: " + (filterAlgorithm ? "2" : "1"));
			String adaptersFString = adaptersF.toString();
			logWriter.println("Read Adapters: " + adaptersFString.substring(1, adaptersFString.length() - 1));
			logWriter.println("Minimum Adapter Overlap: " + minOverlapA);
			logWriter.println("Remove Adapter Algorithm: " + (adapterAlgorithm ? "2" : "1"));
			logWriter.println("Quality Trim Algorithm: " + (trimAlgorithm ? "2" : "1"));
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
			logWriter.println();
			logWriter.println("Filtering...");
			logWriter.println();
			logWriter.flush();
			
			startTime = System.currentTimeMillis();
			filterFile(); //filter the file
			
			logWriter.println();
			logWriter.println("Filtering Completed.");
			logWriter.println();
			logWriter.println("Total Run Time: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
			logWriter.println();
			logWriter.println("Number of Lines Processed: " + DECIMAL_FORMAT.format(totalDNAProcessed * 4));
			logWriter.println("Number of Reads Processed: " + DECIMAL_FORMAT.format(totalDNAProcessed));
			logWriter.println("Number of Reads With Removed Adapter: " + DECIMAL_FORMAT.format(totalRemovedAdapters));
			logWriter.println("% Removed Adapter: " + DECIMAL_FORMAT.format((double)totalRemovedAdapters / (double)totalDNAProcessed));
			logWriter.println("Number of Reads Quality Trimmed: " + DECIMAL_FORMAT.format(totalQualityTrimmed));
			logWriter.println("% Quality Trimmed: " + DECIMAL_FORMAT.format((double)totalQualityTrimmed / (double)totalDNAProcessed));
			logWriter.println("Number of Reads Removed: " + DECIMAL_FORMAT.format(undeterminedDNA));
			logWriter.println("% Reads Removed: " + DECIMAL_FORMAT.format((double)undeterminedDNA / (double)totalDNAProcessed));
			logWriter.println("Total Number of BP: " + DECIMAL_FORMAT.format(baseCount));
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
			if(line.length > 4){ //if there is reversed barcodes
				sampleDNAR.add(line[4].toUpperCase());
				hasReversedBarcode = true;
			}
			if(!flag){
				enzymeF = line[2];
				constEnzymesF = EnzymeList.enzymes.get(enzymeF.toUpperCase());
				if(line.length > 3){ //if there is reversed enzymes
					enzymeR = line[3];
					constEnzymesR = EnzymeList.enzymes.get(enzymeR.toUpperCase());
				}else{ //reversed enzymes will be forwards enzymes
					enzymeR = line[2];
					constEnzymesR = EnzymeList.enzymes.get(enzymeF.toUpperCase());
				}
				flag = true;
			}
		}
		
		reader.close();
	}
	
	private ArrayList<String> demultiplexFile() throws Exception{
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
		BufferedWriter[] writers2 = null; //reversed output
		BufferedWriter[] writers3 = null; //forwards index output
		BufferedWriter[] writers4 = null; //reversed index output
		BufferedWriter undeterminedWriter1 = null; //for reversed undetermined reads
		BufferedWriter undeterminedWriter2 = null; //for reversed undetermined index reads
		if(indexFile != null)
			writers3 = new BufferedWriter[sampleMapF.size() + 1];
		if(!mergePairedEnds && inputFile2 != null){
			writers2 = new BufferedWriter[sampleMapF.size()];
			if(indexFile2 != null)
				writers4 = new BufferedWriter[sampleMapF.size()];
		}
		ArrayList<String> files = new ArrayList<String>();
		long[][] DNACounts = new long[sampleMapF.size() + 1][5];
		
		String[] lines1 = new String[4];
		String[] lines2 = new String[4];
		String[] lines3 = new String[4];
		String[] lines4 = new String[4];
		while((lines1[0] = reader1.readLine()) != null && (lines1[1] = reader1.readLine()) != null &&
				(lines1[2] = reader1.readLine()) != null && (lines1[3] = reader1.readLine()) != null &&
				(inputFile2 == null || ((lines2[0] = reader2.readLine()) != null && (lines2[1] = reader2.readLine()) != null &&
				(lines2[2] = reader2.readLine()) != null && (lines2[3] = reader2.readLine()) != null)) &&
				(indexFile == null || ((lines3[0] = reader3.readLine()) != null && (lines3[1] = reader3.readLine()) != null &&
				(lines3[2] = reader3.readLine()) != null && (lines3[3] = reader3.readLine()) != null)) &&
				(indexFile2 == null || ((lines4[0] = reader4.readLine()) != null && (lines4[1] = reader4.readLine()) != null &&
				(lines4[2] = reader4.readLine()) != null && (lines4[3] = reader4.readLine()) != null))){
			totalDNAProcessed += inputFile2 == null ? 1 : 2;
			baseCount += lines1[1].length() + (inputFile2 == null ? 0 : lines2[1].length());
			boolean flag = false;
			if(totalDNAProcessed % printProcessedInterval == 0){
				logWriter.println("Reads Processed So Far: " + DECIMAL_FORMAT.format(totalDNAProcessed));
				logWriter.println("Run Time So Far: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
				logWriter.println("Current Memory Usage (GB): " + DECIMAL_FORMAT.format(UtilMethods.currentMemoryUsage()));
				logWriter.flush();
			}
			
			//save the quality and sequence lines so the undetermined reads are not trimmed or anything
			String tempF1 = lines1[1];
			String tempF2 = lines1[3];
			String tempR1 = null;
			String tempR2 = null;
			if(inputFile2 != null){
				tempR1 = lines2[1];
				tempR2 = lines2[3];
			}
			
			//check length of reads and the percentage of N
			if(minLength <= lines1[1].length() && lines1[1].length() <= maxLength && (inputFile2 == null || (minLength <= lines2[1].length() && lines2[1].length() <= maxLength)) &&
					(removeDNAWithNPercent > 1.0 || (removeDNAWithNPercent >= 0.0 && UtilMethods.percentN(lines1[1]) <= removeDNAWithNPercent && (inputFile2 == null || UtilMethods.percentN(lines2[1]) <= removeDNAWithNPercent)) ||
							(removeDNAWithNPercent < 0.0 && UtilMethods.countN(lines1[1]) <= 0 && (inputFile2 == null || UtilMethods.countN(lines2[1]) <= 0)))){
				int barcodeIndex = -1;
				int barcodeEnd = -1;
				int enzymeEnd = -1;
				int barcodeEnd2 = -1;
				int enzymeEnd2 = -1;
				int minEdit = Integer.MAX_VALUE;
				int barcodeMatchCount = 0;
				if(indexFile == null){ //check for enzyme and barcode in forwards reads
					for(int i = 0; i < sampleDNAF.size(); i++){
						ArrayList<Pair<Integer>> matches = UtilMethods.searchWithN(lines1[1].substring(0, Math.min(lines1[1].length(), maxOffset + sampleDNAF.get(i).length() + (allowIndels ? (editMaxB < 0.0 ? (int)(-editMaxB * sampleDNAF.get(i).length()) : (int)editMaxB) : 0))), sampleDNAF.get(i), editMaxB, maxOffset, allowIndels, false, minOverlapB, wildcard);
						boolean isMatch = false;
						for(int j = 0; j < matches.size(); j++){
							for(int k = 0; k < constEnzymesF.size(); k++){
								ArrayList<Pair<Integer>> enzymeMatches = UtilMethods.searchWithN(lines1[1].substring(matches.get(j).a/* + randUMILength*/, Math.min(lines1[1].length(), maxOffset + matches.get(j).a + /*randUMILength + */constEnzymesF.get(k).length() + (allowIndels ? (editMaxB < 0.0 ? (int)(-editMaxB * constEnzymesF.get(k).length()) : (int)editMaxB) : 0))), constEnzymesF.get(k), editMaxB, maxOffset, allowIndels, true, Integer.MAX_VALUE, wildcard);
								if(!enzymeMatches.isEmpty()){
									int tempBarcodeEnd2 = -1;
									int tempEnzymeEnd2 = -1;
									if(inputFile2 != null && checkReversedReads){
										int minEdit2 = Integer.MAX_VALUE;
										ArrayList<Pair<Integer>> rMatches = UtilMethods.searchWithN(lines2[1].substring(0, Math.min(lines2[1].length(), maxOffset + (hasReversedBarcode ? sampleDNAR.get(i).length() : sampleDNAF.get(i).length()) + (allowIndels ? (editMaxB < 0.0 ? (int)(-editMaxB * (hasReversedBarcode ? sampleDNAR.get(i).length() : sampleDNAF.get(i).length())) : (int)editMaxB) : 0))), hasReversedBarcode ? sampleDNAR.get(i) : UtilMethods.complement(sampleDNAF.get(i)), editMaxB, maxOffset, allowIndels, false, minOverlapB, wildcard);
										for(int ri = 0; ri < rMatches.size(); ri++){
											for(int rj = 0; rj < constEnzymesR.size(); rj++){
												ArrayList<Pair<Integer>> rEnzymeMatches = UtilMethods.searchWithN(lines2[1].substring(rMatches.get(ri).a/* + randUMILength*/, Math.min(lines2[1].length(), maxOffset + rMatches.get(ri).a + /*randUMILength + */constEnzymesR.get(rj).length() + (allowIndels ? (editMaxB < 0.0 ? (int)(-editMaxB * constEnzymesR.get(rj).length()) : (int)editMaxB) : 0))), constEnzymesR.get(rj), editMaxB, maxOffset, allowIndels, true, Integer.MAX_VALUE, wildcard);
												if(!rEnzymeMatches.isEmpty() && rMatches.get(ri).b <= minEdit2 && (rMatches.get(ri).b < minEdit2 || rMatches.get(ri).a > tempBarcodeEnd2)){
													tempBarcodeEnd2 = rMatches.get(ri).a;
													tempEnzymeEnd2 = rMatches.get(ri).a + /*randUMILength + */rEnzymeMatches.get(rEnzymeMatches.size() - 1).a;
													minEdit2 = rMatches.get(ri).b;
													break;
												}
											}
										}
									}
									
									if(inputFile2 == null || !checkReversedReads || tempEnzymeEnd2 != -1){
										isMatch = true;
										if(matches.get(j).b <= minEdit && (matches.get(j).b < minEdit || matches.get(j).a > barcodeEnd)){
											barcodeIndex = i;
											barcodeEnd = matches.get(j).a;
											enzymeEnd = matches.get(j).a + /*randUMILength + */enzymeMatches.get(enzymeMatches.size() - 1).a;
											if(inputFile2 != null){
												if(checkReversedReads){
													barcodeEnd2 = tempBarcodeEnd2;
													enzymeEnd2 = tempEnzymeEnd2;
												}else{
													barcodeEnd2 = Math.max(0, Math.min(lines2[1].length(), barcodeEnd/* - randUMILength*/));
													enzymeEnd2 = Math.max(0, Math.min(lines2[1].length(), enzymeEnd/* - randUMILength*/));
												}
											}
											minEdit = matches.get(j).b;
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
						ArrayList<Pair<Integer>> matches = UtilMethods.searchWithN(lines3[1].substring(0, Math.min(lines3[1].length(), sampleDNAF.get(i).length() + (allowIndels ? (editMaxB < 0.0 ? (int)(-editMaxB * sampleDNAF.get(i).length()) : (int)editMaxB) : 0))), sampleDNAF.get(i), editMaxB, 0, allowIndels, true, minOverlapB, wildcard);
						if(!matches.isEmpty()){
							ArrayList<Pair<Integer>> rMatches = null;
							if(inputFile2 != null && checkReversedReads){
								rMatches = UtilMethods.searchWithN(lines4[1].substring(randUMILength, Math.min(lines4[1].length(), randUMILength + (hasReversedBarcode ? sampleDNAR.get(i).length() : sampleDNAF.get(i).length()) + (allowIndels ? (editMaxB < 0.0 ? (int)(-editMaxB * (hasReversedBarcode ? sampleDNAR.get(i).length() : sampleDNAF.get(i).length())) : (int)editMaxB) : 0))),
										hasReversedBarcode ? sampleDNAR.get(i) : UtilMethods.complement(sampleDNAF.get(i)), editMaxB, 0, allowIndels, true, Integer.MAX_VALUE, wildcard);
							}
							if((inputFile2 == null || !checkReversedReads || !rMatches.isEmpty()) && matches.get(matches.size() - 1).b <= minEdit && (matches.get(matches.size() - 1).b < minEdit || matches.get(matches.size() - 1).a > barcodeEnd)){
								barcodeIndex = i;
								barcodeEnd = matches.get(matches.size() - 1).a;
								minEdit = matches.get(matches.size() - 1).b;
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
						qualityAcceptable = UtilMethods.toError(lines1[3], 0) <= qualityFilter && (inputFile2 == null || UtilMethods.toError(lines2[3], 0) <= qualityFilter);
					}else{
						qualityAcceptable = UtilMethods.toQScore(lines1[3], 0) >= qualityFilter && (inputFile2 == null || UtilMethods.toQScore(lines2[3], 0) >= qualityFilter);
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
								newSequence1 = lines1[1].substring(enzymeEnd);
								newQuality1 = lines1[3].substring(enzymeEnd);
							}else if(removeEnzyme){
								newSequence1 = lines1[1].substring(0, barcodeEnd/* + randUMILength*/) + lines1[1].substring(enzymeEnd);
								newQuality1 = lines1[3].substring(0, barcodeEnd/* + randUMILength*/) + lines1[3].substring(enzymeEnd);
							}else if(removeBarRand){
								newSequence1 = lines1[1].substring(barcodeEnd/* + randUMILength*/);
								newQuality1 = lines1[3].substring(barcodeEnd/* + randUMILength*/);
							}else{
								newSequence1 = lines1[1];
								newQuality1 = lines1[3];
							}
						}else{
							newSequence1 = lines1[1];
							newQuality1 = lines1[3];
						}
						
						if(inputFile2 != null){
							if(indexFile2 == null){
								if(removeEnzyme && removeBarRand){
									newSequence2 = lines2[1].substring(enzymeEnd2);
									newQuality2 = lines2[3].substring(enzymeEnd2);
								}else if(removeEnzyme){
									newSequence2 = lines2[1].substring(0, barcodeEnd2) + lines2[1].substring(enzymeEnd2);
									newQuality2 = lines2[3].substring(0, barcodeEnd2) + lines2[3].substring(enzymeEnd2);
								}else if(removeBarRand){
									newSequence2 = lines2[1].substring(barcodeEnd2);
									newQuality2 = lines2[3].substring(barcodeEnd2);
								}else{
									newSequence2 = lines2[1];
									newQuality2 = lines2[3];
								}
							}else{
								newSequence2 = lines2[1];
								newQuality2 = lines2[3];
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
						String[] removedAdapters = UtilMethods.removeAdapters(qualityTrimmed[0], qualityTrimmed[1], adaptersF, editMaxA, minOverlapA, allowIndels, adapterAlgorithm, wildcard);
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
							removedAdapters = UtilMethods.removeAdapters(qualityTrimmed[0], qualityTrimmed[1], adaptersR, editMaxA, minOverlapA, allowIndels, adapterAlgorithm, wildcard);
							if(qualityTrimmed[0].length() != removedAdapters[0].length()){
								removedAdapter++;
							}
							newSequence2 = removedAdapters[0];
							newQuality2 = removedAdapters[1];
						}
						
						boolean mergedReads = false;
						//merge paired-end reads
						if(mergePairedEnds && inputFile2 != null){
							String[] merged = UtilMethods.mergeReads(newSequence1, newQuality1, newSequence2, newQuality2, editMaxM, mergeAlgorithm, wildcard);
							if(merged[0].length() != newSequence1.length() + newSequence2.length()){
								mergedReads = true;
							}
							newSequence1 = merged[0];
							newQuality1 = merged[1];
						}
						
						if((!removeUntrimmedReads || trimmedQuality == (inputFile2 == null ? 1 : 2)) && (!removeNoAdapterReads || removedAdapter == (inputFile2 == null ? 1 : 2)) &&
								(!removeNoMergeReads || mergedReads)){
							if(writers1[barcodeIndex] == null){
								if(!outputGZIP && !removeFirstDup && !removeBestDup){
									writers1[barcodeIndex] = new BufferedWriter(new FileWriter(outputDir + "sample_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R1.fastq"), BUFFER_SIZE);
									files.add(outputDir + "sample_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R1.fastq");
									if(indexFile != null){
										writers3[barcodeIndex] = new BufferedWriter(new FileWriter(outputDir + "index_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R1.fastq"), BUFFER_SIZE);
										files.add(outputDir + "index_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R1.fastq");
									}
								}else{
									writers1[barcodeIndex] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(((removeFirstDup || removeBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "sample_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R1.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
									files.add(((removeFirstDup || removeBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "sample_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R1.fastq.gz");
									if(indexFile != null){
										writers3[barcodeIndex] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(((removeFirstDup || removeBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "index_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R1.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
										files.add(((removeFirstDup || removeBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "index_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R1.fastq.gz");
									}
								}
							}
							if(!mergePairedEnds && inputFile2 != null && writers2[barcodeIndex] == null){
								if(!outputGZIP && !removeFirstDup && !removeBestDup){
									writers2[barcodeIndex] = new BufferedWriter(new FileWriter(outputDir + "sample_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R2.fastq"), BUFFER_SIZE);
									files.add(outputDir + "sample_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R2.fastq");
									if(indexFile2 != null){
										writers4[barcodeIndex] = new BufferedWriter(new FileWriter(outputDir + "index_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R2.fastq"), BUFFER_SIZE);
										files.add(outputDir + "index_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R2.fastq");
									}
								}else{
									writers2[barcodeIndex] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(((removeFirstDup || removeBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "sample_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R2.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
									files.add(((removeFirstDup || removeBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "sample_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R2.fastq.gz");
									if(indexFile2 != null){
										writers4[barcodeIndex] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(((removeFirstDup || removeBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "index_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R2.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
										files.add(((removeFirstDup || removeBestDup) ? (outputDir + "temp" + File.separatorChar) : outputDir) + "index_" + sampleMapF.get(sampleDNAF.get(barcodeIndex)) + "_R2.fastq.gz");
									}
								}
							}
							
							DNACounts[barcodeIndex][0] += inputFile2 == null ? 1 : 2;
							DNACounts[barcodeIndex][1] += lines1[1].length() + (inputFile2 == null ? 0 : lines2[1].length());
							DNACounts[barcodeIndex][3] += removedAdapter;
							totalRemovedAdapters += removedAdapter;
							DNACounts[barcodeIndex][2] += mergedReads ? 2 : 0;
							totalReadsMerged += mergedReads ? 2 : 0;
							DNACounts[barcodeIndex][4] += trimmedQuality;
							totalQualityTrimmed += trimmedQuality;
							
							//print to the whatever file the read belongs to
							writers1[barcodeIndex].write(lines1[0]);
							writers1[barcodeIndex].newLine();
							writers1[barcodeIndex].write(newSequence1);
							writers1[barcodeIndex].newLine();
							writers1[barcodeIndex].write(description2);
							writers1[barcodeIndex].newLine();
							writers1[barcodeIndex].write(newQuality1);
							writers1[barcodeIndex].newLine();
							
							if(indexFile != null){
								writers3[barcodeIndex].write(lines3[0]);
								writers3[barcodeIndex].newLine();
								writers3[barcodeIndex].write(lines3[1]);
								writers3[barcodeIndex].newLine();
								writers3[barcodeIndex].write(description2);
								writers3[barcodeIndex].newLine();
								writers3[barcodeIndex].write(lines3[3]);
								writers3[barcodeIndex].newLine();
							}
							
							if(!mergePairedEnds && inputFile2 != null){
								writers2[barcodeIndex].write(lines2[0]);
								writers2[barcodeIndex].newLine();
								writers2[barcodeIndex].write(newSequence2);
								writers2[barcodeIndex].newLine();
								writers2[barcodeIndex].write(description2);
								writers2[barcodeIndex].newLine();
								writers2[barcodeIndex].write(newQuality2);
								writers2[barcodeIndex].newLine();
								
								if(indexFile2 != null){
									writers4[barcodeIndex].write(lines4[0]);
									writers4[barcodeIndex].newLine();
									writers4[barcodeIndex].write(lines4[1]);
									writers4[barcodeIndex].newLine();
									writers4[barcodeIndex].write(description2);
									writers4[barcodeIndex].newLine();
									writers4[barcodeIndex].write(lines4[3]);
									writers4[barcodeIndex].newLine();
								}
							}
							flag = true;
						}
					}
				}
			}
			if(!flag){
				if(writers1[writers1.length - 1] == null){ //initialize writer for the undetermined file if it has not already been initialized
					if(outputGZIP){
						writers1[writers1.length - 1] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "sample_undetermined_R1.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
						if(inputFile2 != null)
							undeterminedWriter1 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "sample_undetermined_R2.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
						if(indexFile != null){
							writers3[writers3.length - 1] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "index_undetermined_R1.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
							if(inputFile2 != null)
								undeterminedWriter2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputDir + "index_undetermined_R2.fastq.gz"), BUFFER_SIZE_GZIP)), BUFFER_SIZE);
						}
					}else{
						writers1[writers1.length - 1] = new BufferedWriter(new FileWriter(outputDir + "sample_undetermined_R1.fastq"), BUFFER_SIZE);
						if(inputFile2 != null)
							undeterminedWriter1 = new BufferedWriter(new FileWriter(outputDir + "sample_undetermined_R2.fastq"), BUFFER_SIZE);
						if(indexFile != null){
							writers3[writers3.length - 1] = new BufferedWriter(new FileWriter(outputDir + "index_undetermined_R1.fastq"), BUFFER_SIZE);
							if(inputFile2 != null)
								undeterminedWriter2 = new BufferedWriter(new FileWriter(outputDir + "index_undetermined_R2.fastq"), BUFFER_SIZE);
						}
					}
					generatedUndeterminedFile = true;
				}
				DNACounts[DNACounts.length - 1][0] += inputFile2 == null ? 1 : 2;
				DNACounts[DNACounts.length - 1][1] += lines1[1].length() + (inputFile2 == null ? 0 : lines2[1].length());
				//print to undetermined file
				writers1[writers1.length - 1].write(lines1[0]);
				writers1[writers1.length - 1].newLine();
				writers1[writers1.length - 1].write(tempF1);
				writers1[writers1.length - 1].newLine();
				writers1[writers1.length - 1].write(description2);
				writers1[writers1.length - 1].newLine();
				writers1[writers1.length - 1].write(tempF2);
				writers1[writers1.length - 1].newLine();
				if(indexFile != null){
					writers3[writers3.length - 1].write(lines3[0]);
					writers3[writers3.length - 1].newLine();
					writers3[writers3.length - 1].write(lines3[1]);
					writers3[writers3.length - 1].newLine();
					writers3[writers3.length - 1].write(description2);
					writers3[writers3.length - 1].newLine();
					writers3[writers3.length - 1].write(lines3[3]);
					writers3[writers3.length - 1].newLine();
				}
				if(inputFile2 != null){
					undeterminedWriter1.write(lines2[0]);
					undeterminedWriter1.newLine();
					undeterminedWriter1.write(tempR1);
					undeterminedWriter1.newLine();
					undeterminedWriter1.write(description2);
					undeterminedWriter1.newLine();
					undeterminedWriter1.write(tempR2);
					undeterminedWriter1.newLine();
					if(indexFile2 != null){
						undeterminedWriter2.write(lines4[0]);
						undeterminedWriter2.newLine();
						undeterminedWriter2.write(lines4[1]);
						undeterminedWriter2.newLine();
						undeterminedWriter2.write(description2);
						undeterminedWriter2.newLine();
						undeterminedWriter2.write(lines4[3]);
						undeterminedWriter2.newLine();
					}
				}
				undeterminedDNA += inputFile2 == null ? 1 : 2;
			}
		}
		
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
			if(undeterminedWriter1 != null)
				undeterminedWriter1.close();
			if(indexFile2 != null){
				reader4.close();
				if(undeterminedWriter2 != null)
					undeterminedWriter2.close();
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
			map.put(Stat.SEQUENCE_COUNT, DECIMAL_FORMAT.format(DNACounts[i][0]));
			map.put(Stat.SEQUENCE_PERCENT, DECIMAL_FORMAT.format((double)DNACounts[i][0] / (double)totalDNAProcessed));
			map.put(Stat.MERGED_COUNT, DECIMAL_FORMAT.format(DNACounts[i][2]));
			map.put(Stat.MERGED_PERCENT, DECIMAL_FORMAT.format((double)DNACounts[i][2] / (double)totalDNAProcessed));
			map.put(Stat.REMOVEADAPTER_COUNT, DECIMAL_FORMAT.format(DNACounts[i][3]));
			map.put(Stat.REMOVEADAPTER_PERCENT, DECIMAL_FORMAT.format((double)DNACounts[i][3] / (double)DNACounts[i][0]));
			map.put(Stat.QUALITYTRIM_COUNT, DECIMAL_FORMAT.format(DNACounts[i][4]));
			map.put(Stat.QUALITYTRIM_PERCENT, DECIMAL_FORMAT.format((double)DNACounts[i][4] / (double)DNACounts[i][0]));
			map.put(Stat.BASE_COUNT, DECIMAL_FORMAT.format(DNACounts[i][1]));
			map.put(Stat.BASE_PERCENT, DECIMAL_FORMAT.format((double)DNACounts[i][1] / (double)baseCount));
		}
		stats.get("Undetermined").put(Stat.SEQUENCE_COUNT, DECIMAL_FORMAT.format(DNACounts[DNACounts.length - 1][0]));
		stats.get("Undetermined").put(Stat.SEQUENCE_PERCENT, DECIMAL_FORMAT.format((double)DNACounts[DNACounts.length - 1][0] / (double)totalDNAProcessed));
		stats.get("Undetermined").put(Stat.BASE_COUNT, DECIMAL_FORMAT.format(DNACounts[DNACounts.length - 1][1]));
		stats.get("Undetermined").put(Stat.BASE_PERCENT, DECIMAL_FORMAT.format((double)DNACounts[DNACounts.length - 1][1] / (double)baseCount));
		
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
					totalDNAProcessed += readPath2 == null ? 1 : 2;
					baseCount += lines1[1].length() + (readPath2 == null ? 0 : lines2[1].length());
					if(totalDNAProcessed % printProcessedInterval == 0){
						logWriter.println("Reads Processed So Far: " + DECIMAL_FORMAT.format(totalDNAProcessed));
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
			
			HashMap<BitSet, Read> map = new HashMap<BitSet, Read>();
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
					totalDNAProcessed += readPath2 == null ? 1 : 2;
					baseCount += lines1[1].length() + (readPath2 == null ? 0 : lines2[1].length());
					if(totalDNAProcessed % printProcessedInterval == 0){
						logWriter.println("Reads Processed So Far: " + DECIMAL_FORMAT.format(totalDNAProcessed));
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
					Read other = map.get(key);
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
							map.put(key, new Read(currQuality, lines1[0], lines1[1], lines1[3], lines3[0], lines3[1], lines3[3]));
						else
							map.put(key, new Read(currQuality, lines1[0], lines1[1], lines1[3], lines3[0], lines3[1], lines3[3], lines2[0], lines2[1], lines2[3], lines4[0], lines4[1], lines4[3]));
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
						map.put(key, new Read(currQuality, lines1[0], lines1[1], lines1[3], lines3[0], lines3[1], lines3[3]));
					else
						map.put(key, new Read(currQuality, lines1[0], lines1[1], lines1[3], lines3[0], lines3[1], lines3[3], lines2[0], lines2[1], lines2[3], lines4[0], lines4[1], lines4[3]));
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
			for(Read val : map.values()){
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
		
		String[] lines1 = new String[4];
		String[] lines2 = new String[4];
		while((lines1[0] = reader1.readLine()) != null && (lines1[1] = reader1.readLine()) != null &&
				(lines1[2] = reader1.readLine()) != null && (lines1[3] = reader1.readLine()) != null &&
				(lines2[0] = reader2.readLine()) != null && (lines2[1] = reader2.readLine()) != null &&
				(lines2[2] = reader2.readLine()) != null && (lines2[3] = reader2.readLine()) != null){
			totalDNAProcessed += 2;
			baseCount += lines1[1].length() + lines2[1].length();
			if(totalDNAProcessed % printProcessedInterval == 0){
				logWriter.println("Reads Processed So Far: " + DECIMAL_FORMAT.format(totalDNAProcessed));
				logWriter.println("Run Time So Far: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
				logWriter.println("Current Memory Usage (GB): " + DECIMAL_FORMAT.format(UtilMethods.currentMemoryUsage()));
				logWriter.flush();
			}
			//merge the two lines
			String[] merged = UtilMethods.mergeReads(lines1[1], lines1[3], lines2[1], lines2[3], editMaxM, mergeAlgorithm, wildcard);
			if(!removeNoMergeReads || merged[0].length() != lines1[1].length() + lines2[1].length()){
				totalReadsMerged += 2;
				writer.write(lines1[0]);
				writer.newLine();
				writer.write(merged[0]);
				writer.newLine();
				writer.write(description2);
				writer.newLine();
				writer.write(merged[1]);
				writer.newLine();
			}else{
				undeterminedDNA += 2;
			}
		}
		
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
		
		String[] lines = new String[4];
		while((lines[0] = reader.readLine()) != null && (lines[1] = reader.readLine()) != null &&
				(lines[2] = reader.readLine()) != null && (lines[3] = reader.readLine()) != null){
			totalDNAProcessed++;
			baseCount += lines[1].length();
			if(totalDNAProcessed % printProcessedInterval == 0){
				logWriter.println("Reads Processed So Far: " + DECIMAL_FORMAT.format(totalDNAProcessed));
				logWriter.println("Run Time So Far: " + UtilMethods.formatElapsedTime(System.currentTimeMillis() - startTime));
				logWriter.println("Current Memory Usage (GB): " + DECIMAL_FORMAT.format(UtilMethods.currentMemoryUsage()));
				logWriter.flush();
			}
			
			//check if quality is good enough
			boolean qualityAcceptable;
			if(filterAlgorithm){
				qualityAcceptable = UtilMethods.toError(lines[3], 0) <= qualityFilter;
			}else{
				qualityAcceptable = UtilMethods.toQScore(lines[3], 0) >= qualityFilter;
			}
			if(lines[1].length() >= minLength && lines[1].length() <= maxLength && qualityAcceptable &&
					(removeDNAWithNPercent > 1.0 || (removeDNAWithNPercent >= 0.0 && UtilMethods.percentN(lines[1]) <= removeDNAWithNPercent) || (removeDNAWithNPercent < 0.0 && UtilMethods.countN(lines[1]) <= 0))){
				String[] temp;
				String[] qualityTrimmed;
				
				//trim 'N' bp
				temp = UtilMethods.trimN(lines[1], lines[3], trimNPercent);
				lines[1] = temp[0];
				lines[3] = temp[1];
				
				int trimmedQuality = 0;
				int removedAdapter = 0;
				
				//quality trim
				if(trimAlgorithm){
					temp = UtilMethods.qualityTrim2(lines[1], lines[3], qualityTrimQScore1, true, qualityTrimLength);
					qualityTrimmed = UtilMethods.qualityTrim2(temp[0], temp[1], qualityTrimQScore2, false, qualityTrimLength);
				}else{
					temp = UtilMethods.qualityTrim1(lines[1], lines[3], qualityTrimQScore1, true, qualityTrimLength);
					qualityTrimmed = UtilMethods.qualityTrim1(temp[0], temp[1], qualityTrimQScore2, false, qualityTrimLength);
				}
				if(lines[1].length() != qualityTrimmed[0].length()){
					trimmedQuality++;
				}
				//remove adapters
				String[] removedAdapters = UtilMethods.removeAdapters(qualityTrimmed[0], qualityTrimmed[1], adaptersF, editMaxA, minOverlapA, allowIndels, adapterAlgorithm, wildcard);
				if(qualityTrimmed[0].length() != removedAdapters[0].length()){
					removedAdapter++;
				}
				lines[1] = removedAdapters[0];
				lines[3] = removedAdapters[1];
				
				if((!removeUntrimmedReads || trimmedQuality == 1) && (!removeNoAdapterReads || removedAdapter == 1)){
					totalQualityTrimmed += trimmedQuality;
					totalRemovedAdapters += removedAdapter;
					
					//print to file
					writer.write(lines[0]);
					writer.newLine();
					writer.write(lines[1]);
					writer.newLine();
					writer.write(description2);
					writer.newLine();
					writer.write(lines[3]);
					writer.newLine();
				}else{
					undeterminedDNA++;
				}
			}else{
				undeterminedDNA++;
			}
		}
		
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
				if(args[i].equals("--keepfirst") || args[i].equals("-k")){
					removeFirstDup = true;
					removeBestDup = false;
				}else if(args[i].equals("--keepbest") || args[i].equals("-K")){
					removeBestDup = true;
					removeFirstDup = false;
				}else if(args[i].equals("--umi")){
					randUMILength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--non")){
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						removeDNAWithNPercent = Double.parseDouble(args[++i]);
					}else{
						removeDNAWithNPercent = -1.0;
					}
				}else if(args[i].equals("--maxoffset")){
					maxOffset = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--qfilter")){
					qualityFilter = Double.parseDouble(args[++i]);
				}else if(args[i].equals("--reads") || args[i].equals("-r")){
					inputFile = new File(args[++i]);
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						inputFile2 = new File(args[++i]);
					}
				}else if(args[i].equals("--index") || args[i].equals("-i")){
					indexFile = new File(args[++i]);
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						indexFile2 = new File(args[++i]);
					}
					randUMILength = 12;
				}else if(args[i].equals("--sample") || args[i].equals("-s")){
					sampleFile = new File(args[++i]);
				}else if(args[i].equals("--output") || args[i].equals("-o")){
					outputDir = args[++i];
					if(!outputDir.endsWith(File.separator) && !outputDir.endsWith("/") && !outputDir.endsWith("\\")){
						outputDir += File.separator;
					}
				}else if(args[i].equals("--gzip") || args[i].equals("-gz")){
					outputGZIP = true;
				}else if(args[i].equals("--cleardir")){
					clearDir = true;
				}else if(args[i].equals("--savetemp")){
					saveTemp = true;
				}else if(args[i].equals("--keepbarcode")){
					removeBarRand = false;
				}else if(args[i].equals("--keepenzyme")){
					removeEnzyme = false;
				}else if(args[i].equals("--maxedita")){
					editMaxA = Double.parseDouble(args[++i]);
				}else if(args[i].equals("--maxeditb")){
					editMaxB = Double.parseDouble(args[++i]);
				}else if(args[i].equals("--maxeditm")){
					editMaxM = Double.parseDouble(args[++i]);
				}else if(args[i].equals("--maxeditap")){
					editMaxA = -Double.parseDouble(args[++i]);
				}else if(args[i].equals("--maxeditbp")){
					editMaxB = -Double.parseDouble(args[++i]);
				}else if(args[i].equals("--maxeditmp")){
					editMaxM = -Double.parseDouble(args[++i]);
				}else if(args[i].equals("--savedup")){
					saveDup = true;
				}else if(args[i].equals("--printprocessed")){
					printProcessedInterval = Long.parseLong(args[++i]);
				}else if(args[i].equals("--printduplicate")){
					printDuplicateInterval = Long.parseLong(args[++i]);
				}else if(args[i].equals("--pairmerge") || args[i].equals("-m")){
					mergePairedEnds = true;
				}else if(args[i].equals("--minlen")){
					minLength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--maxlen")){
					maxLength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--fadapterstart") || args[i].equals("-a")){
					while(i + 1 < args.length && !args[i + 1].startsWith("-")){
						if(args[i + 1].startsWith("^")){
							adaptersF.add(new Adapter(args[i + 1].substring(1), true, true));
						}else{
							adaptersF.add(new Adapter(args[i + 1], true, false));
						}
						i++;
					}
				}else if(args[i].equals("--fadapterend") || args[i].equals("-A")){
					while(i + 1 < args.length && !args[i + 1].startsWith("-")){
						if(args[i + 1].endsWith("$")){
							adaptersF.add(new Adapter(args[i + 1].substring(0, args[i + 1].length() - 1), false, true));
						}else{
							adaptersF.add(new Adapter(args[i + 1], false, false));
						}
						i++;
					}
				}else if(args[i].equals("--radapterstart") || args[i].equals("-z")){
					while(i + 1 < args.length && !args[i + 1].startsWith("-")){
						if(args[i + 1].startsWith("^")){
							adaptersR.add(new Adapter(args[i + 1].substring(1), true, true));
						}else{
							adaptersR.add(new Adapter(args[i + 1], true, false));
						}
						i++;
					}
				}else if(args[i].equals("--radapterend") || args[i].equals("-Z")){
					while(i + 1 < args.length && !args[i + 1].startsWith("-")){
						if(args[i + 1].endsWith("$")){
							adaptersR.add(new Adapter(args[i + 1].substring(0, args[i + 1].length() - 1), false, true));
						}else{
							adaptersR.add(new Adapter(args[i + 1], false, false));
						}
						i++;
					}
				}else if(args[i].equals("--minoverlapa")){
					minOverlapA = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--minoverlapb")){
					minOverlapB = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--altqtrim")){
					trimAlgorithm = true;
				}else if(args[i].equals("--qtrim") || args[i].equals("-q")){
					if(i + 2 >= args.length || args[i + 2].startsWith("-")){
						qualityTrimQScore2 = Integer.parseInt(args[++i]);
					}else{
						qualityTrimQScore1 = Integer.parseInt(args[++i]);
						qualityTrimQScore2 = Integer.parseInt(args[++i]);
					}
				}else if(args[i].equals("--qtrimlen")){
					qualityTrimLength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--altmerge")){
					mergeAlgorithm = true;
				}else if(args[i].equals("--ntrim")){
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						trimNPercent = Double.parseDouble(args[++i]);
					}else{
						trimNPercent = 0.5;
					}
				}else if(args[i].equals("--indels")){
					allowIndels = true;
				}else if(args[i].equals("--altadapter")){
					adapterAlgorithm = true;
				}else if(args[i].equals("--altqfilter")){
					filterAlgorithm = true;
				}else if(args[i].equals("--checkreversed")){
					checkReversedReads = true;
				}else if(args[i].equals("--nowildcard")){
					wildcard = false;
				}else if(args[i].equals("--singlematch")){
					singleBarcodeMatchOnly = true;
				}else if(args[i].equals("--removeuntrimmed")){
					removeUntrimmedReads = true;
				}else if(args[i].equals("--removenoadapter")){
					removeNoAdapterReads = true;
				}else if(args[i].equals("--removenomerge")){
					removeNoMergeReads = true;
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
				if(args[i].equals("--keepfirst") || args[i].equals("-k")){
					removeFirstDup = true;
					removeBestDup = false;
				}else if(args[i].equals("--keepbest") || args[i].equals("-K")){
					removeBestDup = true;
					removeFirstDup = false;
				}else if(args[i].equals("--reads") || args[i].equals("-r")){
					inputFile = new File(args[++i]);
				}else if(args[i].equals("--index") || args[i].equals("-i")){
					indexFile = new File(args[++i]);
					randUMILength = 12;
				}else if(args[i].equals("--output") || args[i].equals("-o")){
					outputDir = args[++i];
					if(!outputDir.endsWith(File.separator) && !outputDir.endsWith("/") && !outputDir.endsWith("\\")){
						outputDir += File.separator;
					}
				}else if(args[i].equals("--overwrite") || args[i].equals("-O")){
					replaceOriginal = true;
				}else if(args[i].equals("--gzip") || args[i].equals("-gz")){
					outputGZIP = true;
				}else if(args[i].equals("--savedup")){
					saveDup = true;
				}else if(args[i].equals("--umi")){
					randUMILength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--printprocessed")){
					printProcessedInterval = Long.parseLong(args[++i]);
				}else if(args[i].equals("--printduplicate")){
					printDuplicateInterval = Long.parseLong(args[++i]);
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
				if(args[i].equals("--reads") || args[i].equals("-r")){
					inputFile = new File(args[++i]);
					inputFile2 = new File(args[++i]);
				}else if(args[i].equals("--output") || args[i].equals("-o")){
					outputDir = args[++i];
					if(!outputDir.endsWith(File.separator) && !outputDir.endsWith("/") && !outputDir.endsWith("\\")){
						outputDir += File.separator;
					}
				}else if(args[i].equals("--gzip") || args[i].equals("-gz")){
					outputGZIP = true;
				}else if(args[i].equals("--maxeditm")){
					editMaxM = Double.parseDouble(args[++i]);
				}else if(args[i].equals("--maxeditmp")){
					editMaxM = -Double.parseDouble(args[++i]);
				}else if(args[i].equals("--printprocessed")){
					printProcessedInterval = Long.parseLong(args[++i]);
				}else if(args[i].equals("--altmerge")){
					mergeAlgorithm = true;
				}else if(args[i].equals("--nowildcard")){
					wildcard = false;
				}else if(args[i].equals("--removenomerge")){
					removeNoMergeReads = true;
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
				if(args[i].equals("--reads") || args[i].equals("-r")){
					inputFile = new File(args[++i]);
				}else if(args[i].equals("--output") || args[i].equals("-o")){
					outputDir = args[++i];
					if(!outputDir.endsWith(File.separator) && !outputDir.endsWith("/") && !outputDir.endsWith("\\")){
						outputDir += File.separator;
					}
				}else if(args[i].equals("--gzip") || args[i].equals("-gz")){
					outputGZIP = true;
				}else if(args[i].equals("--maxedita")){
					editMaxA = Double.parseDouble(args[++i]);
				}else if(args[i].equals("--maxeditap")){
					editMaxA = -Double.parseDouble(args[++i]);
				}else if(args[i].equals("--printprocessed")){
					printProcessedInterval = Long.parseLong(args[++i]);
				}else if(args[i].equals("--qfilter")){
					qualityFilter = Double.parseDouble(args[++i]);
				}else if(args[i].equals("--overwrite") || args[i].equals("-O")){
					replaceOriginal = true;
				}else if(args[i].equals("--minlen")){
					minLength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--maxlen")){
					maxLength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--fadapterstart") || args[i].equals("-a")){
					while(i + 1 < args.length && !args[i + 1].startsWith("-")){
						if(args[i + 1].startsWith("^")){
							adaptersF.add(new Adapter(args[i + 1].substring(1), true, true));
						}else{
							adaptersF.add(new Adapter(args[i + 1], true, false));
						}
						i++;
					}
				}else if(args[i].equals("--fadapterend") || args[i].equals("-A")){
					while(i + 1 < args.length && !args[i + 1].startsWith("-")){
						if(args[i + 1].endsWith("$")){
							adaptersF.add(new Adapter(args[i + 1].substring(0, args[i + 1].length() - 1), false, true));
						}else{
							adaptersF.add(new Adapter(args[i + 1], false, false));
						}
						i++;
					}
				}else if(args[i].equals("--minoverlapa")){
					minOverlapA = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--altqtrim")){
					trimAlgorithm = true;
				}else if(args[i].equals("--qtrim") || args[i].equals("-q")){
					if(i + 2 >= args.length || args[i + 2].startsWith("-")){
						qualityTrimQScore2 = Integer.parseInt(args[++i]);
					}else{
						qualityTrimQScore1 = Integer.parseInt(args[++i]);
						qualityTrimQScore2 = Integer.parseInt(args[++i]);
					}
				}else if(args[i].equals("--qtrimlen")){
					qualityTrimLength = Integer.parseInt(args[++i]);
				}else if(args[i].equals("--non")){
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						removeDNAWithNPercent = Double.parseDouble(args[++i]);
					}else{
						removeDNAWithNPercent = -1.0;
					}
				}else if(args[i].equals("--ntrim")){
					if(i + 1 < args.length && !args[i + 1].startsWith("-")){
						trimNPercent = Double.parseDouble(args[++i]);
					}else{
						trimNPercent = 0.5;
					}
				}else if(args[i].equals("--indels")){
					allowIndels = true;
				}else if(args[i].equals("--altadapter")){
					adapterAlgorithm = true;
				}else if(args[i].equals("--altqfilter")){
					filterAlgorithm = true;
				}else if(args[i].equals("--nowildcard")){
					wildcard = false;
				}else if(args[i].equals("--removeuntrimmed")){
					removeUntrimmedReads = true;
				}else if(args[i].equals("--removenoadapter")){
					removeNoAdapterReads = true;
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
		}
		if(logWriter != null){
			logWriter.close();
		}
	}
}
